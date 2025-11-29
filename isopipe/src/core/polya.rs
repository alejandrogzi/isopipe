use memchr::memchr;
use std::{
    collections::HashMap,
    fs::{create_dir_all, File, OpenOptions},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use crate::{config::*, consts::*, executor::job::Job, isotools, rm};

/// Run polya mod [3 steps]
///
/// # Arguments
///
/// * `step` - The pipeline step to run
/// * `config` - The configuration to use
/// * `input_dir` - The input directory
/// * `output_dir` - The output directory
///
/// # Returns
///
/// A vector of jobs to run
///
/// # Examples
///
/// ```rust, no_run
/// use isopipe::core::polya;
/// use isopipe::config::Config;
/// use std::path::PathBuf;
///
/// let config = Config::default();
/// let input_dir = PathBuf::from("/path/to/input");
/// let output_dir = PathBuf::from("/path/to/output");
///
/// let jobs = polya(&step, &config, &input_dir, &output_dir);
/// ```
pub fn polya(
    step: &PipelineStep,
    config: &Config,
    input_dir: &Path,
    step_output_dir: &Path,
) -> Vec<Job> {
    let binary = isotools!(ISO_POLYA);
    let mut jobs = Vec::new();

    let parts = &step_output_dir.join(POLYA_PARTS);
    let keep_temp = config
        .get_step_custom_field(step, KEEP_TEMP)
        .parse::<bool>()
        .unwrap_or(false);

    let _ = create_dir_all(parts);

    if !keep_temp {
        log::info!(
            "INFO [STEP 6]: Removing minimap2 chunked directory: {}",
            input_dir.join(CHUNKS).display()
        );
        rm!(input_dir.join(CHUNKS));
    } else {
        log::info!("INFO [STEP 6]: Keeping temporary files -> keep_temp set to true!",);
    }

    let args = config.get_step_args(
        step,
        vec![INPUT_DIR, OUTPUT_DIR, MEMORY, TIME, TOGA, KEEP_TEMP],
    );

    // WARN: input_dir needs to be suffixed by /bam
    for (batch, entry) in std::fs::read_dir(input_dir.join(BAM))
        .expect("Failed to read assets directory")
        .flatten()
        .filter(|entry| {
            entry
                .path()
                .file_name()
                .and_then(|ext| ext.to_str())
                .map(|name| name.ends_with(BAM))
                .unwrap_or(false)
        })
        .enumerate()
    {
        let bam = entry.path();
        log::debug!(
            "DEBUG [STEP 6]: processing {} with --batch {batch}...",
            bam.display()
        );

        // INFO: for cannonical runs this yields -> chunk*{run}.singletons
        // INFO: and outputs chunk*{run}.singletons.good.bed
        let bind = bam.with_extension("");
        let prefix = bind
            .file_name()
            .unwrap_or_else(|| panic!("ERROR: could not build prefix for {:?}", bam));

        if std::fs::metadata(&bam).unwrap().len() == 0 {
            log::warn!("WARN: {} its empty!", bam.display());
            continue;
        }

        let mut cmd = format!(
            "{} {} --bam {} {args} --prefix {} --outdir {} --batch {}",
            binary.display(),
            SEGMENT,
            bam.display(),
            prefix.to_string_lossy(),
            &parts.display(),
            batch + 1,
        );

        if bam
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("")
            .ends_with("singletons.bam")
        {
            cmd = format!("{} --singleton", cmd,); // INFO: SG tag
        }

        if !keep_temp {
            cmd = format!("{} && rm {} {}.bai", &cmd, bam.display(), bam.display());
        }

        log::debug!("DEBUG [STEP 6]: running command: {}", cmd);
        jobs.push(Job::from(cmd));
    }

    log::info!("INFO [STEP 6]: Pre-processing completed -> Running...");

    jobs
}

/// Merges polya .bed results into a single .bed per category.
/// Where categories are: singletons, accepts and rejections.
/// After, it chunks each category into per-chromosome files.
///
/// # Arguments
///
/// * `input_dir` - The input directory
///
/// # Example
///
/// ```rust, no_run
/// let input_dir = PathBuf::from("/path/to/input");
/// merge(input_dir);
/// ```
pub fn merge(input_dir: &PathBuf, config: &Config, step: &PipelineStep) {
    let parts = &input_dir.join(POLYA_PARTS);

    log::info!(
        "INFO [MERGE]: Merging polyA segmentation results in {}...",
        &parts.display()
    );

    let keep_temp = config
        .get_step_custom_field(step, KEEP_TEMP)
        .parse::<bool>()
        .unwrap_or(false);

    if keep_temp {
        log::info!("INFO [MERGE]: Keeping temporary files -> keep_temp set to true!",);
    }

    if input_dir.join(POLYA_PARTS).exists() {
        log::info!("INFO [MERGE]: chunked directory found, merging parts...");

        let files = scan_dir(&input_dir.join(POLYA_PARTS), BED).unwrap_or_else(|| {
            log::error!(
                "ERROR: could not find any files with extension {BED:?} in {}!",
                input_dir.join(POLYA_PARTS).display()
            );
            std::process::exit(1);
        });

        // INFO: partition paths into their respective categories
        let (singletons, accepts, rejections): (Vec<_>, Vec<_>, Vec<_>) = files.into_iter().fold(
            (Vec::new(), Vec::new(), Vec::new()),
            |(mut _s, mut g, mut b), path| {
                match &path.file_name().unwrap().to_string_lossy() {
                    p if p.ends_with(BED_SGN_ACCEPT) => g.push(path), // WARN: since 3/10/25 -> sending singletons to good [they are already with SG]
                    p if p.ends_with(BED_ACCEPT) => g.push(path),
                    p if p.ends_with(BED_REJECT) => b.push(path),
                    _ => {}
                }
                (_s, g, b)
            },
        );

        log::debug!("DEBUG: collected all .bed types -> SGN_ACCEPT: {singletons:?}, BED_ACCEPT: {accepts:?}, BED_REJECT: {rejections:?}");

        for (file, group) in ALN_POLYA_FILES
            .iter()
            .zip([singletons, accepts, rejections].iter())
        {
            // INFO: does not make any sense cat rejected files
            if *file == ALN_POLYA_REJECT && !keep_temp {
                log::info!("INFO [MERGE]: skipping rejected files from merging step...");
                continue;
            }

            if group.is_empty() {
                log::warn!(
                    "WARN [MERGE]: No {} files found in {}",
                    file,
                    input_dir.join(POLYA_PARTS).display()
                );
                continue;
            }

            log::info!(
                "INFO [MERGE]: Trying to merge {} files to per chromsome {} file",
                group.len(),
                file
            );

            let _ = __split_by_chr(group, input_dir, file);
        }

        if !keep_temp {
            rm!(input_dir.join(POLYA_PARTS));
        }
    } else {
        log::warn!(
            "WARN [MERGE]: parts/ directory not found under {}, skipping merge step and grabbing .bed files directly from input directory!",
            input_dir.display()
        );

        let files = scan_dir(input_dir, BED).unwrap_or_else(|| {
            log::error!(
                "ERROR: could not find any files with extension {BED:?} in {}!",
                input_dir.display()
            );
            std::process::exit(1);
        });
        log::debug!("DEBUG: collected all .bed types -> {files:?}");

        log::warn!(
            "WARN [MERGE]: forcing chunking on non-cannonical files or already merged files in {}!",
            input_dir.display()
        );

        for bed in files {
            let suffix = bed
                .file_name()
                .unwrap_or_else(|| panic!("ERROR: could not build suffix for {:?}", bed))
                .to_string_lossy();

            let sp = __split_by_chr(&[bed.clone()], input_dir, &suffix);

            if let Err(e) = sp {
                log::error!(
                    "ERROR [MERGE]: failed to split {} into chromosomes: {}",
                    bed.display(),
                    e
                );
            } else {
                log::info!(
                    "INFO [MERGE]: successfully split {} into chunked chromosomes!",
                    bed.display()
                );

                if !keep_temp {
                    rm!(bed);
                }
            }
        }
    }
}

/// Splits one or more input files into multiple output files, where each output
/// file contains records pertaining to a specific chromosome.
///
/// This function reads lines from each input file. For each line, it identifies
/// the chromosome name (assumed to be the first tab-separated field, or the
/// entire line if no tab is present). It then writes the line to a dedicated
/// output file named after the chromosome and a specified suffix, located within
/// the given output directory. If an output file for a chromosome does not exist,
/// it is created; otherwise, lines are appended to the existing file.
///
/// # Arguments
///
/// * `files` - A slice of `PathBuf` representing the paths to the input files to be split.
/// * `dir` - A `PathBuf` representing the directory where the chromosome-specific output files will be created.
/// * `suffix` - A string slice that will be appended to the chromosome name to form the output file names (e.g., `chr1_suffix`).
///
/// # Returns
///
/// A `std::io::Result<()>` which is:
/// - `Ok(())`: If all files are successfully processed and written.
/// - `Err(std::io::Error)`: If any I/O error occurs during file operations
///   (e.g., opening, reading, writing, flushing).
///
/// # Panics
///
/// This function will panic if:
/// - It fails to convert the extracted chromosome bytes to a UTF-8 string,
///   suggesting malformed input.
/// - It fails to open or create an output file for a specific chromosome,
///   likely due to permissions or an invalid path.
///
/// # Example
///
/// ```rust, no_run
/// use std::fs::{self, File};
/// use std::io::Write;
///
/// // Create a temporary directory for output
/// let temp_output_dir = std::env::temp_dir().join("test_split_chr");
/// let output_path = temp_output_dir.join("output");
///
/// std::fs::create_dir_all(&output_path).unwrap();
///
/// // Create dummy input files
/// let mut file1 = File::create(output_path.join("input1.bed")).unwrap();
/// file1
///     .write_all(b"chrA\tdata1\nchrB\tdata2\nchrA\tdata3\n")
///     .unwrap();
///
/// let mut file2 = File::create(output_path.join("input2.bed")).unwrap();
/// file2.write_all(b"chrC\tdata4\nchrB\tdata5\n").unwrap();
///
/// let files_to_split = vec![
///     output_path.join("input1.bed"),
///     output_path.join("input2.bed"),
/// ];
///
/// let suffix = "split.bed";
///
/// let result = __split_by_chr(&files_to_split, &output_path, suffix);
/// assert!(result.is_ok());
///
/// // Verify output files
/// let chr_a_file_content = fs::read_to_string(output_path.join("chrA_split.bed")).unwrap();
/// assert_eq!(chr_a_file_content, "chrA\tdata1\nchrA\tdata3\n");
///
/// let chr_b_file_content = fs::read_to_string(output_path.join("chrB_split.bed")).unwrap();
/// assert_eq!(chr_b_file_content, "chrB\tdata2\nchrB\tdata5\n");
///
/// let chr_c_file_content = fs::read_to_string(output_path.join("chrC_split.bed")).unwrap();
/// assert_eq!(chr_c_file_content, "chrC\tdata4\n");
///
/// std::fs::remove_dir_all(temp_output_dir).unwrap();
/// ```
pub fn __split_by_chr(files: &[PathBuf], dir: &Path, suffix: &str) -> std::io::Result<()> {
    let mut writers: HashMap<String, BufWriter<File>> = HashMap::new();

    for path in files {
        let f = File::open(path)?;
        let mut rdr = BufReader::new(f);
        let mut line = Vec::with_capacity(1 << 16);

        while {
            line.clear(); // clear buffer before each read
            rdr.read_until(b'\n', &mut line)?
        } != 0
        {
            let chr = {
                // INFO: find first TAB; fallback to the whole line (trim \n)
                let end = memchr(b'\t', &line).unwrap_or_else(|| {
                    line.iter()
                        .rposition(|&b| b != b'\n')
                        .map(|i| i + 1)
                        .unwrap_or(0)
                });

                // SAFETY: chromosome ids are usually ASCII; if not panic
                std::str::from_utf8(&line[..end])
                    .unwrap_or_else(|_| {
                        panic!(
                            "ERROR: could not parse chromosome from line {:?} in {:?}",
                            line, path
                        )
                    })
                    .trim()
                    .to_string()
            };

            let w = writers.entry(chr.clone()).or_insert_with(|| {
                let outfile = dir.join(format!("{}_{}", &chr, suffix));
                let file = OpenOptions::new()
                    .create(true)
                    // .write(true)
                    .append(true) // <- append to existing file [because SG is now in good]
                    // .truncate(true) // <- overwrite any existing file from previous runs
                    .open(outfile)
                    .unwrap_or_else(|_| {
                        panic!(
                            "ERROR: could not open file for writing {}",
                            dir.join(format!("{}_{}", &chr, suffix)).display()
                        )
                    });

                BufWriter::new(file)
            });

            w.write_all(&line)?;
            if !line.ends_with(b"\n") {
                w.write_all(b"\n")?;
            }
        }
    }

    for w in writers.values_mut() {
        w.flush()?;
    }

    if writers.is_empty() {
        log::warn!("WARN [SPLIT]: No chromosomes found in input files, no output files created.");
        panic!("ERROR: No chromosomes found in input files, no output files created.");
    } else {
        log::info!(
            "INFO [SPLIT]: Successfully split input files into {} chromosome-specific files.",
            writers.len()
        );

        Ok(())
    }
}
