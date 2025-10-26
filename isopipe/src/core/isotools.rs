use config::{CHUNK_SIZE, FUSION_FAKES, FUSION_FREE, FUSION_REVIEW};
use iso_polya::{
    cli::AparentArgs,
    core::apa::{calculate_polya, write_bed, RAM_PER_SITE},
};
use packbed::par_reader;

use std::{
    collections::HashMap,
    fmt::Debug,
    fs::remove_dir_all,
    io::{BufRead, BufReader, Write},
    path::{Path, PathBuf},
};

use crate::{config::*, consts::*, executor::job::Job};
use crate::{executor::manager::ParallelExecutor, isotools};

/// This function processes input directories to identify specific prediction files
/// (`tmp_predictions.bed`) and generates a series of commands (jobs) for the
/// `isotools iso-nmd` utility. It iterates through a hierarchical directory
/// structure, extracts chromosome and chunk information from directory names,
/// and constructs a command for each relevant file.
///
/// # Arguments
/// * `step` - A reference to a `PipelineStep` object, representing the current step in the pipeline.
///   It's used to retrieve custom configuration fields, such as `PARALLEL_MODE`.
/// * `config` - A reference to a `Config` object, which holds the overall configuration for the
///   pipeline. This is used to access step-specific custom fields.
/// * `input_dir` - A reference to a `PathBuf` indicating the root directory where the function
///   should search for input data. This directory is expected to contain subdirectories structured
///   as `{step_orf}/seqs_{suffix}` or `toga`.
/// * `step_output_dir` - A reference to a `PathBuf` specifying the directory where the output of
///   the `iso-nmd` commands will be stored.
///
/// # Returns
/// A `Vec<Job>`: A vector containing `Job` objects, where each `Job` represents a command to be
/// executed by the pipeline. These commands are generated based on the `tmp_predictions.bed`
/// files found in the input directory structure.
#[allow(clippy::map_entry)]
pub fn iso_nmd(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &Path,
) -> Vec<Job> {
    let _mode = ParallelMode::from_str(&config.get_step_custom_field(step, PARALLEL_MODE));
    let mut jobs = Vec::new();

    let keep_temp = config
        .get_step_custom_field(step, KEEP_TEMP)
        .parse::<bool>()
        .unwrap_or(false);

    if keep_temp {
        log::info!("INFO [STEP 9]: Keeping temporary files -> keep_temp set to true!");
    }

    let mut predictions = HashMap::new();
    let mut out_predictions = HashMap::new();

    // INFO: path would look like: {step_orf}/seqs_{suffix} or toga
    for entry in std::fs::read_dir(input_dir)
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", input_dir))
        .flatten()
        .filter(|e| e.path().is_dir())
    {
        let subdir = entry.path(); // INFO: seqs_{suffix}
        let suffix = subdir
            .file_name()
            .unwrap_or_else(|| panic!("ERROR: could not get file name from {:?}", subdir))
            .to_str()
            .unwrap_or_else(|| panic!("ERROR: could not convert file name to str"));

        // INFO: structure {step_orf}/seqs_{suffix}/{chr}:{chunk}
        for chunk in std::fs::read_dir(&subdir)
            .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", subdir))
            .flatten()
            .filter(|e| e.path().is_dir())
        {
            let chunked_dir = chunk.path(); // INFO: {chr}:{chunk}

            // WARN: need to check that tmp.predictions.bed exists and is not empty!
            if !chunked_dir.join("tmp.predictions.bed").exists() {
                log::warn!("WARN: skipping chunk {chunked_dir:?} because it does not have tmp.predictions.bed!");
                continue;
            } else if std::fs::metadata(chunked_dir.join("tmp.predictions.bed"))
                .unwrap_or_else(|e| {
                    panic!("ERROR: could not get tmp.predictions.bed for {chunked_dir:?} but seems to exist -> {e}")
                })
                .len()
                == 0
            {
                log::warn!("WARN: skipping chunk {chunked_dir:?} because it is empty!");
                continue;
            }

            let chr = chunked_dir
                .file_name()
                .unwrap_or_else(|| panic!("ERROR: could not get file name from {:?}", chunked_dir))
                .to_str()
                .unwrap_or_else(|| panic!("ERROR: could not convert file name to str"))
                .split(':')
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: could not get chromosome from {:?}", chunked_dir)
                });
            let chr_outdir = step_output_dir.join(chr).join(suffix);
            std::fs::create_dir_all(&chr_outdir).unwrap_or_else(|e| {
                panic!(
                    "ERROR: Failed to create directory -> {} -> {e}",
                    chr_outdir.display()
                )
            });

            for file in std::fs::read_dir(&chunked_dir)
                .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", chunked_dir))
                .flatten()
            {
                let file = file.path();

                if file.is_file() {
                    // INFO: matching specific name
                    if file.ends_with("tmp.predictions.bed") {
                        log::debug!("DEBUG: found predictions.bed for {chr:?} -> {file:?}");

                        let prediction = file.with_extension("tsv");

                        predictions
                            .entry(format!("{}.{}", chr, suffix))
                            .or_insert_with(Vec::new)
                            .push(prediction.clone());

                        let prefix =
                            format!("tmp_{}", chunked_dir.file_name().unwrap().to_str().unwrap());
                        let mut cmd = format!(
                            "{} --ref {} --outdir {} --prefix {}",
                            isotools!(ISO_NMD).display(),
                            file.display(),
                            chr_outdir.display(),
                            prefix,
                        );

                        // INFO: if out_predictions does not have a key, set it and do in-place cleaning
                        if !out_predictions.contains_key(&format!("{}.{}", chr, suffix)) {
                            out_predictions.insert(
                                format!("{}.{}", chr, suffix),
                                chr_outdir.join(format!("{}.predictions.{}.tsv", chr, suffix)),
                            );

                            // INFO: in-place duplicated header cleaning
                            cmd = format!(
                                    "{cmd} && gawk -i inplace 'NR==1 || $0 != header {{print}} NR==1 {{header=$0}}' {}",
                                    out_predictions
                                        .get(&format!("{}.{}", chr, suffix))
                                        .unwrap()
                                        .display()
                                );
                        }

                        // INFO: takes care of step8 tmp files except predictions
                        if !keep_temp {
                            let rest = format!(
                                r" && find {} \( -name '*fa' -or -name '*reduced*' -or -name '*index' \) -exec rm {{}} ';'",
                                file.parent().unwrap().display()
                            );

                            cmd += &rest;
                        }

                        log::debug!("DEBUG: executing cmd: {cmd}");
                        jobs.push(Job::from(cmd));
                    }
                }
            }
        }

        if !predictions.is_empty() && !out_predictions.is_empty() {
            for (key, target) in out_predictions.iter() {
                let childs = predictions.get(key).unwrap_or_else(|| {
                    panic!("ERROR: could not get predictions for key {key:?} in {predictions:?}")
                });

                log::debug!(
                    "DEBUG: merging predictions for {childs:?} for {suffix:?} in {entry:?} -> {target:?}"
                );
                cat(childs, target).unwrap_or_else(|e| {
                    panic!("ERROR: could not concatenate predictions for {childs:?} to {target:?} -> {e}")
                });
            }
        } else {
            log::warn!("WARN: No predictions found in {:?} -> skipping...", entry);
        }
    }

    log::info!("INFO [STEP 9]: Pre-processing completed -> Running...");
    jobs
}

/// Run isotools iso-fusion
///
/// # Arguments
/// * `step` - The pipeline step
/// * `config` - The configuration
/// * `input_dir` - The input directory
/// * `step_output_dir` - The output directory
///
/// # Returns
/// A vector of jobs
///
/// # Example
/// ```rust, ignore
/// let jobs = iso_fusion(&step, &config, &input_dir, &step_output_dir);
/// ```
pub fn iso_fusion(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &Path,
) -> Vec<Job> {
    let mut non_cannonical = true;
    let mut jobs = Vec::new();

    let refs = config.get_step_custom_field(step, TOGA);
    let keep_rejected = config
        .get_step_custom_field(step, KEEP_REJECTED)
        .parse::<bool>()
        .unwrap_or(false);

    // INFO: synchronizing IDs across all .bed files in input_dir
    __sync_ids(input_dir, keep_rejected);

    // INFO: input_dir would have per-chr .bed files with polya suffixes
    for file in std::fs::read_dir(input_dir)
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", input_dir))
        .flatten()
        .filter(|entry| {
            entry
                .path()
                .file_name()
                .and_then(|ext| ext.to_str())
                .map(|name| name.ends_with(BED))
                .unwrap_or(false)
        })
    {
        if file
            .path()
            .file_name()
            .and_then(|ext| ext.to_str())
            .unwrap()
            .ends_with(ALN_POLYA_REJECT)
            && !keep_rejected
        {
            continue;
        }

        // WARN: will need to change for the corrected minimap step suffix!
        let query = file.path();

        if !std::path::Path::new(&query).exists() {
            log::warn!("WARN: {} does not exist, skipping...", query.display());
            continue;
        }

        // INFO: file stem should be 'chrX_all.aligned.accept' -> dir name
        let prefix = step_output_dir.join(
            file.path()
                .file_stem()
                .unwrap_or_else(|| panic!("ERROR: could not build prefix for {:?}", file.path())),
        );

        std::fs::create_dir_all(&prefix).unwrap_or_else(|e| {
            panic!(
                "ERROR: Failed to create directory -> {} -> {e}",
                prefix.display()
            )
        });

        if file
            .path()
            .file_name()
            .and_then(|ext| ext.to_str())
            .unwrap()
            .ends_with(ALN_POLYA_SGN)
        {
            let cmd = format!(
                "{} --ref {} --query {} --prefix {} --tag {} --colorize {} --recover", // INFO: recover is mandatory to recognize real fusioned loci
                isotools!(ISO_FUSION).display(),
                refs,
                query.display(),
                prefix.display(),
                SINGLETONS,
                SGN_COLOR
            );

            log::debug!("DEBUG: processing {file:?} with cmd: {cmd:?}...");
            jobs.push(Job::from(cmd));
        } else {
            let cmd = format!(
                "{} --ref {} --query {} --prefix {} --recover",
                isotools!(ISO_FUSION).display(),
                refs,
                query.display(),
                prefix.display(),
            );

            log::debug!("DEBUG: processing {file:?} with cmd: {cmd:?}...");
            jobs.push(Job::from(cmd));
        }

        non_cannonical = false;
    }

    log::info!("INFO [STEP 7]: Pre-processing completed -> Running...");

    if non_cannonical {
        let nc_jobs = __build_non_cannonical_fusions(input_dir, step_output_dir, refs);
        return nc_jobs;
    }

    jobs
}

/// Synchronize read IDs across all BED files in a directory by replacing them with sequential numbers.
///
/// This function processes all BED files in the given directory, replacing the read IDs
/// with sequentially numbered IDs (R1, R2, R3, etc.) to ensure unique identifiers across
/// all files. Files are processed per chromosome and can optionally include rejected reads.
///
/// # Arguments
/// * `input_dir` - Directory containing BED files to process
/// * `keep_rejected` - Whether to include files with rejected reads in processing
///
/// # Panics
/// Panics if files cannot be opened, read, written, or renamed during processing
///
/// # Example
/// ```
/// sync_ids("/path/to/bed/files", true);
/// ```
pub fn __sync_ids<P: AsRef<Path> + Debug + Copy>(input_dir: P, keep_rejected: bool) {
    log::info!("INFO: start synchronizing read IDs in {input_dir:?}...");

    // INFO: chrX -> [chrX_all.aligned.accept.bed, chrX_all.aligned.singletons.bed, Option<reject>]
    let chroms = __collect_pairs(input_dir, keep_rejected);
    for (chr, files) in chroms {
        log::debug!("DEBUG: synchronizing read IDs for {chr:?} in {files:?}...");
        let mut counter = 0usize;

        for file in files {
            log::debug!("DEBUG: processing file {file:?}...");
            let reader = BufReader::new(
                std::fs::File::open(&file)
                    .unwrap_or_else(|e| panic!("ERROR: could not open file {:?} -> {e}", &file)),
            );

            let tmp = file.with_extension("tmp");
            log::debug!("DEBUG: using temporary file {tmp:?}...");

            let mut writer = std::fs::File::create(&tmp).unwrap_or_else(|e| {
                panic!("ERROR: could not create temporary file {:?} -> {e}", &tmp)
            });

            for line in reader.lines() {
                let mut l = line.unwrap_or_else(|e| {
                    panic!(
                        "ERROR: could not read line {} in {:?} -> {e}",
                        counter + 1,
                        &file
                    )
                });

                // INFO: replace R<number> with R<counter>
                l = __rename_id(&l, counter + 1);
                writeln!(writer, "{}", l).unwrap_or_else(|e| {
                    panic!(
                        "ERROR: could not write line {} to temporary file {:?} -> {e}",
                        counter + 1,
                        &tmp
                    )
                });
                counter += 1;
            }

            log::debug!("DEBUG: synchronized {counter} read IDs in {tmp:?}, now renaming from {tmp:?} to {file:?}...");
            std::fs::rename(&tmp, &file).unwrap_or_else(|e| {
                panic!(
                    "ERROR: could not rename temporary file {:?} to {:?} -> {e}",
                    tmp, file
                )
            });
        }
    }
}

/// Collect and group BED files by chromosome for synchronization processing.
///
/// Scans the input directory for BED files and groups them by chromosome name.
/// Files are expected to follow the naming pattern: `chrX_all.aligned.*.bed`
/// where X is the chromosome identifier.
///
/// # Arguments
/// * `input_dir` - Directory to scan for BED files
/// * `keep_rejected` - Whether to include files containing rejected reads
///
/// # Returns
/// HashMap mapping chromosome names to vectors of file paths
///
/// # Panics
/// Panics if the directory cannot be read or file names cannot be parsed
///
/// # Example
/// ```
/// let chrom_files = __collect_pairs("/path/to/files", false);
/// ```
pub fn __collect_pairs<P: AsRef<Path> + Debug + Copy>(
    input_dir: P,
    keep_rejected: bool,
) -> HashMap<String, Vec<PathBuf>> {
    // INFO: input_dir would have per-chr .bed files with polya suffixes
    let mut chroms = HashMap::new();
    for file in std::fs::read_dir(input_dir)
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", input_dir))
        .flatten()
        .filter(|entry| {
            entry
                .path()
                .file_name()
                .and_then(|ext| ext.to_str())
                .map(|name| name.ends_with(BED))
                .unwrap_or(false)
        })
    {
        if file
            .path()
            .file_name()
            .and_then(|ext| ext.to_str())
            .unwrap()
            .ends_with(ALN_POLYA_REJECT)
            && !keep_rejected
        {
            continue;
        }

        // WARN: will need to change for the corrected minimap step suffix!
        let query = file.path();

        if !std::path::Path::new(&query).exists() {
            log::warn!("WARN: {} does not exist, skipping...", query.display());
            continue;
        }

        // INFO: file stem should now be 'chrX_all.aligned.accept'
        let chr = query
            .file_stem()
            .unwrap_or_else(|| panic!("ERROR: could not get file name from {:?}", query))
            .to_str()
            .unwrap_or_else(|| panic!("ERROR: could not convert file name to str"))
            .split("_all")
            .next()
            .unwrap_or_else(|| panic!("ERROR: could not get chromosome from {:?}", query));

        // INFO: chrX -> [chrX_all.aligned.accept.bed, chrX_all.aligned.singletons.bed]
        chroms
            .entry(chr.to_string())
            .or_insert(Vec::new())
            .push(query.clone());
    }

    log::debug!("DEBUG: collected chroms -> {chroms:?}");
    chroms
}

/// Rename the read ID in a BED file line with a sequential identifier.
///
/// Replaces the original read name in the fourth column of a BED line
/// with a formatted sequential ID while preserving any suffix after
/// the first underscore in the original read name.
///
/// # Arguments
/// * `line` - The BED file line to process
/// * `new_id` - The sequential ID number to use (will be formatted as "R{new_id}")
///
/// # Returns
/// Modified line with updated read ID, or original line if format is unexpected
///
/// # Example
/// ```
/// let modified = __rename_id("chr1\t100\t200\toriginal_read_extra\t0\t+", 42);
/// // Returns: "chr1\t100\t200\tR42_extra\t0\t+"
/// ```
fn __rename_id(line: &str, new_id: usize) -> String {
    // INFO: find 3rd tab (before the read name column)
    let mut tabs = line.match_indices('\t');
    let third_tab = tabs.nth(2).map(|(i, _)| i);

    if let Some(start) = third_tab {
        // INFO: find the next tab after the read name
        let end = line[start + 1..]
            .find('\t')
            .map(|j| start + 1 + j)
            .unwrap_or(line.len());

        let read_name = &line[start + 1..end];

        // INFO: find the first underscore in the read name
        if let Some(pos) = read_name.find('_') {
            let mut out = String::with_capacity(line.len() + 8);
            out.push_str(&line[..start + 1]); // INFO: keep up to and including 3rd tab
            out.push_str(&format!("R{}", new_id));
            out.push_str(&read_name[pos..]); // INFO: append rest
            out.push_str(&line[end..]); // INFO: append the rest of the BED line
            return out;
        }
    }

    line.to_string()
}

/// Aggregate fusions from all categories into a single file.
///
/// # Arguments
/// * `step_output_dir` - The output directory
///
/// # Returns
/// A vector of jobs
///
/// # Example
/// ```
/// let jobs = aggregate_fusions(&step_output_dir);
/// ```
pub fn agg_fusions(
    dir: &PathBuf,
    executor: &mut ParallelExecutor,
    config: &Config,
    step: &PipelineStep,
) {
    let mode = ParallelMode::from_str(&config.get_step_custom_field(step, PARALLEL_MODE));

    let jobs: Vec<Job> = match mode {
        ParallelMode::Chromosome => {
            log::info!(
                "INFO [STEP 8]: Aggregating fusions per-chr -> {}...",
                dir.display()
            );

            let mut jobs = Vec::new();

            // INFO: loop through all subdirectories in {dir}
            for entry in std::fs::read_dir(dir)
                .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", dir))
                .flatten()
                .filter(|e| e.path().is_dir())
            {
                let subdir = entry.path();

                // INFO: structure should be {step_fusion}/chr{chr}_all.aligned.accept/fusions*
                // INFO: expected: fusions.bed / fusions.free.bed; optional: fusions.review + fusions.fakes
                for file in std::fs::read_dir(&subdir)
                    .unwrap()
                    .flatten()
                    .filter(|e| e.path().is_file())
                {
                    let file = file.path();

                    if file.ends_with(FUSION_FAKES) || file.ends_with(FUSION_REVIEW) {
                        let target = file
                            .parent()
                            .unwrap_or_else(|| {
                                panic!("ERROR: could not get parent from {:?}", file)
                            })
                            .join(FUSION_FREE);

                        let cmd = format!(
                            "cat {} >> {} && rm {}",
                            file.display(),
                            target.display(),
                            file.display()
                        );

                        jobs.push(Job::from(cmd))
                    }
                }
            }

            jobs
        }
        ParallelMode::Genome => {
            let jobs = FUSION_TYPES
                .iter()
                .map(|ty| {
                    log::info!(
                        "INFO [STEP 8]: Aggregating {} in single files -> {}...",
                        ty,
                        dir.display()
                    );

                    let (output, pattern) = if *ty == "fusions" {
                        (
                            format!("{0}/{1}.bed", dir.display(), ty),
                            format!("{0}/*/{1}.bed", dir.display(), ty),
                        )
                    } else {
                        // WARN: will aggregate free/fakes/review into fusions.free.bed!
                        // INFO: fakes have :FK tag, review has :RW tag
                        (
                            format!("{0}/fusions.free.bed", dir.display()),
                            format!("{0}/*/*.{1}.bed", dir.display(), ty),
                        )
                    };

                    Job::from(format!("cat {} >> {}", pattern, output))
                })
                .collect();

            jobs
        }
    };

    executor
        .add_jobs(jobs)
        .and_send(config, "agg-fusions", dir.clone(), 1, 8, None, None);
}

/// Build non-cannonical fusions
///
/// # Arguments
///
/// * `input_dir` - The input directory
/// * `step_output_dir` - The output directory
/// * `parts` - The parts to be used
///
/// # Example
///
/// ```rust, no_run
/// let input_dir = PathBuf::from("/path/to/input_dir");
/// let step_output_dir = PathBuf::from("/path/to/step_output_dir");
/// let parts = vec![String::from("--part1"), String::from("--part2")];
///
/// __build_non_cannonical_fusions(&input_dir, &step_output_dir, parts);
/// ```
fn __build_non_cannonical_fusions(
    input_dir: &PathBuf,
    step_output_dir: &Path,
    refs: String,
) -> Vec<Job> {
    log::warn!(
        "WARN: No cannonical jobs found to run for isotools fusion in {} -> trying to grab any .bed!",
        input_dir.display()
    );

    let mut jobs = Vec::new();

    for (idx, entry) in std::fs::read_dir(input_dir)
        .expect("Failed to read assets directory")
        .flatten()
        .filter(|entry| {
            entry
                .path()
                .file_name()
                .and_then(|ext| ext.to_str())
                .map(|name| name.ends_with(BED))
                .unwrap_or(false)
        })
        .enumerate()
    {
        let query = entry.path();
        let prefix = step_output_dir.join(format!("reads_{}", idx));

        std::fs::create_dir_all(&prefix).unwrap_or_else(|e| {
            panic!(
                "ERROR: Failed to create directory -> {} -> {e}",
                prefix.display()
            )
        });

        let cmd = format!(
            "{} --ref {} --query {} --prefix {}",
            isotools!(ISO_FUSION).display(),
            refs,
            query.display(),
            prefix.display(),
        );

        jobs.push(Job::from(cmd));
    }

    if jobs.is_empty() {
        log::error!(
            "ERROR: No .bed files found in {} to run non-cannonical fusion!",
            input_dir.display()
        );
        std::process::exit(1);
    } else {
        log::info!(
            "INFO [STEP 7]: Found {} non-cannonical fusion files in {}",
            jobs.len(),
            input_dir.display(),
        );
    }

    jobs
}

/// Run iso-polya aparent
///
/// # Arguments
///
/// * `executor` - The parallel executor
/// * `config` - The configuration
/// * `input_dir` - The input directory
/// * `step_output_dir` - The output directory
/// * `step` - The pipeline step
///
/// # Example
///
/// ```rust, no_run
/// let executor = ParallelExecutor::new();
/// let config = Config::new();
/// let input_dir = PathBuf::from("/path/to/input_dir");
/// let step_output_dir = PathBuf::from("/path/to/step_output_dir");
/// let step = PipelineStep::new();
///
/// iso_polya_aparent(&executor, &config, &input_dir, &step_output_dir, &step);
/// ```
fn iso_polya_aparent(step_output_dir: &Path, bed: &str, twobit: &str, prefix: &str) -> Vec<Job> {
    log::info!("INFO [STEP 10a]: Pre-processing for APARENT...");

    let mut jobs = Vec::new();

    let outdir = step_output_dir.join(prefix);

    std::fs::create_dir_all(&outdir).unwrap_or_else(|e| {
        panic!(
            "ERROR: Failed to create directory -> {} -> {e}",
            outdir.display()
        )
    });

    let args = vec![
        String::from("aparent"),
        String::from("--bed"),
        bed.to_string(),
        String::from("--twobit"),
        twobit.to_string(),
        String::from("--outdir"),
        outdir.display().to_string(),
    ];

    let accumulator = calculate_polya(AparentArgs::from(args))
        .expect("ERROR: Failed to calculate polyA apparent");

    for path in accumulator.paths.iter() {
        let cmd = format!(
            "source {} && python3 {} -p {}",
            TAI_VENV,
            APARENT_PY,
            path.as_str()
        );

        jobs.push(Job::from(cmd));
    }

    jobs
}

/// Run isotools on final ORF-called data
/// [iso-classify, iso-intron, iso-polya, iso-utr]
///
/// # Arguments
///
/// * `step` - The pipeline step
/// * `config` - The configuration
/// * `input_dir` - The input directory
/// * `step_output_dir` - The output directory
/// * `executor` - The parallel executor
///
/// # Returns
///
/// A vector of jobs
///
/// # Example
///
/// ```rust
/// let jobs = polish(&step, &config, &input_dir, &step_output_dir, &executor);
/// ```
pub fn polish(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &Path,
    executor: &mut ParallelExecutor,
) -> Vec<Job> {
    let mut jobs = Vec::new();

    let args = config.get_step_args(
        step,
        vec![
            INPUT_DIR,
            OUTPUT_DIR,
            MEMORY,
            TIME,
            GENOME,
            NUM_THREADS,
            KEEP_TEMP,
        ],
    );

    let keep_temp = config
        .get_step_custom_field(step, KEEP_TEMP)
        .parse::<bool>()
        .unwrap_or(false);

    if keep_temp {
        log::info!("INFO [STEP 10]: Keeping temporary files -> keep_temp set to true!");
    }

    let twobit = config.get_step_custom_fields(step, vec![GENOME])[0].clone();

    // INFO: cleaning input_dir to avoid calling cmds on empty files/dirs
    let _clean = format!(
        "find {} -type f -name '*.bed' -empty -delete && find {} -type d -empty -delete",
        input_dir.display(),
        input_dir.display()
    );
    log::debug!("DEBUG: cleaning input directory with cmd: {_clean}");

    std::process::Command::new("bash")
        .arg("-c")
        .arg(_clean)
        .output()
        .unwrap_or_else(|e| panic!("ERROR: Failed to clean input directory -> {e}"));

    // INFO: looping again to run polish with merged aparent
    __pre_polish(input_dir, step_output_dir, executor, config, &twobit);

    // INFO: looping again to run polish with merged aparent -> per-chr
    for entry in std::fs::read_dir(input_dir)
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", input_dir))
        .flatten()
        .filter(|e| e.path().is_dir())
    {
        let subdir = entry.path(); // INFO: {chr}
        let chr =
            subdir.file_name().unwrap().to_str().unwrap_or_else(|| {
                panic!("ERROR: could not convert chr name to str -> {subdir:?}")
            });

        // INFO: again, we only want to run isotools on {chr}/seqs_free
        let bed = subdir.join("seqs_free").join(format!("{}.reads.bed", chr));
        let outdir = step_output_dir.join(chr);

        let preds = subdir.join("*/*predictions*tsv"); // INFO: move predictions from free/fusions

        let fsn = subdir.join("seqs_fusions").join(format!("{chr}.reads.bed"));
        let reads = subdir.join("seqs_free").join(format!("{chr}.reads.bed"));

        let mut cleaning = format!(
            "cat {} > {} && gawk -i inplace 'NR==1 || $0 != header {{print}} NR==1 {{header=$0}}' {}",
            preds.display(),
            step_output_dir.join(chr).join("predictions.tsv").display(),
            step_output_dir.join(chr).join("predictions.tsv").display(),
        );

        // WARN: not asserting existence because iso-nmd will create a file either way -> WRONG
        // WARN: existence of nmds should be checked because if NMD collection is empty no file is created
        // let nmd = subdir.join("*/*nmd.bed"); // INFO: move nmds from free/fusions
        let free_nmd = subdir.join("seqs_free").join(format!("{}.nmd.bed", chr));
        let fusion_nmd = subdir.join("seqs_fusions").join(format!("{}.nmd.bed", chr));
        let nmd_target = step_output_dir.join(chr).join("nmd.bed");

        if free_nmd.exists() {
            log::debug!("DEBUG: found seqs_free nmd file -> {free_nmd:?} -> will concatenate to {nmd_target:?}");
            cleaning += &format!(" && cat {} >> {}", free_nmd.display(), nmd_target.display());
        }

        if fusion_nmd.exists() {
            log::debug!("DEBUG: found seqs_fusions nmd file -> {fusion_nmd:?} -> will concatenate to {nmd_target:?}");
            cleaning += &format!(
                " && cat {} >> {}",
                fusion_nmd.display(),
                nmd_target.display()
            );
        }

        // INFO: subdirs should have then -> nmd.bed, Option<fusions.bed>, Option<raw_reads.bed>
        // INFO: nmd is assumed to exist either way -> otherwise would be catched by dir deleting cmd above
        // INFO: an optionally -> {chr}.predictions.seqs_fusions.tsv, {chr}.predictions.seqs_free.tsv
        let mut decisions = format!("source {}", TAI_VENV);

        if free_nmd.exists() || fusion_nmd.exists() {
            decisions += &format!(
                " && {} --reads {} --predictions {} --name nmd.bed --nmd --outdir {}",
                PRETTY_PY,
                step_output_dir.join(chr).join("nmd.bed").display(),
                step_output_dir.join(chr).join("predictions.tsv").display(),
                step_output_dir.join(chr).display(),
            );
        }

        if fsn.exists() {
            cleaning += &format!(
                " && cp {} {}",
                fsn.display(),
                step_output_dir.join(chr).join("fusions.bed").display()
            );

            decisions += &format!(
                "  && {} --reads {} --predictions {} --name fusions.bed --outdir {}",
                PRETTY_PY,
                step_output_dir.join(chr).join("fusions.bed").display(),
                step_output_dir.join(chr).join("predictions.tsv").display(),
                step_output_dir.join(chr).display()
            )
        } else {
            log::warn!(
                "WARN: {:?} does not exist, skipping...",
                step_output_dir.join(chr).join("fusions.bed")
            );
        }

        // INFO: aparent is only ran per-seqs_free-{chr} -> not running isotools on seqs_fusions!
        let apa = if let Some(file) = merge_aparent(step_output_dir.join(chr), "tmp") {
            file
        } else {
            log::warn!(
                "WARN: could not find any APARENT output for {chr:?} in {} -> skipping isotools polish for this chr!",
                step_output_dir.join(chr).display()
            );

            // INFO: still need to execute decision on nmd and Option<fusions> + cleaning
            let cmd = format!("{} && {}", cleaning, decisions);
            log::debug!("DEBUG: executing non-APARENT + non-raw_reads dir: {cmd}");

            jobs.push(Job::from(cmd));
            continue;
        };

        if reads.exists() {
            cleaning += &format!(
                " && cp {} {}",
                reads.display(),
                step_output_dir.join(chr).join("raw_reads.bed").display()
            );

            let cmd = format!(
                " && {} --reads {} --predictions {} --introns {} --descriptor {} --outdir {}",
                PRETTY_PY,
                step_output_dir.join(chr).join("raw_reads.bed").display(),
                step_output_dir.join(chr).join("predictions.tsv").display(),
                step_output_dir
                    .join(chr)
                    .join("reference_introns.tsv")
                    .display(),
                step_output_dir
                    .join(chr)
                    .join("global_descriptor.tsv")
                    .display(),
                step_output_dir.join(chr).display(),
            );

            decisions += &cmd;
        }

        let mut cmd = format!(
            "{} run --query {} --aparent {} --twobit {} {args} --outdir {} && {} && {}",
            isotools!(ISOTOOLS).display(),
            bed.display(),
            apa.display(),
            twobit,
            outdir.display(),
            cleaning,
            decisions
        );

        if !keep_temp {
            // INFO: need to specify cleanup to include it in the cmd!
            let tmp = subdir.join("*"); // INFO: remove anything tmp

            let rest = format!(
                " && find {} -type f -name 'tmp*' -exec rm {{}} ';'",
                tmp.display()
            );
            cmd += &rest;
        }

        log::debug!("DEBUG: executing cmd: {cmd}");
        jobs.push(Job::from(cmd));
    }

    jobs
}

/// Pre-polish input read files to run isotools polish
/// on the raw reads and APARENT output
///
/// # Arguments
///
/// * `input_dir` - The input directory
/// * `step_output_dir` - The output directory
/// * `executor` - The parallel executor
/// * `config` - The configuration
/// * `twobit` - The twobit file
///
/// # Example
///
/// ```rust, no_run
/// let input_dir = PathBuf::from("/path/to/input_dir");
/// let step_output_dir = PathBuf::from("/path/to/step_output_dir");
/// let executor = ParallelExecutor::new();
/// let config = Config::new();
/// let twobit = String::from("/path/to/twobit");
///
/// __pre_polish(&input_dir, &step_output_dir, &executor, &config, &twobit);
/// ```
fn __pre_polish<P: AsRef<Path> + Debug + Copy>(
    input_dir: P,
    step_output_dir: &Path,
    executor: &mut ParallelExecutor,
    config: &Config,
    twobit: &str,
) {
    // let mem = CHUNK_SIZE as f32 * RAM_PER_SITE; // INFO: will be converted to MB by executor
    let mut inner_jobs = Vec::new();

    // INFO: path would look like: {step_nmd}/{chr} -> looping for each chr
    for entry in std::fs::read_dir(input_dir)
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", input_dir))
        .flatten()
        .filter(|e| e.path().is_dir())
    {
        let subdir = entry.path(); // INFO: {chr}

        let chr =
            subdir.file_name().unwrap().to_str().unwrap_or_else(|| {
                panic!("ERROR: could not convert chr name to str -> {subdir:?}")
            });

        // INFO: {chr}/seqs_{suffix}
        for seqs in std::fs::read_dir(&subdir)
            .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", subdir))
            .flatten()
            .filter(|e| e.path().is_dir())
        {
            let seqs_dir = seqs.path();
            let suffix = seqs_dir
                .file_name()
                .unwrap_or_else(|| panic!("ERROR: could not get file name from {:?}", seqs_dir))
                .to_str()
                .unwrap_or_else(|| panic!("ERROR: could not convert file name to str"));

            let beds = format!("find {} -type f -name '{}' -print0 | xargs -0 cat > {} && find {} -type f -empty -name '{}' -delete",
                seqs_dir.display(),
                "tmp*reads.bed",
                seqs_dir.join(format!("{chr}.reads.bed")).display(),
                seqs_dir.display(),
                "*reads.bed"
            );

            let nmds = format!("find {} -type f -name '{}' -print0 | xargs -0 cat > {} && find {} -type f -empty -name '{}' -delete",
                seqs_dir.display(),
                "tmp*nmd.bed",
                seqs_dir.join(format!("{chr}.nmd.bed")).display(),
                seqs_dir.display(),
                "*nmd.bed"
            );

            let cat = format!("{beds} && {nmds}");
            log::debug!("DEBUG: concatenating bed files with cmd: {cat}");

            std::process::Command::new("bash")
                .arg("-c")
                .arg(cat)
                .current_dir(&seqs_dir)
                .output()
                .expect("ERROR: Failed to concatenate bed files");

            if suffix == "seqs_free" {
                let bed = seqs_dir.join(format!("{chr}.reads.bed"));

                if !bed.exists() {
                    log::warn!(
                        "WARN: could not find {} -> skipping APARENT for {chr:?}!",
                        bed.display()
                    );
                    continue;
                }

                let apa_jobs =
                    iso_polya_aparent(step_output_dir, bed.to_str().unwrap(), twobit, chr);

                log::debug!("DEBUG: APARENT jobs for {chr:?} -> {:?}", apa_jobs.len());
                inner_jobs.extend(apa_jobs);
            }
        }
    }

    log::info!("INFO [STEP 10a]: Pre-processing completed -> Running APARENT...");

    executor.add_jobs(inner_jobs).and_send(
        config,
        "aparent",
        step_output_dir.to_path_buf(),
        1,
        4,
        None,
        None,
    );

    log::info!("INFO [STEP 10b]: Pre-processing completed -> Polishing...");
}

/// Merge chunks of APARENT data
///
/// # Arguments
///
/// * `outdir` - The output directory
///
/// # Example
///
/// ```rust, no_run
/// let outdir = PathBuf::from("/path/to/output_dir");
/// merge_aparent(&outdir);
/// ```
fn merge_aparent(outdir: PathBuf, _prefix: &str) -> Option<PathBuf> {
    let assets = outdir.join(APARENT_CHUNKS);

    if !assets.exists() {
        log::warn!(
            "WARN: could not find any APARENT chunks in {:?} -> skipping...",
            assets
        );
        return None;
    }

    let mut beds = Vec::new();

    for entry in std::fs::read_dir(assets.clone())
        .unwrap_or_else(|e| panic!("ERROR: could not read directory -> {:?}. {e}", assets))
        .flatten()
    {
        let path = entry.path();
        if let Some(ext) = path.extension() {
            if let Some("bed") = ext.to_str() {
                beds.push(path);
            }
        }
    }

    let bed = par_reader(beds).expect("ERROR: Failed to merge bed files");
    let bed_dest = outdir.join(APARENT_OUTPUT);

    if bed.is_empty() {
        log::warn!(
            "WARN: could not find any .bed files in {:?} to merge for APARENT -> will erase this dir!",
            assets
        );
        remove_dir_all(&outdir)
            .unwrap_or_else(|e| panic!("ERROR: {e} -> could not remove {outdir:?}"));

        return None;
    }

    log::debug!("DEBUG: will write all chunked .bed to -> {bed_dest:?}");
    write_bed(bed_dest.clone(), bed);

    log::info!("INFO: Merged chunks and cleaning...");
    for entry in std::fs::read_dir(&assets)
        .expect("ERROR: Failed to read assets directory")
        .flatten()
    {
        let path = entry.path();
        if path
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .starts_with("tmp")
        {
            let _ = std::fs::remove_file(path);
        }
    }

    remove_dir_all(&assets).unwrap_or_else(|e| panic!("ERROR: {e} -> could not remove {assets:?}"));

    log::info!("SUCCESS: APPARENT finished successfully!");
    Some(bed_dest)
}

/// Run isotools iso-split
///
/// # Arguments
/// * `step` - The pipeline step
/// * `config` - The configuration
/// * `input_dir` - The input directory
/// * `step_output_dir` - The output directory
///
/// # Returns
/// A vector of jobs
///
/// # Example
/// ```
/// let jobs = iso_fusion(&step, &config, &input_dir, &step_output_dir);
/// ```
pub fn iso_split(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &Path,
    executor: &mut ParallelExecutor,
) -> PathBuf {
    log::info!(
        "INFO [ISO-SPLIT]: Preparing to split files in {} and writing to {}/chunks",
        input_dir.display(),
        step_output_dir.display()
    );

    let chunks = crate::numerical!(config.get_step_custom_field(step, CHUNK) => usize)
        .expect("Missing or invalid chunk size");
    let outdir = step_output_dir.join(CHUNKS);

    // WARN: if any of the entries is FASTQ_GZ we assume the whole run is non-cannonical
    // WARN: would not make any sense combining both -> user should make two different runs and merge
    let entries = std::fs::read_dir(input_dir)
        .expect("ERROR: Failed to read input directory")
        .flatten()
        .map(|entry| entry.path())
        .filter(|path| {
            path.file_name()
                .and_then(|name| name.to_str())
                .map(|filename| {
                    filename.ends_with(FASTQ_GZ)
                        || filename.ends_with(FQ_GZ)
                        || filename.ends_with(HQ_FASTA_GZ)
                        || filename.ends_with(HQ_FA_GZ)
                        || filename.ends_with(SGN_FASTA_GZ)
                        || filename.ends_with(SGN_FA_GZ)
                })
                .unwrap_or(false)
        })
        .collect::<Vec<_>>();

    if entries.is_empty() {
        log::error!(
            "ERROR: could not find any .fa/.fq in {}!",
            input_dir.display()
        );
        std::process::exit(1);
    } else {
        log::info!(
            "INFO [ISO-SPLIT]: Found {} files in {}",
            entries.len(),
            input_dir.display(),
        );
    }

    split_reads(entries, &chunks, &outdir, executor, config, step, input_dir);

    outdir
}

/// Generic processor for a collection of input
/// FASTA/FASTQ files
///
/// # Arguments
///
/// * `files` - A slice of `PathBuf`s representing the input files.
/// * `chunks` - The number of chunks to split each file into.
/// * `outdir` - The output directory where the chunks should be written.
fn split_reads(
    files: Vec<PathBuf>,
    chunks: &usize,
    outdir: &Path,
    executor: &mut ParallelExecutor,
    config: &Config,
    step: &PipelineStep,
    input_dir: &Path,
) {
    log::info!(
        "INFO [ISO-SPLIT]: Running iso-split on {} files!",
        files.len()
    );

    let mut jobs = Vec::new();
    let mut global_threads = 0;

    files.iter().for_each(|file| {
        let threads = if file
            .file_name()
            .and_then(|name| name.to_str())
            .map(|name| name.ends_with(FASTQ_GZ))
            .unwrap_or(false)
        {
            1
        } else {
            16
        };

        // WARN: if any of the entries is FASTQ_GZ we assume the whole run is non-cannonical
        // WARN: would not make any sense combining both -> user should make two different runs and merge
        // WARN: here enforcing the amount of threads dependent on extension!
        if global_threads < 1 {
            global_threads = threads
        }

        let bind = file.with_extension("");
        let suffix = bind
            .file_stem()
            .unwrap_or_else(|| panic!("ERROR: could not build prefix for {:?}", file));

        let cmd = format!(
            "{} --file {} --chunks {} --outdir {} --suffix {} --threads {}",
            isotools!(ISO_SPLIT).display(),
            file.display(),
            chunks,
            outdir.display(),
            suffix.to_string_lossy(),
            threads,
        );

        jobs.push(Job::from(cmd));
    });

    executor.add_jobs(jobs).and_send(
        config,
        "iso-split",
        input_dir.to_path_buf(),
        global_threads,
        8,
        Some(config.get_package_from_step(step)),
        None,
    );
}
