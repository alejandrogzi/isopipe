use crate::executor::manager::ParallelExecutor;
use crate::{config::*, consts::*, executor::job::Job};

use config::{OverlapType, Sequence, Strand, SCALE};
use iso_polya::utils::get_sequences;
use packbed::{unpack, GenePred};
use rayon::prelude::*;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

/// Run ORF prediction on fusion files
///
/// # Arguments
///
/// * `step` - Pipeline step
/// * `config` - Pipeline configuration
/// * `input_dir` - Path to the input directory
/// * `step_output_dir` - Path to the output directory
///
/// # Returns
///
/// * Vec<Job> - List of jobs to be executed
///
/// # Example
///
/// ```rust, no_run
/// let jobs = orf(
///     &step,
///     &config,
///     &input_dir,
///     &step_output_dir,
/// );
/// ```
pub fn orf(
    step: &PipelineStep,
    config: &Config,
    executor: &mut ParallelExecutor,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
) -> Vec<Job> {
    let mut jobs = Vec::new();
    let args = config.get_step_custom_fields(step, vec![GENOME, ORFIPY, ORF_MIN_LEN, DATABASE]);

    let mode = ParallelMode::from_str(&config.get_step_custom_field(step, PARALLEL_MODE));

    let twobit = PathBuf::from(args[0].clone());
    let chunk_size = crate::numerical!(config.get_step_custom_field(step, CHUNK) => usize)
        .unwrap_or_else(|e| panic!("ERROR: could not convert chunk to num -> {e}!"));

    log::info!(
        "INFO: Merging TOGA predictions in a single file here: {}...",
        step_output_dir.display()
    );
    __merge_toga(step_output_dir, config, step);

    // INFO: looping through all [per-chr or merged -> depend on parallel_mode opt]
    // INFO: each [subdir/main] should have -> free + fakes + review [other color + tag] and fusions
    match mode {
        ParallelMode::Chromosome => {
            unbounded_extract(
                config,
                executor,
                input_dir,
                &twobit,
                step_output_dir,
                chunk_size,
            );

            process_bed(
                None,
                &twobit,
                step_output_dir,
                chunk_size,
                &args,
                &mut jobs,
                &mode,
            );
        }
        ParallelMode::Genome => {
            for file in FUSION_FILES.iter().take(2) {
                let bed = input_dir.join(file);
                process_bed(
                    Some(bed),
                    &twobit,
                    step_output_dir,
                    chunk_size,
                    &args,
                    &mut jobs,
                    &mode,
                );
            }
        }
    }

    return jobs;
}

/// Orchestrates the "unbounded" extraction process for a set of BED files,
/// generating jobs for parallel execution.
///
/// This function is designed to handle multiple BED files organized in a specific directory
/// structure. It iterates through subdirectories of `input_dir` (expected to be named `accept`
/// or `reject`), and for each subdirectory, it finds BED files to be processed. For each
/// found BED file, it constructs a command-line job to run a sequence extraction tool.
/// The constructed jobs are then added to a `ParallelExecutor` for parallel processing.
///
/// The extraction tool (`EXTRACT_RELEASE`) is configured with the following options:
/// - Path to the 2bit genome file.
/// - Path to the input BED file.
/// - An output directory for the results.
/// - `--index` flag to enable indexed extraction mode.
/// - `--suffix` to label the output files based on the subdirectory (`FREE` for `accept`,
///   `FUSIONS` for `reject`).
/// - `--chunk-size` to specify the number of records per chunk.
///
/// # Arguments
///
/// * `config` - A reference to the `Config` struct.
/// * `executor` - A mutable reference to a `ParallelExecutor` which manages the execution of jobs.
/// * `input_dir` - A `PathBuf` pointing to the directory containing `accept` and `reject` subdirectories.
/// * `twobit` - A `PathBuf` pointing to the 2bit genome file.
/// * `step_output_dir` - A `PathBuf` for the output directory of the extraction step.
/// * `chunk_size` - An integer specifying the desired chunk size for the extraction.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to read the `input_dir` directory.
/// - A subdirectory name does not end with `accept` or `reject`.
/// - It fails to read a subdirectory.
/// - It fails to create a `Job` from the command string.
fn unbounded_extract(
    config: &Config,
    executor: &mut ParallelExecutor,
    input_dir: &PathBuf,
    twobit: &PathBuf,
    step_output_dir: &PathBuf,
    chunk_size: usize,
) {
    let mut jobs = Vec::new();

    // INFO: structure should be {step_fusion}/chr{chr}_all.aligned.{accept/reject}
    for entry in std::fs::read_dir(input_dir)
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", input_dir))
        .flatten()
        .filter(|e| e.path().is_dir())
    {
        let subdir = entry.path();
        let suffix = get_suffix_from_path(&subdir);

        // INFO: structure should be {step_fusion}/chr{chr}_all.aligned.accept/fusions*
        // INFO: expected: fusions.bed / fusions.free.bed
        for bed in std::fs::read_dir(&subdir)
            .unwrap()
            .flatten()
            .filter(|e| e.path().is_file())
        {
            // INFO: collect paths and collect extract cmds
            // INFO: end path would look like: {step_orf}/seqs_{suffix}/{chr}:{chunk}/{name}{fa/bed}
            let cmd = format!(
                "{} --twobit {} --bed {} -o {} --index --suffix {} --chunk-size {}",
                EXTRACT_RELEASE,
                &twobit.display(),
                bed.path().display(),
                step_output_dir.display(),
                suffix,
                chunk_size
            );

            jobs.push(Job::from(cmd));
        }
    }

    // INFO: running inner_jobs
    executor
        .add_jobs(jobs)
        .and_send(config, EXTRACT, step_output_dir.clone(), 8, 32, None);
}

/// Processes a BED file by extracting sequences, chunking them, and generating
/// `blast` and `tai` jobs for ORF detection.
///
/// This function takes a BED file and a 2bit genome file, extracts sequences
/// corresponding to the BED entries, and then splits these extracted sequences
/// and their corresponding BED entries into smaller "chunks." For each chunk,
/// it constructs two command-line jobs: one for running `orf blast` and another
/// for `orf tai`. These jobs are then added to a provided list of `Job` objects.
/// Empty or non-existent BED files are logged as warnings and skipped.
///
/// # Arguments
///
/// * `bed` - A `PathBuf` representing the path to the input BED file.
/// * `twobit` - A `PathBuf` representing the path to the 2bit genome file
///              used for sequence extraction.
/// * `step_output_dir` - A `PathBuf` representing the base output directory
///                       for the current pipeline step. Chunked files will be
///                       placed within subdirectories of this path.
/// * `chunk_size` - An `usize` specifying the desired number of entries per chunk.
/// * `args` - A reference to a `Vec<String>` containing command-line arguments.
///            Specifically, `args[1]` is expected to be the `orfipy` executable path,
///            `args[2]` is the minimum ORF length, and `args[3]` is the BLAST database path.
/// * `jobs` - A mutable reference to a `Vec<Job>` where the generated
///            `blast` and `tai` jobs will be appended.
///
/// # Panics
///
/// This function will panic if:
/// - It cannot extract the file stem from the `bed` path.
/// - It cannot determine a suffix from the `bed` file name.
/// - The `extract` function (responsible for chunking) fails.
/// - It cannot determine the parent directory for a chunked BED file.
///
/// # Example
///
/// ```rust, no_run
/// process_bed(
///     input_bed_path.clone(),
///     &twobit_path,
///     &output_dir,
///     chunk_size,
///     &cli_args,
///     &mut job_list,
/// );
/// ```
fn process_bed(
    bed: Option<PathBuf>,
    twobit: &PathBuf,
    step_output_dir: &PathBuf,
    chunk_size: usize,
    args: &Vec<String>,
    jobs: &mut Vec<Job>,
    mode: &ParallelMode,
) {
    match mode {
        ParallelMode::Chromosome => parallel_processing(step_output_dir, args, jobs),
        ParallelMode::Genome => cannonical_processing(
            bed.unwrap(),
            twobit,
            step_output_dir,
            chunk_size,
            args,
            jobs,
        ),
    }
}

/// Prepares parallel jobs for BLAST and Translation AI (TAI) analysis on a set of
/// previously indexed and chunked files.
///
/// This function is designed to be run after an indexing and chunking step has been
/// completed. It iterates through a directory structure (e.g., `{step_orf}/seqs_{suffix}/{chr}:{chunk}`)
/// to find the necessary files for each chunk: a FASTA file, a reduced BED file, and an
/// index file.
///
/// For each chunk directory, it constructs two command-line jobs: one for the `blast`
/// subcommand and one for the `tai` subcommand of the `ORF_RELEASE` tool. It dynamically
/// builds the command strings by appending the paths to the FASTA, reduced BED, and index
/// files found within each chunk directory.
///
/// The `blast` command includes arguments for e-value, output directory, minimum ORF length,
/// and the protein database. The `tai` command similarly specifies the output directory.
/// Both commands are then enriched with the paths to the FASTA, reduced BED, and index files.
///
/// The resulting jobs are added to a `jobs` vector for parallel execution.
///
/// # Arguments
///
/// * `step_output_dir` - A `PathBuf` pointing to the root directory containing the chunked
///                       output from a previous step (e.g., `extract`).
/// * `args` - A reference to a `Vec<String>` containing command-line arguments,
///            specifically for `orfipy`'s e-value, minimum ORF length, and database path.
/// * `jobs` - A mutable reference to a `Vec<Job>` where the generated BLAST and TAI jobs are stored.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to read the `step_output_dir` or any of its subdirectories.
/// - It fails to construct a `Job` from the command string.
fn parallel_processing(step_output_dir: &PathBuf, args: &Vec<String>, jobs: &mut Vec<Job>) {
    // INFO: need to loop again to run blast and tai
    // INFO: end path would look like: {step_orf}/seqs_{suffix}/{chr}:{chunk}/{name}{fa/bed/reduced/index}
    for entry in std::fs::read_dir(step_output_dir)
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", step_output_dir))
        .flatten()
        .filter(|e| e.path().is_dir())
    {
        let subdir = entry.path(); // INFO: seqs_{suffix}

        // INFO: structure should be {step_orf}/seqs_{suffix}/{chr}:{chunk}/{fa/bed/reduced/index}
        // INFO: expected: reduced_bed/fa/index triplet
        for chunk in std::fs::read_dir(&subdir)
            .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", subdir))
            .flatten()
            .filter(|e| e.path().is_dir())
        {
            let chunked_dir = chunk.path(); // INFO: {chr}:{chunk}

            let mut blast = format!(
                "{} blast -e {} --outdir {} --orf-min-len {} --db {} ",
                ORF_RELEASE,
                args[1], // INFO: orfipy
                chunked_dir.display(),
                args[2], // INFO: orf_min_len,
                args[3], // INFO: database
            );

            let mut tai = format!("{} tai --outdir {} ", ORF_RELEASE, chunked_dir.display(),);

            for file in std::fs::read_dir(&chunked_dir)
                .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", chunked_dir))
                .flatten()
                .filter(|e| e.path().is_file())
            {
                let file = file.path();

                if file.ends_with(REDUCED_BED) {
                    let part = format!("--alignments {} ", file.display());
                    tai += &part;
                    blast += &part;
                } else if file.ends_with(FA) {
                    let part = format!("--fasta {} ", file.display());
                    tai += &part;
                    blast += &part;
                } else if file.ends_with(INDEX) {
                    let part = format!("--index {} ", file.display());
                    tai += &part;
                    blast += &part;
                } else {
                    continue;
                };
            }

            jobs.push(Job::from(blast));
            jobs.push(Job::from(tai));
        }
    }
}

/// Processes a single BED file in "canonical" mode by chunking it and creating
/// parallel jobs for BLAST and TAI (Translation AI) analysis.
///
/// This function is a core part of the canonical processing pipeline. It first validates
/// the input BED file to ensure it exists and is not empty. It then extracts a suffix
/// from the filename, which is used for organizing output files.
///
/// The function's main task is to call `bounded_extract` to split the large BED file
/// into smaller, more manageable chunks. This is the second chunking step in the
/// pipeline, specifically designed for handling fusion-related files. For each
/// resulting chunk (a pair of FASTA and BED files), it generates two command-line
/// jobs:
/// 1. A **BLAST job**: Executes the `ORF_RELEASE` tool with the `blast` subcommand.
///    This job searches for open reading frames (ORFs) and aligns them against a
///    specified protein database.
/// 2. A **TAI job**: Executes the `ORF_RELEASE` tool with the `tai` subcommand.
///    This job performs a Translation AI analysis on the ORFs.
///
/// Both of these jobs are added to a `jobs` vector to be executed in parallel
/// by a `ParallelExecutor`.
///
/// # Arguments
///
/// * `bed` - A `PathBuf` to the input BED file to be processed.
/// * `twobit` - A reference to a `PathBuf` pointing to the `.2bit` genome file.
/// * `step_output_dir` - A reference to a `PathBuf` representing the base output directory for this step.
/// * `chunk_size` - The maximum number of records to include in each chunk.
/// * `args` - A reference to a `Vec<String>` containing command-line arguments for the child processes,
///            specifically `orfipy`'s e-value, minimum ORF length, and the path to the protein database.
/// * `jobs` - A mutable reference to a `Vec<Job>` where the generated BLAST and TAI jobs are stored.
///
/// # Panics
///
/// This function will panic if:
/// - It cannot get the file stem from the `bed` path.
/// - It cannot extract the file suffix from the basename.
/// - It cannot determine the parent directory for a chunked BED file.
/// - The `bounded_extract` or `Job::from` functions fail.
fn cannonical_processing(
    bed: PathBuf,
    twobit: &PathBuf,
    step_output_dir: &PathBuf,
    chunk_size: usize,
    args: &Vec<String>,
    jobs: &mut Vec<Job>,
) {
    if !bed.exists() || std::fs::metadata(&bed).unwrap().len() == 0 {
        log::warn!("WARN: {} does not exist or its empty!", bed.display());
        return;
    }

    let basename = bed
        .file_stem()
        .unwrap_or_else(|| panic!("ERROR: could not get file stem from -> {:?}", bed))
        .to_string_lossy();
    let suffix = basename
        .split(".")
        .last()
        .unwrap_or_else(|| panic!("ERROR: could not get suffix from file -> {:?}", bed));

    // INFO: inflection point -> chunking fusion files [2nd chunking step in the pipeline]
    let paths = bounded_extract(&bed, &twobit, step_output_dir, chunk_size, suffix);

    for (chunked_fa, chunked_bed) in paths {
        let chunked_dir = &chunked_bed.parent().unwrap_or_else(|| {
            panic!(
                "ERROR: could not get parent directory for chunked_bed -> {}",
                chunked_bed.display()
            )
        });

        let blast = format!(
            "{} blast -e {} --fasta {} --alignments {} --outdir {} --orf-min-len {} --db {}",
            ORF_RELEASE,
            args[1], // INFO: orfipy
            &chunked_fa.display(),
            &chunked_bed.display(),
            chunked_dir.display(),
            args[2], // INFO: orf_min_len,
            args[3], // INFO: database
        );

        let tai = format!(
            "{} tai --fasta {} --alignments {} --outdir {}",
            ORF_RELEASE,
            chunked_fa.display(),
            chunked_bed.display(),
            chunked_dir.display(),
        );

        jobs.push(Job::from(blast));
        jobs.push(Job::from(tai));
    }
}

/// Extract sequences for fusion-free predicted
/// transcripts from a .2bit file into a .fa file
/// chunked
///
/// # Arguments
///
/// * `reads` - Path to the reads file
/// * `twobit` - Path to the .2bit file
/// * `step_output_dir` - Path to the output directory
///
/// # Returns
///
/// * None
///
/// # Example
///
/// ```rust
/// let jobs = orf(
///     &step,
///     &config,
///     &input_dir,
///     &step_output_dir,
/// );
/// ```
pub fn bounded_extract(
    reads: &PathBuf,
    twobit: &PathBuf,
    step_output_dir: &PathBuf,
    chunk_size: usize,
    suffix: &str,
) -> Vec<(PathBuf, PathBuf)> {
    log::info!(
        "INFO: Extracting mapped read sequences [{}] from .2bit file...",
        reads.display()
    );

    let tmp_dir = step_output_dir.join(format!("{}_{}", SEQS, suffix));
    std::fs::create_dir_all(&tmp_dir).unwrap_or_else(|e| {
        panic!(
            "ERROR: could not creat temporary directory in {} -> {e}",
            &tmp_dir.display()
        )
    });

    let bed = unpack::<GenePred, _>(vec![reads.clone()], OverlapType::Exon, false)
        .unwrap_or_else(|e| panic!("ERROR: could not unpack reads -> {}. {e}", reads.display()));
    let (genome, _) = get_sequences(twobit.clone()).unwrap_or_else(|| {
        panic!(
            "ERROR: could not get sequences from .2bit -> {}",
            twobit.display()
        )
    });

    // INFO: define the chunk size for parallel processing
    // INFO: if chunk size > bed records -> symlink
    let paths: Vec<(PathBuf, PathBuf)> = bed
        .into_par_iter()
        .flat_map(|(chrom, records)| {
            records
                .chunks(chunk_size)
                .enumerate()
                .map(move |(chunk_id, chunk)| (format!("{}:{}", chrom, chunk_id), chunk.to_vec()))
                .collect::<Vec<_>>()
        })
        .map(|(chunk_id, transcripts)| {
            let chunk_path = tmp_dir.join(&chunk_id);
            std::fs::create_dir_all(&chunk_path).unwrap();

            let chunk_fa = chunk_path.join(format!("tmp_chunk_{}.fa", &chunk_id));
            let chunk_bed = chunk_path.join(format!("tmp_chunk_{}.bed", &chunk_id));

            let mut writer_fa = BufWriter::new(File::create(&chunk_fa).unwrap());
            let mut writer_bed = BufWriter::new(File::create(&chunk_bed).unwrap());

            let chr = chunk_id.split(':').next().unwrap_or(&chunk_id);
            for tx in transcripts {
                let seq = match tx.strand {
                    Strand::Forward => Sequence::new(
                        genome.get(chr).expect("ERROR: missing chromosome")
                            [tx.start as usize..tx.end as usize]
                            .as_ref(),
                    ),
                    Strand::Reverse => Sequence::new(
                        genome.get(chr).expect("ERROR: missing chromosome")
                            [(SCALE - tx.end) as usize..(SCALE - tx.start) as usize]
                            .as_ref(),
                    )
                    .reverse_complement(),
                };

                writeln!(writer_fa, ">{}\n{}", tx.name, seq).unwrap();
                writeln!(writer_bed, "{}", tx.line).unwrap();
            }

            (chunk_fa, chunk_bed)
        })
        .collect();

    return paths;
}

/// Merges TOGA predictions into a single file by invoking the `orf toga` command.
///
/// This function constructs and executes a shell command to run the `orf toga`
/// subcommand. It retrieves the path to the TOGA results from the provided
/// `Config` and specifies the output directory for the merged file.
///
/// # Arguments
///
/// * `step_output_dir` - A `PathBuf` representing the output directory for the current pipeline step.
/// * `config` - A reference to a `Config` struct, used to retrieve the path to TOGA results.
/// * `step` - A reference to a `PipelineStep` enum, used to identify the current step
///            and retrieve its custom fields from the `config`.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to retrieve the TOGA path from the `config`.
/// - The `shell` function (which executes the command) encounters an error.
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
///
/// struct Config;
/// impl Config {
///     fn get_step_custom_field(&self, step: &PipelineStep, field_name: &str) -> PathBuf {
///         // In a real scenario, this would return a path based on config and step
///         PathBuf::from("/path/to/toga_raw_results")
///     }
/// }
///
/// enum PipelineStep {
///     MergeToga,
/// }
///
/// fn shell(cmd: String, msg: &str, tool: &str) {
///     println!("Executing: {}", cmd);
///     println!("Message: {}", msg);
///     println!("Tool: {}", tool);
/// }
///
/// const TOGA: &str = "toga_field_name"; // Placeholder for the actual field name in Config
/// const ORF_RELEASE: &str = "orf"; // Placeholder for the actual orf executable path
///
/// let output_dir = PathBuf::from("/tmp/pipeline_output/merge_toga_step");
/// let app_config = Config;
/// let current_step = PipelineStep::MergeToga;
///
/// // This would execute a command similar to:
/// // orf --path /path/to/toga_raw_results --outdir /tmp/pipeline_output/merge_toga_step
/// __merge_toga(&output_dir, &app_config, &current_step);
/// ```
fn __merge_toga(step_output_dir: &PathBuf, config: &Config, step: &PipelineStep) {
    let toga = config.get_step_custom_field(step, TOGA);
    let msg = "INFO: Merging TOGA predictions in a single file...";
    let tool = "iso-orf";

    let cmd = format!(
        "{} toga --path {} --outdir {}",
        ORF_RELEASE,
        toga,
        step_output_dir.display()
    );

    shell(cmd, msg, tool);
}

/// Classifies a file path based on its final directory or filename suffix.
///
/// This function is a utility for pipeline steps that process files organized
/// into directories ending with either "accept" or "reject". It determines
/// the type of file by checking its suffix and returns a corresponding
/// static string (`FREE` or `FUSIONS`). This helps in dynamically naming
/// output files or directing logic based on the input type.
///
/// The function operates on the final component of the path. If the name
/// ends with "accept", it returns the string associated with free reads (`FREE`).
/// If it ends with "reject", it returns the string for fusions (`FUSIONS`).
/// If the suffix is not recognized, the program will panic.
///
/// # Arguments
///
/// * `path` - A reference to a `PathBuf` pointing to a file or directory.
///
/// # Returns
///
/// A `&'static str` which will be either `FREE` or `FUSIONS`.
///
/// # Panics
///
/// This function will panic if:
/// - The path has no file name component.
/// - The file name cannot be converted to a UTF-8 string.
/// - The file name's suffix does not end with "accept" or "reject".
///
/// # Example
///
/// ```rust
/// use std::path::{Path, PathBuf};
/// use std::process;
///
/// // Assume FREE and FUSIONS are defined constants somewhere in the module.
/// const FREE: &str = "free";
/// const FUSIONS: &str = "fusions";
///
/// // Dummy function to demonstrate how `get_suffix_from_path` would be used.
/// fn get_suffix_from_path(path: &Path) -> &'static str {
///     path.file_name()
///         .and_then(|name| name.to_str())
///         .map(|name_str| {
///             if name_str.ends_with("accept") {
///                 FREE
///             } else if name_str.ends_with("reject") {
///                 FUSIONS
///             } else {
///                 eprintln!("ERROR: directory suffix not recognized: {:?}", path);
///                 std::process::exit(1);
///             }
///         }).unwrap_or_else(|| panic!("ERROR: could recognize suffix!"))
/// }
///
/// let path = Path::new("/foo/bar/chrX_all.aligned.reject");
///
/// match get_suffix_from_path(path) {
///     "free" => println!("This is a path for free reads."),
///     "fusions" => println!("This is a path for fusion reads."),
///     _ => unreachable!(), // This case is handled by the panic in the function
/// }
///
/// let accept_path = Path::new("some/directory/accept");
/// let result = get_suffix_from_path(accept_path);
/// assert_eq!(result, FREE);
///
/// let reject_path = Path::new("another/dir/reject");
/// let result = get_suffix_from_path(reject_path);
/// assert_eq!(result, FUSIONS);
/// ```
pub fn get_suffix_from_path(path: &PathBuf) -> &'static str {
    path.file_name()
        .and_then(|name| name.to_str())
        .map(|name_str| {
            if name_str.ends_with("accept") {
                FREE
            } else if name_str.ends_with("reject") {
                FUSIONS
            } else {
                log::error!("ERROR: subdir {path:?} suffix could not be recognized -> should end it 'accept' or 'reject'!");
                std::process::exit(1);
            }
        }).unwrap_or_else(|| panic!("ERROR: could recognize suffix!"))
}
