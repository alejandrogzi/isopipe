use crate::executor::manager::ParallelExecutor;
use crate::{config::*, consts::*, executor::job::Job};

use config::{OverlapType, Sequence, Strand, FUSION_FREE, SCALE};
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
    global_output_dir: &PathBuf,
) -> Vec<Job> {
    let mut jobs = Vec::new();
    let args = config.get_step_custom_fields(step, vec![GENOME, ORF_MIN_LEN, DATABASE]);

    let mode = ParallelMode::from_str(&config.get_step_custom_field(step, PARALLEL_MODE));

    let twobit = PathBuf::from(args[0].clone());
    let chunk_size = crate::numerical!(config.get_step_custom_field(step, CHUNK) => usize)
        .unwrap_or_else(|e| panic!("ERROR: could not convert chunk to num -> {e}!"));

    log::info!(
        "INFO: Merging TOGA predictions in a single file here: {}...",
        step_output_dir.display()
    );

    let toga_merged = __merge_toga(step_output_dir, config, step);

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

    executor
        .add_jobs(jobs)
        .execute(config, step, global_output_dir.clone(), Some("prep"));

    predict(step_output_dir, &mode, &toga_merged)
}

/// Predicts Open Reading Frames (ORFs)
///
/// This function processes input data organized in a specific directory structure,
/// generating a list of `Job` commands that can be executed to perform the prediction.
/// The operation mode can be either `Chromosome` or `Genome`, dictating how the input
/// directories are traversed and processed.
///
/// # Arguments
///
/// * `step_output_dir` - A `PathBuf` reference pointing to the base directory where
///   the step output data (e.g., `seqs_{suffix}` or `toga` directories) is located.
///   This directory is expected to contain subdirectories for chromosomes and chunks.
/// * `mode` - A reference to a `ParallelMode` enum, specifying whether to run
///   the prediction per chromosome or per entire genome.
///   - `ParallelMode::Chromosome`: Processes data chunk by chunk within chromosome-specific directories.
///   - `ParallelMode::Genome`: (Currently unimplemented) Intended to process the entire genome data.
/// * `toga_merged` - A `PathBuf` reference to the merged TOGA (Tool for Genome Annotation) file.
///   This file is crucial for the prediction process and its existence is validated.
///
/// # Returns
///
/// A `Vec<Job>` containing a collection of `Job` structs, each representing a command
/// to be executed for ORF prediction. These commands are formatted to run an external
/// Python script (`PREDICT_PY`) with specific arguments derived from the input files.
///
/// # Panics
///
/// This function will panic under the following conditions:
///
/// * If `step_output_dir` or any of its expected subdirectories (chromosome subdirs, chunk subdirs)
///   cannot be read.
/// * If the `toga_merged` path does not exist.
/// * If required input files (`.bed` alignments, `blast` results, or `tai` results)
///   are not found within the expected `chunked_dir` for a `Chromosome` mode job.
fn predict(step_output_dir: &PathBuf, mode: &ParallelMode, toga_merged: &PathBuf) -> Vec<Job> {
    let mut jobs = Vec::new();

    match mode {
        ParallelMode::Chromosome => {
            log::info!("INFO: Running ORF prediction in parallel mode: Chromosome");

            // INFO: path would look like: {step_orf}/seqs_{suffix} or toga
            for entry in std::fs::read_dir(step_output_dir)
                .unwrap_or_else(|e| {
                    panic!("ERROR: could read directory -> {:?}. {e}", step_output_dir)
                })
                .flatten()
                .filter(|e| e.path().is_dir())
            {
                let subdir = entry.path(); // INFO: seqs_{suffix} or toga

                // INFO: structure {step_orf}/seqs_{suffix}/{chr}:{chunk}
                for chunk in std::fs::read_dir(&subdir)
                    .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", subdir))
                    .flatten()
                    .filter(|e| e.path().is_dir())
                {
                    let chunked_dir = chunk.path(); // INFO: {chr}:{chunk}

                    let mut alignments: Option<PathBuf> = None;
                    let mut blast: Option<PathBuf> = None;
                    let mut tai: Option<PathBuf> = None;

                    for file in std::fs::read_dir(&chunked_dir)
                        .unwrap_or_else(|e| {
                            panic!("ERROR: could read directory -> {:?}. {e}", chunked_dir)
                        })
                        .flatten()
                    {
                        let file = file.path();

                        if file.is_file() {
                            // INFO: matching .bed and not _reduced.bed
                            if let Some(ext) = file.extension().and_then(|e| e.to_str()) {
                                if ext == "bed"
                                    && !file
                                        .file_name()
                                        .unwrap_or_default()
                                        .to_string_lossy()
                                        .contains("_reduced")
                                {
                                    alignments = Some(file.clone());
                                }
                            }
                        } else if file.is_dir() {
                            let dir = file.file_name().unwrap_or_default().to_string_lossy();

                            // INFO: looks like -> {step_orf}/seqs_{suffix}/{chr}:{chunk}/{orf, tai}
                            if dir == ORF {
                                // INFO: grab .result
                                blast = Some(file.join("orfs.pep.dedup.dmd.result"));
                            } else if dir == TAI {
                                // INFO: grab .result
                                for tai_entry in std::fs::read_dir(&file)
                                    .unwrap_or_else(|e| {
                                        panic!("ERROR: could read directory -> {:?}. {e}", file)
                                    })
                                    .flatten()
                                {
                                    let res_path = tai_entry.path();
                                    if res_path.is_file() {
                                        tai = Some(res_path.clone());
                                    }
                                }
                            } else {
                                continue;
                            }
                        } else {
                            continue;
                        }
                    }

                    if !toga_merged.exists() {
                        panic!("ERROR: toga_merged path is empty!");
                    }

                    let cmd =
                        format!(
                        "source {} && {} --blast {} --tai {} --toga {} --alignments {} --outdir {} --prefix tmp_",
                        TAI_VENV, PREDICT_PY, blast.unwrap_or_else(|| {
                            panic!(
                                "ERROR: could not find blast file in chunked dir -> {}",
                                chunked_dir.display()
                            )
                        }).display(),
                        tai.unwrap_or_else(|| {
                            panic!(
                                "ERROR: could not find tai file in chunked dir -> {}",
                                chunked_dir.display()
                            )
                        }).display(),
                        toga_merged.display(),
                        alignments.unwrap_or_else(|| {
                            panic!(
                                "ERROR: could not find alignments file in chunked dir -> {}",
                                chunked_dir.display()
                            )
                        }).display(),
                        chunked_dir.display()
                    );

                    jobs.push(Job::from(cmd));
                }
            }
        }
        ParallelMode::Genome => {
            log::info!("INFO: Running ORF prediction in parallel mode: Genome");
            todo!()
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

    let suffixes = vec![config::FUSIONS, FUSION_FREE];
    let matched_suffixes = vec![crate::consts::FUSIONS, FREE];

    // INFO: structure should be {step_fusion}/chr{chr}_all.aligned.{accept/reject}
    for entry in std::fs::read_dir(input_dir)
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", input_dir))
        .flatten()
        .filter(|e| e.path().is_dir())
    {
        let subdir = entry.path();

        // INFO: structure should be {step_fusion}/chr{chr}_all.aligned.accept/fusions*
        // INFO: expected: fusions.bed / fusions.free.bed
        for bed in std::fs::read_dir(&subdir)
            .unwrap()
            .flatten()
            .filter(|e| e.path().is_file())
        {
            let bed = bed.path();
            let suffix =
                get_suffix_from_path(&bed, &suffixes, &matched_suffixes).unwrap_or_else(|| {
                    panic!(
                        "ERROR: path {:?} does not end with any recognized suffix: {:?}",
                        bed, matched_suffixes
                    )
                }); // INFO: safe to unwrap and error because we always expect matches only!

            // INFO: collect paths and collect extract cmds
            // INFO: end path would look like: {step_orf}/seqs_{suffix}/{chr}:{chunk}/{name}{fa/bed}
            let cmd = format!(
                "{} --twobit {} --bed {} -o {} --index --suffix {} --chunk-size {}",
                EXTRACT_RELEASE,
                &twobit.display(),
                bed.display(),
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
        .and_send(config, EXTRACT, step_output_dir.clone(), 8, 16, None, None);
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
    let suffixes = vec![REDUCED_BED, FA, INDEX];

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
                "{} blast --outdir {} --orf-min-len {} --db {} ",
                ORF_RELEASE,
                chunked_dir.display(),
                args[1], // INFO: orf_min_len,
                args[2], // INFO: database
            );

            let mut tai = format!("{} tai --outdir {} ", ORF_RELEASE, chunked_dir.display(),);

            for file in std::fs::read_dir(&chunked_dir)
                .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", chunked_dir))
                .flatten()
                .filter(|e| e.path().is_file())
            {
                let file = file.path();
                let ext = get_suffix_from_path(&file, &suffixes, &suffixes);

                if ext == Some(REDUCED_BED) {
                    let part = format!("--alignments {} ", file.display());
                    tai += &part;
                    blast += &part;
                } else if ext == Some(FA) {
                    let part = format!("--fasta {} ", file.display());
                    tai += &part;
                    blast += &part;
                } else if ext == Some(INDEX) {
                    let part = format!("--index {} ", file.display());
                    tai += &part;
                    blast += &part;
                } else {
                    continue;
                };
            }

            let cmd = format!("{} && {}", tai, blast);

            jobs.push(Job::from(cmd));
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

        let cmd = format!(
            "{} tai --fasta {} --alignments {} --outdir {} && {} blast -e {} --fasta {} --alignments {} --outdir {} --orf-min-len {} --db {}",
            ORF_RELEASE,
            chunked_fa.display(),
            chunked_bed.display(),
            chunked_dir.display(),
            ORF_RELEASE,
            args[1], // INFO: orfipy
            &chunked_fa.display(),
            &chunked_bed.display(),
            chunked_dir.display(),
            args[2], // INFO: orf_min_len,
            args[3], // INFO: database
        );

        jobs.push(Job::from(cmd));
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
fn __merge_toga(step_output_dir: &PathBuf, config: &Config, step: &PipelineStep) -> PathBuf {
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

    return step_output_dir.join(TOGA).join("toga_merged.tsv");
}

/// Retrieves a return value by matching the file's suffix against a predefined list.
///
/// This function provides a flexible way to classify a file path based on its name's
/// suffix. It takes a vector of possible suffixes to match and a corresponding vector
/// of return values. It iterates through these pairs and returns the value associated
/// with the first matching suffix found at the end of the file name. This is a robust
/// alternative to hard-coded checks, allowing for easy expansion of supported suffixes.
///
/// # Arguments
///
/// * `path` - A reference to a `PathBuf` representing the file or directory.
/// * `match_suffixes` - A vector of string slices (`&[&str]`) containing the suffixes to check against.
/// * `return_values` - A vector of string slices (`&[&'a str]`) containing the values to return for each corresponding suffix.
///
/// # Returns
///
/// A `&'a str` which is the value associated with the matched suffix.
///
/// # Panics
///
/// This function will panic if:
/// - The lengths of `match_suffixes` and `return_values` are not equal.
/// - The path has no file name component or the file name is not valid UTF-8.
/// - No matching suffix is found in the `match_suffixes` vector. The program
///   will exit with a non-zero status code in this case.
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
/// use std::process;
///
/// // Define suffixes and return values
/// const FREE: &str = "free";
/// const FUSIONS: &str = "fusions";
/// let match_suffixes = vec!["accept", "reject"];
/// let return_values = vec![FREE, FUSIONS];
///
/// // Example with a matching suffix
/// let path_accept = PathBuf::from("/foo/bar/chrX_all.aligned.accept");
/// let result = get_suffix_from_path(&path_accept, &match_suffixes, &return_values);
/// assert_eq!(result, FREE);
///
/// // Example with another matching suffix
/// let path_reject = PathBuf::from("another/dir/chrX_all.aligned.reject");
/// let result_reject = get_suffix_from_path(&path_reject, &match_suffixes, &return_values);
/// assert_eq!(result_reject, FUSIONS);
///
/// // Example that would cause a panic (mismatched lengths)
/// // let bad_return_values = vec![FREE];
/// // let _ = get_suffix_from_path(&path_accept, &match_suffixes, &bad_return_values);
///
/// // Example that would exit the program (no matching suffix)
/// // let path_unknown = PathBuf::from("some/file.txt");
/// // let _ = get_suffix_from_path(&path_unknown, &match_suffixes, &return_values);
/// ```
pub fn get_suffix_from_path<'a>(
    path: &PathBuf,
    match_suffixes: &Vec<&str>,
    return_values: &Vec<&'a str>,
) -> Option<&'a str> {
    assert_eq!(
        match_suffixes.len(),
        return_values.len(),
        "match_suffixes and return_values must have the same length"
    );

    let name_str = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or_else(|| {
            panic!(
                "Could not extract valid UTF-8 file name from path: {:?}",
                path
            )
        });

    for (suffix, value) in match_suffixes.iter().zip(return_values.iter()) {
        if name_str.ends_with(suffix) {
            return Some(value);
        }
    }

    return None;
}
