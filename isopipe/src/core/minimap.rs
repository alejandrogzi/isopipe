use crate::{config::*, consts::*, executor::job::Job};

use std::{
    fs::{create_dir_all, File},
    io::BufWriter,
    path::{Path, PathBuf},
};
use twobit::{convert, TwoBitFile};

/// Run minimap2
///
/// # Arguments
/// * `step` - The pipeline step to run
/// * `config` - The configuration for the pipeline
/// * `input_dir` - The directory containing the input files
/// * `step_output_dir` - The directory to write the output files to
///
///
/// # Note
///
/// 'input_dir' points to CHUNKS, the output of the chunking part
/// considering only hg and singletons
///
/// # Returns
/// A vector of jobs to run
///
/// # Example
/// ```
/// let jobs = minimap2(
///     &step,
///     &config,
///     &input_dir,
///     &step_output_dir,
/// );
/// ```
pub fn minimap2(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &Path,
    executor: &mut crate::core::ParallelExecutor,
) -> Vec<Job> {
    let mut jobs = Vec::new();

    let args = config.get_step_args(
        step,
        vec![
            INPUT_DIR, OUTPUT_DIR, MEMORY, TIME, GENOME, CHUNK, KEEP_TEMP,
        ],
    );

    let keep_temp = config
        .get_step_custom_field(step, KEEP_TEMP)
        .parse::<bool>()
        .unwrap_or(false);

    if keep_temp {
        log::info!("INFO [STEP 5]: Keeping temporary files -> keep_temp set to true!");
    }

    let genome_index = index_genome(step, config, step_output_dir, executor);

    let _ = create_dir_all(step_output_dir.join(SAM));
    let _ = create_dir_all(step_output_dir.join(BAM));

    for entry in std::fs::read_dir(input_dir)
        .expect("ERROR: failed to read input directory")
        .flatten()
        .map(|entry| entry.path())
        .filter(|path| {
            path.file_name()
                .and_then(|name| name.to_str())
                .map(|filename| GZ_EXTENSIONS.iter().any(|ext| filename.ends_with(ext)))
                .unwrap_or(false)
        })
    {
        let basename = entry.with_extension(""); // INFO: discard .gz

        // INFO: send anything .f[a/q] inside chunks/
        let sam = &step_output_dir
            .join(SAM)
            .join(basename.with_extension(SAM).file_name().unwrap());
        let bam = &step_output_dir
            .join(BAM)
            .join(basename.with_extension(BAM).file_name().unwrap());

        let compression = if !keep_temp {
            format!(
                "&& {SAMTOOLS} view -@ 8 -b {} | {SAMTOOLS} sort -@ 8 -o {} && rm {} && {SAMTOOLS} index -@ 8 {}",
                sam.display(),
                bam.display(),
                sam.display(),
                bam.display()
            )
        } else {
            format!(
                "&& {SAMTOOLS} view -@ 8 -b {} | {SAMTOOLS} sort -@ 8 -o {} && {SAMTOOLS} index -@ 8 {}",
                sam.display(),
                bam.display(),
                bam.display()
            )
        };

        let job = Job::new()
            .task(*step)
            .arg(&args)
            .arg(format!("-o {}", sam.display()))
            .arg(genome_index.display())
            .arg(entry.display())
            .arg(compression);

        jobs.push(job);
    }

    if jobs.is_empty() {
        log::warn!(
            "WARN: No cannonical jobs found to run for minimap2 in {} -> trying to grab any .fasta.gz!",
            input_dir.display()
        );
        build_non_cannonical(
            input_dir,
            step_output_dir,
            step,
            args,
            genome_index,
            &mut jobs,
            keep_temp,
        );
    }

    log::info!("INFO [STEP 5]: Pre-processing completed -> Running...");

    jobs
}

/// Creates a FASTA file from a 2bit file in-place
///
/// # Arguments
/// * `genome` - The path to the 2bit file
///
/// # Returns
/// The path to the created FASTA file
///
/// # Example
/// ```
/// let genome = "path/to/genome.2bit";
/// let fasta = twobit_to_fa(genome);
/// assert_eq!(fasta, "path/to/genome.fa");
/// ```
fn twobit_to_fa(genome: String, step_output_dir: &Path) -> String {
    let mut twobit = TwoBitFile::open(&genome)
        .unwrap_or_else(|e| panic!("ERROR: Failed to open 2bit file -> {} -> {e}", genome));
    let fasta = step_output_dir.join(GENOME_FA);

    let mut writer = BufWriter::new(File::create(&fasta).unwrap_or_else(|e| {
        panic!(
            "ERROR: Failed to create FASTA file -> {} -> {e}",
            fasta.display()
        )
    }));

    convert::fasta::to_fasta(&mut twobit, &mut writer)
        .unwrap_or_else(|e| panic!("ERROR: Failed to convert {} to FASTA -> {e}", genome));

    fasta.display().to_string()
}

/// Returns the path to the genome file
///
/// # Arguments
/// * `config` - The configuration object
/// * `step` - The pipeline step
///
/// # Returns
/// The path to the genome file
///
/// # Example
/// ```
/// let config = Config::new();
/// let step = PipelineStep::new();
/// let genome = get_genome(&config, &step);
/// assert_eq!(genome, "path/to/genome.fa");
/// ```
fn get_genome(config: &Config, step: &PipelineStep, step_output_dir: &Path) -> String {
    let fields = config.get_step_custom_fields(step, vec![GENOME]);
    let file = fields
        .first()
        .unwrap_or_else(|| panic!("ERROR: {} field not found!", GENOME));

    if file.ends_with(TWOBIT) {
        twobit_to_fa(file.clone(), step_output_dir)
    } else {
        file.to_string()
    }
}

/// Builds a list of non-canonical files
///
/// # Arguments
///
/// * `input_dir` - The input directory
///
/// # Returns
///
/// A vector of paths to the non-canonical files
///
/// # Example
///
/// ```rust, no_run
/// let input_dir = PathBuf::from("path/to/input");
/// let files = __build_non_cannonical(&input_dir);
///
/// assert_eq!(files.len(), 2);
/// ```
fn build_non_cannonical(
    input_dir: &PathBuf,
    step_output_dir: &Path,
    step: &PipelineStep,
    args: String,
    genome: PathBuf,
    jobs: &mut Vec<Job>,
    keep_temp: bool,
) {
    // INFO: only considering .fasta.gz + fa.gz; others are already considered as cannonical!
    for entry in std::fs::read_dir(input_dir)
        .expect("ERROR: failed to read input directory!")
        .flatten()
        .map(|entry| entry.path())
        .filter(|path| {
            path.file_name()
                .and_then(|name| name.to_str())
                .map(|filename| filename.ends_with(FASTA_GZ) || filename.ends_with(FA_GZ))
                .unwrap_or(false)
        })
    {
        log::info!("INFO: Found non-cannonical file {}", entry.display());

        let bind = entry.with_extension("");
        let basename = bind
            .file_stem()
            .unwrap_or_else(|| panic!("ERROR: Failed to get file stem from {}", entry.display()))
            .to_string_lossy();

        // INFO: each one to sam/ or bam/
        let sam = step_output_dir
            .join(SAM)
            .join(format!("aligned.{}.{}", basename, SAM));
        let bam = step_output_dir
            .join(BAM)
            .join(format!("aligned.{}.{}", basename, BAM));

        let compression = if !keep_temp {
            format!(
                "&& {SAMTOOLS} view -@ 8 -b {} | {SAMTOOLS} sort -@ 8 -o {} && rm {} && {SAMTOOLS} index -@ 8 {}",
                sam.display(),
                bam.display(),
                sam.display(),
                bam.display()
            )
        } else {
            format!(
                "&& {SAMTOOLS} view -@ 8 -b {} | {SAMTOOLS} sort -@ 8 -o {} && {SAMTOOLS} index -@ 8 {}",
                sam.display(),
                bam.display(),
                bam.display()
            )
        };

        let job = Job::new()
            .task(*step)
            .arg(&args)
            .arg(format!("-o {}", sam.display()))
            .arg(genome.display())
            .arg(entry.display())
            .arg(compression);

        jobs.push(job);
    }

    if jobs.is_empty() {
        log::error!(
            "ERROR: No non-cannonical files found in {}!",
            input_dir.display()
        );
        std::process::exit(1);
    } else {
        log::info!("INFO: sending {} non-cannonical jobs...", jobs.len());
    }
}

/// Indexes a genome file using minimap2
///
/// This function first retrieves the genome file from the configuration and then builds an index for it using the `minimap2` tool. The resulting index is saved in the step's output directory.
///
/// # Arguments
///
/// * `step` - The current pipeline step
/// * `config` - The pipeline configuration
/// * `step_output_dir` - The output directory for the current step
/// * `executor` - The parallel executor for running jobs
///
/// # Returns
///
/// The path to the newly created genome index file
///
/// # Example
///
/// ```rust, no_run
/// let mut executor = ParallelExecutor::new();
/// let step = PipelineStep::IndexGenome;
/// let output_dir = PathBuf::from("path/to/output");
/// let genome_index_path = index_genome(&step, &config, &output_dir, &mut executor);
///
/// assert!(genome_index_path.exists());
/// ```
pub fn index_genome(
    step: &PipelineStep,
    config: &Config,
    step_output_dir: &Path,
    executor: &mut crate::core::ParallelExecutor,
) -> PathBuf {
    // INFO: genome -> step_output_dir.join("genome.fa")
    let genome = get_genome(config, step, step_output_dir);
    let index = PathBuf::from(genome.clone()).with_extension("index");

    // INFO: run minimap indexing and return new index path
    let cmd = format!("{MINIMAP2} -d {} {genome}", index.display());
    let job = Job::from(cmd);

    let package = config.get_package_from_step(step);

    executor
        .add_job(job)
        // .execute(config, step, step_output_dir.clone(), Some("minimizer"));
        .and_send(
            config,
            &step.to_unique_str(),
            step_output_dir.to_path_buf(),
            8,
            20,
            Some(package),
            Some("minimizer"),
        );

    index
}
