use crate::{config::*, consts::*, executor::job::Job};

use std::{fs::File, io::BufWriter, path::PathBuf};
use twobit::{convert, TwoBitFile};

/// Run minimap2
///
/// # Arguments
/// * `step` - The pipeline step to run
/// * `config` - The configuration for the pipeline
/// * `input_dir` - The directory containing the input files
/// * `step_output_dir` - The directory to write the output files to
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
    step_output_dir: &PathBuf,
) -> Vec<Job> {
    let mut jobs = Vec::new();

    let args = config.get_step_args(step, vec![INPUT_DIR, OUTPUT_DIR, MEMORY, TIME, GENOME]);
    let genome = get_genome(config, step, step_output_dir);

    for category in CLUSTERING_CATEGORIES {
        let alignment = step_output_dir.join(format!("{}.{}.{}", CU_ALN, category, SAM));
        let reads = input_dir.join(format!("{}.{}.{}", CLUSTERED, category, FASTA_GZ));

        if !reads.exists()
            || reads
                .metadata()
                .expect(&format!(
                    "ERROR: failed to get metadata from {}",
                    reads.display()
                ))
                .len()
                == 0
        {
            continue;
        }

        let job = Job::new()
            .task(*step)
            .arg(&args)
            .arg(&format!("-o {}", alignment.display()))
            .arg(&genome)
            .arg(reads.display());

        jobs.push(job);
    }

    log::info!("INFO [STEP 5]: Pre-processing completed -> Running...");

    return jobs;
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
fn twobit_to_fa(genome: String, step_output_dir: &PathBuf) -> String {
    let mut twobit = TwoBitFile::open(&genome).expect("ERROR: Failed to open 2bit file");
    let fasta = step_output_dir.join(GENOME_FA);

    let mut writer =
        BufWriter::new(File::create(&fasta).expect("ERROR: Failed to create FASTA file"));

    let _ = convert::fasta::to_fasta(&mut twobit, &mut writer)
        .expect("ERROR: Failed to convert 2bit to FASTA");

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
fn get_genome(config: &Config, step: &PipelineStep, step_output_dir: &PathBuf) -> String {
    let fields = config.get_step_custom_fields(step, vec![GENOME]);
    let file = fields
        .get(0)
        .expect(format!("ERROR: {} field not found!", GENOME).as_str());

    let genome = if file.ends_with(TWOBIT) {
        twobit_to_fa(file.clone(), step_output_dir)
    } else {
        file.to_string()
    };

    genome
}
