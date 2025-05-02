use crate::{config::*, consts::*, executor::job::Job};
use std::path::PathBuf;

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
    let args = config.get_step_args(step, vec![INPUT_DIR, OUTPUT_DIR, MEMORY, TIME, GENOME]);
    let fields = config.get_step_custom_fields(step, vec![GENOME]);

    let reads = input_dir.join(CLUSTERED_FA);
    let alignment = step_output_dir.join(CU_ALN_SAM);

    let jobs = vec![Job::new()
        .task(*step)
        .arg(&args)
        .arg(&format!("-o {}", alignment.display()))
        .arg(&fields[0])
        .arg(reads.display())];

    log::info!("INFO [STEP 5]: Pre-processing completed -> Running...");

    return jobs;
}
