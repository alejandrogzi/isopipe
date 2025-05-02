use crate::{config::*, consts::*, executor::job::Job};
use std::path::PathBuf;

/// Run lima
///
/// # Arguments
/// * `step` - The pipeline step being processed.
/// * `config` - The configuration for the pipeline.
/// * `input_dir` - The directory containing the input BAM files.
/// * `step_output_dir` - The directory where the output files will be written.
/// * `prefix` - The prefix to be used for the output files.
///
/// # Returns
/// A vector of jobs to be executed.
///
/// # Example
/// ```
/// let jobs = __lima(
///     &step,
///     &config,
///     &input_dir,
///     &step_output_dir,
///     "prefix".to_string(),
/// );
/// ```
pub fn lima(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
    prefix: String,
) -> Vec<Job> {
    let mut jobs = Vec::new();

    let fields = config.get_step_custom_fields(step, vec![PRIMERS]);
    let args = config.get_step_args(
        step,
        vec![INPUT_DIR, PREFIX, OUTPUT_DIR, MEMORY, TIME, PRIMERS],
    );

    for entry in std::fs::read_dir(input_dir)
        .expect("Failed to read assets directory")
        .flatten()
        .filter(|entry| {
            entry
                .path()
                .extension()
                .and_then(|ext| ext.to_str())
                .map(|ext| ext.eq_ignore_ascii_case(BAM))
                .unwrap_or(false)
        })
    {
        let bam = entry.path();
        let identifier = bam
            .file_stem()
            .expect("ERROR: failed to get file stem")
            .to_str()
            .expect("ERROR: failed to convert path to str")
            .split('.')
            .last()
            .expect("ERROR: failed to get last element from .bam name!");

        let out_bam = step_output_dir.join(format!("{}.fl.{}.bam", prefix, identifier));

        let job = Job::new()
            .task(PipelineStep::Lima)
            .arg(&args)
            .arg(bam.to_str().expect("ERROR: failed to convert path to str"))
            .arg(&fields[0])
            .arg(
                out_bam
                    .to_str()
                    .expect("ERROR: failed to convert path to str"),
            );

        jobs.push(job)
    }

    log::info!("INFO [STEP 2]: Pre-processing completed -> Running...");

    return jobs;
}
