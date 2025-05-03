use crate::{
    config::*,
    consts::*,
    executor::{job::Job, manager::ParallelExecutor},
};
use std::path::PathBuf;

/// Generates a .pbi for a set of BAM files in parallel
///
/// # Arguments
///
/// * `bam` - The paths to the BAM files.
/// * `config` - The configuration for the pipeline.
/// * `executor` - The executor for parallel jobs.
/// * `step_output_dir` - The output directory for the step.
///
/// # Example
///
/// ```rust, no_run
/// pbindex(PathBuf::from("example.bam"), &config, &mut executor, &PathBuf::from("output"));
/// ```
pub fn pbindex(
    bams: Vec<PathBuf>,
    config: &Config,
    executor: &mut ParallelExecutor,
    step_output_dir: &PathBuf,
) {
    log::info!(
        "INFO [PBINDEX]: Generating .pbi indexes for {} BAM files...",
        bams.len()
    );

    let mut jobs = Vec::new();
    let package = config.get_custom_package(PBINDEX);

    bams.iter().for_each(|bam| {
        let cmd = format!("pbindex {}", bam.display());
        let job = Job::from(cmd);

        jobs.push(job);
    });

    executor
        .add_jobs(jobs)
        .and_send(config, PBINDEX, step_output_dir.clone(), 1, 8, package);
}
