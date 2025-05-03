use crate::{config::*, consts::*, executor::job::Job};
use std::path::PathBuf;

/// Run isoseq3 refine
///
/// # Arguments
/// * `step` - The pipeline step being processed.
/// * `config` - The configuration object.
/// * `input_dir` - The directory containing the input BAM files.
/// * `step_output_dir` - The directory where the output files will be written.
/// * `prefix` - The prefix to be used for the output files.
///
/// # Returns
/// A vector of jobs to be executed.
///
/// # Example
/// ```rust, no_run
/// let jobs = __refine(
///     &step,
///     &config,
///     &input_dir,
///     &step_output_dir,
///     "prefix".to_string()
/// );
/// ```
pub fn refine(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
) -> Vec<Job> {
    let mut jobs = Vec::new();

    let args = config.get_step_args(step, vec![INPUT_DIR, OUTPUT_DIR, MEMORY, TIME, PRIMERS]);
    let fields = config.get_step_custom_fields(step, vec![PRIMERS]);

    // INFO: format of files: {prefix}.{name}.ccs.merged.fl.{primers}.bam
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
        let basename = bam
            .file_stem()
            .expect("ERROR: failed to get file stem")
            .to_string_lossy();

        let out_bam = step_output_dir.join(format!("{}.flnc.bam", basename));

        let job = Job::new()
            .task(*step)
            .arg(bam.to_str().expect("ERROR: failed to convert path to str"))
            .arg(&fields[0])
            .arg(
                out_bam
                    .to_str()
                    .expect("ERROR: failed to convert path to str"),
            )
            .arg(&args);

        jobs.push(job)
    }

    log::info!("INFO [STEP 3]: Pre-processing completed -> Running...");

    return jobs;
}

/// Run isoseq3 cluster
///
/// # Arguments
/// * `step` - The pipeline step being processed.
/// * `config` - The configuration for the pipeline.
/// * `input_dir` - The directory containing the input files.
/// * `step_output_dir` - The directory where the output files will be written.
/// * `prefix` - The prefix to use for the output files.
///
/// # Returns
/// A vector of jobs to be executed.
///
/// # Examples
/// ```rust, no_run
/// let jobs = cluster(
///     &step,
///     &config,
///     &input_dir,
///     &step_output_dir,
///     "prefix".to_string(),
/// );
/// ```
pub fn cluster(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
) -> Vec<Job> {
    let refine_fofn = format!("{}/*flnc.{}", input_dir.display(), BAM);
    let all_fofn = format!("{}/{}", step_output_dir.display(), FOFN);

    shell(
        format!("ls {} > {}", refine_fofn, all_fofn),
        "INFO: Grouping flnc reads...",
        CLUSTER,
    );

    let args = config.get_step_args(step, vec![INPUT_DIR, OUTPUT_DIR, MEMORY, TIME, LOG_FILE]);
    let out_bam = format!("{}/{}", step_output_dir.display(), CLUSTERED_BAM);
    let fields = config.get_step_custom_fields(step, vec![LOG_FILE]);

    let jobs = vec![Job::new()
        .task(*step)
        .arg(&all_fofn)
        .arg(&out_bam)
        .arg(&args)
        .arg(format!("--log-file {}/{}", &step_output_dir.display(), fields[0]).as_str())];

    log::info!("INFO [STEP 4]: Pre-processing completed -> Running...");

    return jobs;
}
