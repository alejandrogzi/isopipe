use crate::{config::*, consts::*, executor::job::Job};
use std::path::PathBuf;

/// Run ccs
///
/// # Arguments
/// * `step` - The pipeline step being processed.
/// * `config` - The configuration settings for the pipeline.
/// * `input_dir` - The directory containing input files.
/// * `step_output_dir` - The directory where output files will be written.
/// * `prefix` - The prefix to be used for output files.
///
/// # Returns
/// A vector of jobs to be executed.
///
/// # Example
/// ```
/// let jobs = __ccs(
///     &step,
///     &config,
///     &input_dir,
///     &step_output_dir,
///     "sample".to_string(),
/// );
/// ```
pub fn ccs(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
    prefix: String,
) -> Vec<Job> {
    let mut jobs = Vec::new();

    let fields = config.get_step_custom_fields(step, vec![CHUNK, REPORT]);
    let args = config.get_step_args(
        step,
        vec![INPUT_DIR, PREFIX, OUTPUT_DIR, CHUNK, MEMORY, TIME, REPORT],
    );

    // WARN: ignoring prefix + .subreads ending -> forcing to isolate samples
    for (idx, entry) in std::fs::read_dir(input_dir)
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
        .enumerate()
    {
        let idx = idx + 1; // WARN: 1-indexed
        let bam = entry.path();
        let out_bam = step_output_dir.join(format!("{}.ccs.{}.bam", prefix, idx));

        // INFO: order/index known ahead from get_step_custom_fields
        let chunks = format!("--chunk {}/{}", idx, fields[0]);
        let report = format!(
            "--report-file {}/{}_{}.txt",
            step_output_dir.display(),
            fields[1],
            idx
        );

        // WARN: need to check if bam has a .pbi file -> if not, run pbindex
        let mut pbi = bam.clone();
        pbi.set_extension("bam.pbi");
        if !pbi.exists() {
            log::warn!(
                "WARN: pbi file not found for {}, generating index...",
                bam.display()
            );

            generate_pb_index(&bam, &config);
        }

        let job = Job::new()
            .task(PipelineStep::Ccs)
            .arg(bam.to_str().expect("ERROR: failed to convert path to str"))
            .arg(
                out_bam
                    .to_str()
                    .expect("ERROR: failed to convert path to str"),
            )
            .arg(&chunks)
            .arg(&args)
            .arg(&report);

        jobs.push(job)
    }

    log::info!("INFO [STEP 1]: Pre-processing completed -> Running...");

    return jobs;
}

/// Generates a .pbi for a BAM file.
///
/// # Arguments
///
/// * `bam` - The path to the BAM file.
/// * `config` - The configuration for the pipeline.
///
/// # Example
///
/// ```rust, no_run
/// generate_pb_index(PathBuf::from("example.bam"), &config);
/// ```
fn generate_pb_index(bam: &PathBuf, config: &Config) {
    let msg = format!("Generating PBINDEX for {}", bam.display());
    let package = config
        .packages
        .get(PBINDEX)
        .expect("PBINDEX package not found");

    let cmd = format!(
        "module load {}/{} && pbindex {}",
        PBINDEX,
        package,
        bam.display()
    );

    shell(cmd, &msg, PBINDEX);
}
