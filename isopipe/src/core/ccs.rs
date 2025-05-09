use crate::{
    config::*,
    consts::*,
    core::pbindex,
    executor::{job::Job, manager::ParallelExecutor},
};
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
    executor: &mut ParallelExecutor,
) -> Vec<Job> {
    let mut jobs = Vec::new();
    let mut require_pbi = Vec::new();

    let fields = config.get_step_custom_fields(step, vec![CHUNK, REPORT]);
    let args = config.get_step_args(
        step,
        vec![
            INPUT_DIR, PREFIX, OUTPUT_DIR, CHUNK, MEMORY, TIME, REPORT, NUM_CORES,
        ],
    );

    for (_, entry) in std::fs::read_dir(input_dir)
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
        let chunk_size = fields[0]
            .parse::<usize>()
            .expect("ERROR: Failed to parse chunk size");
        let bam = entry.path();

        for chunk_idx in 0..chunk_size {
            let chunk_idx = chunk_idx + 1;

            let out_bam = step_output_dir.join(format!(
                "{}.{}.ccs.{}.bam",
                prefix,
                bam.file_stem()
                    .expect(&format!(
                        "ERROR: failed to get name from bam: {}",
                        bam.display()
                    ))
                    .to_string_lossy(),
                chunk_idx
            ));

            let chunks = format!("--chunk {}/{}", chunk_idx, chunk_size);
            let report = format!(
                "--report-file {}/{}_{}.txt",
                step_output_dir.display(),
                fields[1],
                chunk_idx
            );

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

        // WARN: need to check if bam has a .pbi file -> if not, run pbindex
        let mut pbi = bam.clone();
        pbi.set_extension("bam.pbi");
        if !pbi.exists() {
            log::warn!(
                "WARN: pbi file not found for {}, generating index...",
                bam.display()
            );

            require_pbi.push(bam.clone());
        }
    }

    if !require_pbi.is_empty() {
        pbindex::pbindex(require_pbi, &config, executor, step_output_dir);
    }

    log::info!("INFO [STEP 1]: Pre-processing completed -> Running...");

    return jobs;
}
