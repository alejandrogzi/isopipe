use crate::{config::*, consts::*, executor::job::Job};
use std::path::PathBuf;

/// Run isoseq refine
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
/// ```
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
    prefix: String,
) -> Vec<Job> {
    let mut jobs = Vec::new();

    let args = config.get_step_args(step, vec![INPUT_DIR, PREFIX, OUTPUT_DIR, MEMORY, TIME]);
    let fields = config.get_step_custom_fields(step, vec![PRIMERS]);

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
            .collect::<Vec<&str>>();
        let identifier = identifier.get(2).expect("ERROR: failed to get identifier");
        let out_bam = step_output_dir.join(format!("{}.flnc.{}.bam", prefix, identifier));

        let job = Job::new()
            .task(PipelineStep::Refine)
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

pub fn cluster(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
    prefix: String,
) -> Vec<Job> {
    // isoseq3
    // cluster
    // ${P_out_isoC}/ALL.flnc.fofn ${P_out_isoC}/ALL.CuP.bam
    // --singletons
    // -verbose
    // --split-bam $splitBAM
    // --num-threads $defVars{'nThreadsIsoSeq3'}
    // --log-file ${P_out_isoC}/ALL.CuP.log

    todo!()

    // let args = config
    //     .params()
    //     .get(&PipelineStep::Cluster)
    //     .expect("ERROR: ccs not found in config.toml!")
    //     .flat(Some(vec![INPUT_DIR, PREFIX, OUTPUT_DIR, MEMORY, TIME]));

    // let mut file_count = 0;

    // for entry in std::fs::read_dir(input_dir)
    //     .expect("Failed to read assets directory")
    //     .flatten()
    //     .filter(|entry| {
    //         entry
    //             .path()
    //             .extension()
    //             .expect("ERROR: no extension found!")
    //             == BAM
    //     })
    // {
    //     file_count += 1;

    //     let job = Job::new()
    //         .task(PipelineStep::Cluster)
    //         .arg(&args)
    //         .arg(&output_dir);

    //     jobs.push(job);

    //     if file_count > 1 {
    //         log::error!(
    //             "ERROR: more than one .sam file found in input_dir. This is a bug!"
    //         );
    //         std::process::exit(1);
    //     }
    // }
}
