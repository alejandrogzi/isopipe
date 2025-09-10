use crate::{config::*, consts::*, executor::job::Job};
use std::path::PathBuf;

/// Load and concatenate decision files, convert to bigBed format, and optionally upload
///
/// # Arguments
/// * `step` - The pipeline step to run
/// * `config` - The configuration for the pipeline
/// * `input_dir` - The directory containing the input files (from polish step)
/// * `step_output_dir` - The directory to write the output files to
/// * `executor` - Parallel executor for running concatenation jobs
///
/// # Note
///
/// 'input_dir' points to the polish step output, containing chromosome-specific
/// decision directories with categorized BED files. The function concatenates
/// files across all chromosomes for each category, then converts them to bigBed
/// format for visualization.
///
/// Categories processed: pass.bed, trash.bed, rt_reads.bed, retention.bed,
/// intraprimming.bed, truncation.bed, nmd.bed, fusions.bed
///
/// # Returns
/// A vector of jobs to run for bigBed conversion and optional upload
///
/// # Example
/// ```rust, ignore
/// let jobs = load(
///     &step,
///     &config,
///     &input_dir,
///     &step_output_dir,
///     &mut executor,
/// );
/// ```
pub fn load(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
    executor: &mut crate::core::ParallelExecutor,
) -> Vec<Job> {
    let mut jobs = Vec::new();

    let chrom_sizes = config.get_step_custom_field(step, CHROM_SIZES);
    let server = config.get_step_custom_field(step, SERVER);
    let user = config.get_step_custom_field(step, USER);
    let target = config
        .get_step_custom_field(step, TARGET)
        .parse::<PathBuf>()
        .unwrap_or_else(|_| {
            panic!(
                "ERROR: could not parse target directory -> {}",
                config.get_step_custom_field(step, TARGET)
            )
        });
    let web = config
        .get_step_custom_field(step, WEB)
        .parse::<PathBuf>()
        .unwrap_or_else(|_| {
            panic!(
                "ERROR: could not parse web directory -> {}",
                config.get_step_custom_field(step, WEB)
            )
        });
    let upload = config
        .get_step_custom_field(step, UPLOAD_PUBLIC)
        .parse::<bool>()
        .unwrap_or(false);

    let package = config.get_package_from_step(step);

    // INFO: create bed and bb dirs
    for subdir in ["bed", "bb"] {
        let path = step_output_dir.join(subdir);
        std::fs::create_dir_all(&path).unwrap_or_else(|e| {
            panic!(
                "ERROR: could not create {} directory -> {:?}. {e}",
                subdir, path
            )
        });
    }

    let categories = [
        "pass.bed",
        "trash.bed",
        "rt_reads.bed",
        "retention.bed",
        "intraprimming.bed",
        "truncation.bed",
        "nmd.bed",
        "fusions.bed",
    ];

    // INFO: path would look like: {step_polish}/{chr}/decision/ -> looping for each chr
    let inner_jobs = categories
        .iter()
        .map(|target| {
            let pattern = input_dir.join("*").join("decision").join(target);
            let out = target; // INFO: output named same as input basename

            Job::from(format!(
                "cat {} >> {}",
                pattern.display(),
                step_output_dir.join("bed").join(out).display()
            ))
        })
        .collect();

    executor.add_jobs(inner_jobs).and_send(
        config,
        &step.to_unique_str(),
        step_output_dir.clone(),
        1,
        4,
        Some(package),
        Some("cat"),
    );

    // INFO: each .bed in step11_load should be converted into bigBed
    for bed in std::fs::read_dir(&step_output_dir.join("bed"))
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", step_output_dir))
        .flatten()
        .filter(|e| e.path().is_file())
    {
        let input = bed.path();

        let bind = input.with_extension("bb");
        let bb = bind
            .file_name()
            .unwrap_or_else(|| panic!("ERROR: could not get file name from {bed:?}"));

        let mut cmd = format!(
            "{BED_TO_BIG_BED} -tab -sort -as={SCHEMA} -type=bed12+25 {} {} {}",
            input.display(),
            chrom_sizes,
            step_output_dir.join("bb").join(bb).display()
        );

        if upload {
            cmd = format!(
                "{cmd} && ssh {user}@{server} mkdir {} && rsync -av {} {user}@{server}:{} && ln -sf {} {}",
                target.join(ISOPIPE).display(),
                step_output_dir.join("bb").join(bb).display(),
                target.join(ISOPIPE).join(bb).display(),
                target.join(ISOPIPE).join(bb).display(),
                web.join(bb).display(),
            );
        }

        log::debug!("DEBUG: executing cmd: {cmd:?}");
        jobs.push(Job::from(cmd));
    }

    jobs
}
