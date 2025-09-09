use crate::{config::*, consts::*, executor::job::Job};
use std::path::PathBuf;

pub fn load(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
    executor: &mut crate::core::ParallelExecutor,
) -> Vec<Job> {
    let mut jobs = Vec::new();
    let chrom_sizes = config.get_step_custom_field(step, CHROM_SIZES);
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

        let cmd = format!(
            "{BED_TO_BIG_BED} -tab -sort -as={SCHEMA} -type=bed12+25 {} {} {}",
            input.display(),
            chrom_sizes,
            step_output_dir.join("bb").join(bb).display()
        );

        jobs.push(Job::from(cmd));
    }

    jobs
}
