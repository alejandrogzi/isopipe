use crate::{
    config::*,
    consts::*,
    core::pbindex,
    executor::{job::Job, manager::ParallelExecutor},
};
use std::path::PathBuf;

/// Merge BAM files in a directory using samtools
/// and index the merged BAMs
///
/// # Arguments
/// * `input_dir` - The directory containing the BAM files to merge.
/// * `executor` - The executor to use for running the merge command.
/// * `config` - The configuration to use for running the merge command.
///
/// # Example
/// ```
/// use std::path::PathBuf;
/// use isopipe::core::samtools;
/// use isopipe::executor::manager::ParallelExecutor;
/// use isopipe::config::Config;
///
/// let input_dir = PathBuf::from("/path/to/bam/files");
/// let mut executor = ParallelExecutor::new();
/// let config = Config::default();
///
/// samtools::merge(&input_dir, &mut executor, &config);
/// ```
pub fn merge(input_dir: &PathBuf, executor: &mut ParallelExecutor, config: &Config) {
    const THREADS: u32 = 8;
    const MEMORY: u32 = 8;

    let mut jobs = Vec::new();
    let mut pbi = Vec::new();
    let package = config.get_custom_package(SAMTOOLS);

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
        // INFO: format of chunked bams: {prefix}.{name}.ccs.{chunk}.bam
        // INFO: basename -> {prefix}.{name}.ccs
        let chunk = PathBuf::from(
            entry
                .path()
                .file_stem()
                .expect(&format!("Failed to get file name for {:?}", entry.path())),
        );
        let basename = chunk
            .file_stem()
            .expect(&format!("Failed to get file name for {:?}", entry.path()))
            .to_string_lossy();

        // INFO: format of merged: {prefix}.{name}.ccs.merged.bam
        let merged = input_dir.join(format!("{}.{}", basename, MERGED_BAM));

        if !merged.exists() {
            // INFO: format of wildcard: {prefix}.{name}.ccs.*.bam
            let wildcard = input_dir.join(format!("{}.*.{}", basename, BAM));

            // INFO: if merged file does not existe, we merge and delete the unmerged files
            let cmd = format!(
                "samtools merge -@{} {} {} && rm {}",
                THREADS,
                wildcard.display(),
                merged.display(),
                wildcard.display(),
            );

            pbi.push(merged);

            let job = Job::from(cmd);
            jobs.push(job);
        }
    }

    executor.add_jobs(jobs).and_send(
        config,
        SAMTOOLS,
        input_dir.clone(),
        THREADS,
        MEMORY,
        package,
    );

    pbindex::pbindex(pbi, config, executor, input_dir);
}
