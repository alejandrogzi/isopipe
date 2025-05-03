use crate::{
    config::*,
    consts::*,
    core::pbindex,
    executor::{job::Job, manager::ParallelExecutor},
};
use std::collections::HashMap;
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
    const THREADS: u32 = 16;
    const MEMORY: u32 = 8;

    let mut jobs = Vec::new();
    let mut pbi = Vec::new();

    let package = config.get_custom_package(SAMTOOLS);
    let groups = scan_groups(input_dir);

    for (group, bams) in groups {
        if bams.len() > 1 {
            let merged = input_dir.join(format!("{}.ccs.{}", group, MERGED_BAM));

            if !merged.exists() {
                // INFO: format of wildcard: {prefix}.{name}.ccs.*.bam
                let wildcard = input_dir.join(format!("{}*{}", group, BAM));

                // INFO: if merged file does not existe, we merge and delete the unmerged files
                let cmd = format!(
                    "samtools merge -@{} {} {} && rm {}",
                    THREADS,
                    merged.display(),
                    wildcard.display(),
                    wildcard.display(),
                );

                pbi.push(merged);

                let job = Job::from(cmd);
                jobs.push(job);
            }
        }
    }

    if jobs.is_empty() {
        return;
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

/// Scan groups in the input directory and return a HashMap of group
/// names to their corresponding BAM files.
///
/// # Arguments
/// * `input_dir` - The input directory to scan for BAM files.
///
/// # Returns
/// A HashMap of group names to their corresponding BAM files.
///
/// # Examples
/// ```
/// let input_dir = PathBuf::from("/path/to/input");
/// let groups = scan_groups(&input_dir);
/// assert_eq!(groups.len(), 2);
/// assert_eq!(groups["group1"].len(), 3);
/// assert_eq!(groups["group2"].len(), 2);
/// ```
fn scan_groups(input_dir: &PathBuf) -> HashMap<String, Vec<PathBuf>> {
    let mut groups: HashMap<String, Vec<PathBuf>> = HashMap::new();

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
        // INFO: basename -> {prefix}.{name}.{*}.bam -> {name}
        let bam = entry.path();
        let basename = bam
            .to_string_lossy()
            .split(".")
            .nth(1)
            .expect(&format!(
                "ERROR: Failed to get basename from {}",
                entry.path().display()
            ))
            .to_string();

        groups.entry(basename).or_insert(Vec::new()).push(bam);
    }

    groups
}
