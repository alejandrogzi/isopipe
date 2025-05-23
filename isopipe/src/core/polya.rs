use crate::{cat, config::*, consts::*, executor::job::Job, isotools, rm};
use std::{fs::create_dir_all, path::PathBuf};

/// Run polya mod [3 steps]
///
/// # Arguments
///
/// * `step` - The pipeline step to run
/// * `config` - The configuration to use
/// * `input_dir` - The input directory
/// * `output_dir` - The output directory
///
/// # Returns
///
/// A vector of jobs to run
///
/// # Examples
///
/// ```rust, no_run
/// use isopipe::core::polya;
/// use isopipe::config::Config;
/// use std::path::PathBuf;
///
/// let config = Config::default();
/// let input_dir = PathBuf::from("/path/to/input");
/// let output_dir = PathBuf::from("/path/to/output");
///
/// let jobs = polya(&step, &config, &input_dir, &output_dir);
/// ```
pub fn polya(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
) -> Vec<Job> {
    let binary = isotools!(ISO_POLYA);
    let mut jobs = Vec::new();

    let parts = &step_output_dir.join(POLYA_PARTS);
    let _ = create_dir_all(parts);

    rm!(input_dir.join(CHUNKS));

    let args = config.get_step_args(
        step,
        vec![INPUT_DIR, OUTPUT_DIR, MEMORY, TIME, TOGA, ASSEMBLY],
    );

    // WARN: input_dir needs to be suffixed by /bam
    for entry in std::fs::read_dir(input_dir.join(BAM))
        .expect("Failed to read assets directory")
        .flatten()
        .filter(|entry| {
            entry
                .path()
                .file_name()
                .and_then(|ext| ext.to_str())
                .map(|name| name.ends_with(BAM))
                .unwrap_or(false)
        })
    {
        let bam = entry.path();

        // INFO: for cannonical runs this yields -> chunk*{run}.singletons
        // INFO: and outputs chunk*{run}.singletons.good.bed
        let bind = bam.with_extension("");
        let prefix = bind
            .file_stem()
            .unwrap_or_else(|| panic!("ERROR: could not build prefix for {:?}", bam));

        if std::fs::metadata(&bam).unwrap().len() == 0 {
            log::warn!("WARNING: {} its empty!", bam.display());
            continue;
        }

        let cmd = format!(
            "{} {} --bam {} {args} --prefix {} --outdir {} && rm {} {}.bai",
            binary.display(),
            SEGMENT,
            bam.display(),
            prefix.to_string_lossy(),
            &parts.display(),
            bam.display(),
            bam.display()
        );

        jobs.push(Job::from(cmd));
    }

    log::info!("INFO [STEP 6]: Pre-processing completed -> Running...");

    return jobs;
}

/// Merges polya .bed results into a single .bed per category.
/// Where categories are: singletons, accepts and rejections
///
/// # Arguments
///
/// * `input_dir` - The input directory
///
/// # Example
///
/// ```rust, no_run
/// let input_dir = PathBuf::from("/path/to/input");
/// merge(input_dir);
/// ```
pub fn merge(input_dir: &PathBuf) {
    log::info!(
        "INFO [MERGE]: Merging polyA segmentation results in {}...",
        &input_dir.join(POLYA_PARTS).display()
    );

    let files = std::fs::read_dir(input_dir.join(POLYA_PARTS))
        .expect("Failed to read assets directory")
        .flatten()
        .map(|entry| entry.path())
        .filter(|entry| {
            entry
                .file_name()
                .and_then(|ext| ext.to_str())
                .map(|name| name.ends_with(BED))
                .unwrap_or(false)
        })
        .collect::<Vec<_>>();

    if files.is_empty() {
        log::error!("ERROR: could not find any .bed in {}!", input_dir.display());
        std::process::exit(1);
    }

    // INFO: partition paths into their respective categories
    let (singletons, accepts, rejections): (Vec<_>, Vec<_>, Vec<_>) = files.into_iter().fold(
        (Vec::new(), Vec::new(), Vec::new()),
        |(mut s, mut g, mut b), path| {
            match &path {
                p if p.ends_with(BED_SGN_ACCEPT) => s.push(path),
                p if p.ends_with(BED_ACCEPT) => g.push(path),
                p if p.ends_with(BED_REJECT) => b.push(path),
                _ => {}
            }
            (s, g, b)
        },
    );

    for (file, group) in ALN_POLYA_FILES
        .iter()
        .zip([singletons, accepts, rejections].iter())
    {
        maybe_cat(group, file);
    }

    rm!(input_dir.join(POLYA_PARTS));
}

/// Wrapper fn around cat!() to merge paths into a target
///
/// # Arguments
///
/// * `paths` - Collection of paths to write
/// * `target` - Path to output
///
/// # Example
///
/// ```rust, no_run
/// let target = PathBuf::from("/path/to/output");
/// let paths = vec!["file1","file2"];
/// maybe_cat(paths, target);
/// ```
fn maybe_cat(paths: &Vec<PathBuf>, target: impl AsRef<std::path::Path> + std::fmt::Debug) {
    if !paths.is_empty() {
        log::info!(
            "INFO [MERGE]: Merging {} files into {:?}...",
            paths.len(),
            target
        );
        cat!(paths, target);
    }
}
