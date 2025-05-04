use iso_fusion::lib_iso_fusion;

use crate::{config::*, consts::*, executor::job::Job};
use std::path::PathBuf;
use std::sync::Arc;

/// Run isotools iso-fusion
///
/// # Arguments
/// * `step` - The pipeline step
/// * `config` - The configuration
/// * `input_dir` - The input directory
/// * `step_output_dir` - The output directory
///
/// # Returns
/// A vector of jobs
///
/// # Example
/// ```
/// let jobs = iso_fusion(&step, &config, &input_dir, &step_output_dir);
/// ```
pub fn iso_fusion(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
) -> Vec<Job> {
    let parts = config
        .get_step_args(
            step,
            vec![INPUT_DIR, PREFIX, OUTPUT_DIR, MEMORY, TIME, PREFIX],
        )
        .split(" ")
        .map(String::from)
        .collect::<Vec<String>>();

    for category in CLUSTERING_CATEGORIES {
        if *category == "lq" {
            continue;
        }

        let mut args = Vec::new();

        let query = format!(
            "{}/{}.{}.{}",
            input_dir.display(),
            CU_ALN,
            category,
            CORR_MINIMAP_GOOD_BED
        );
        args.extend(vec![String::from("--query"), query]);

        let prefix = step_output_dir.join(format!("{}", category));
        args.extend(vec![String::from("--prefix"), prefix.display().to_string()]);

        std::fs::create_dir_all(&prefix).expect(&format!(
            "ERROR: Failed to create directory {}",
            prefix.display()
        ));

        args.extend(parts.clone());
        let _ = lib_iso_fusion(Arc::new(args));
    }

    let jobs = aggregate_fusions(step_output_dir);
    log::info!("INFO [STEP 7]: Pre-processing completed -> Running...");

    return jobs;
}

/// Aggregate fusions from all categories into a single file.
///
/// # Arguments
/// * `step_output_dir` - The output directory
///
/// # Returns
/// A vector of jobs
///
/// # Example
/// ```
/// let jobs = aggregate_fusions(&step_output_dir);
/// ```
fn aggregate_fusions(step_output_dir: &PathBuf) -> Vec<Job> {
    FUSION_TYPES
        .iter()
        .map(|ty| {
            let pattern = if *ty == "fusions" {
                format!("{0}/*/{1}.bed", step_output_dir.display(), ty)
            } else {
                format!("{0}/*/*.{1}.bed", step_output_dir.display(), ty)
            };
            let output = format!("{0}/fusions.{1}.bed", step_output_dir.display(), ty);
            Job::from(format!("cat {} > {}", pattern, output))
        })
        .collect()
}
