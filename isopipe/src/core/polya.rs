use crate::{
    config::*,
    consts::*,
    executor::{job::Job, manager::__get_assets_dir},
};
use std::path::PathBuf;

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
/// ```
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
    output_dir: &PathBuf,
) -> Vec<Job> {
    let mut jobs = Vec::new();

    let args = config.get_step_args(
        step,
        vec![INPUT_DIR, OUTPUT_DIR, MEMORY, TIME, TOGA, ASSEMBLY],
    );
    let fields = config.get_step_custom_fields(step, vec![TOGA, ASSEMBLY]);
    let assets = __get_assets_dir();

    let filter = assets.join(FILTER_MINIMAP);
    let correct = assets.join(CORRECT_MINIMAP);

    for category in CLUSTERING_CATEGORIES {
        if *category == "lq" {
            continue;
        }

        // INFO: format -> all.clustered.aligned.{hq,lq,singletons}.sam
        let filename = PathBuf::from(format!("{}.{}.{}", CU_ALN, category, SAM));
        let alignment = input_dir.join(&filename);

        // INFO: will output all.clustered.aligned.{hq,lq,singletons}.{good,bad}.sam
        // INFO: script.perl {].sam --perID 96 --clip3 50 --polyAReadSuffix 30 --outdir {}/first_pass
        let first_pass = format!(
            "{} {} {} --outdir {}",
            filter.display(),
            alignment.display(),
            args,
            output_dir.join(POLYA_FIRST_PASS).display()
        );

        // INFO: script.py {toga} {}.good.sam {assembly} {].corrected.sam
        let corrected_sam = output_dir.join(filename.with_extension(CORR_MINIMAP_SAM));
        let correct_step = format!(
            "python3 {} {} {} {} {}",
            correct.display(),
            fields
                .get(0)
                .expect(&format!("ERROR: Could not find TOGA -> {:?}", fields)),
            output_dir
                .join(POLYA_FIRST_PASS)
                .join(filename.with_extension(POLYA_GOOD_SAM))
                .display(),
            fields
                .get(1)
                .expect(&format!("ERROR: Could not find assembly -> {:?}", fields)),
            corrected_sam.display()
        );

        // INFO: script.perl {}.corrected.sam --polyAReadSuffix 30 --outdir {}
        let second_pass = format!(
            "{} {} -polyAReadSuffix 30 --outdir {}",
            filter.display(),
            corrected_sam.display(),
            output_dir.display()
        );

        let convert = format!(
            "{} {} -i {} -bed12 > {}",
            BEDTOOLS,
            BAMTOBED,
            output_dir
                .join(filename.with_extension(CORR_MINIMAP_GOOD_SAM))
                .display(),
            output_dir
                .join(filename.with_extension(CORR_MINIMAP_GOOD_BED))
                .display()
        );

        jobs.push(Job::from(format!(
            "{} && {} && {} && {}",
            first_pass, correct_step, second_pass, convert
        )));
    }

    log::info!("INFO [STEP 6]: Pre-processing completed -> Running...");

    return jobs;
}
