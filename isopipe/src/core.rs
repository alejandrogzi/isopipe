use crate::{config::*, executor::manager::ParallelExecutor};
use std::path::PathBuf;

pub mod ccs;
pub mod isoseq;
pub mod isotools;
pub mod lima;
pub mod minimap;
pub mod orf;
pub mod pbindex;
pub mod polya;
pub mod samtools;

pub fn run(
    config: Config,
    global_output_dir: PathBuf,
    mut executor: ParallelExecutor,
) -> Result<(), Box<dyn std::error::Error>> {
    log::info!("SUCCESS: All dependecies are loaded, starting pipeline...");

    config.steps().iter().for_each(|step| {
        run_step(step, &config, &global_output_dir, &mut executor);
    });

    Ok(())
}

pub fn run_step(
    step: &PipelineStep,
    config: &Config,
    global_output_dir: &PathBuf,
    executor: &mut ParallelExecutor,
) {
    let prefix = config.get_data_prefix();
    let (input_dir, step_output_dir) = config.get_step_dirs(step, global_output_dir);

    let jobs = match step {
        PipelineStep::Ccs => {
            log::info!("INFO [STEP 1]: Pre-processing for ccs started...");
            ccs::ccs(step, config, &input_dir, &step_output_dir, prefix, executor)
        }
        PipelineStep::Lima => {
            log::info!("INFO [STEP 2]: Pre-processing for lima started...");

            samtools::merge(&input_dir, executor, config, prefix);
            lima::lima(step, config, &input_dir, &step_output_dir)
        }
        PipelineStep::Refine => {
            log::info!("INFO [STEP 3]: Pre-processing for isoseq::refine started...");
            isoseq::refine(step, config, &input_dir, &step_output_dir)
        }
        PipelineStep::Cluster => {
            log::info!("INFO [STEP 4]: Pre-processing for isoseq::cluster started...");
            isoseq::cluster(step, config, &input_dir, &step_output_dir)
        }
        PipelineStep::Minimap => {
            log::info!("INFO [STEP 5]: Pre-processing for minimap started...");

            let input_dir =
                isotools::iso_split(step, config, &input_dir, &step_output_dir, executor);
            minimap::minimap2(step, config, &input_dir, &step_output_dir)
        }
        PipelineStep::Polya => {
            log::info!("INFO [STEP 6]: Pre-processing for polya started...");

            samtools::index(config, &input_dir, executor);
            polya::polya(step, config, &input_dir, &step_output_dir)
        }
        PipelineStep::Fusion => {
            log::info!("INFO [STEP 7]: Pre-processing for iso-fusion started...");

            polya::merge(&input_dir);
            isotools::iso_fusion(step, config, &input_dir, &step_output_dir)
        }
        PipelineStep::Orf => {
            log::info!("INFO [STEP 8]: Pre-processing for orf started...");
            orf::orf(step, config, &input_dir, &step_output_dir)
        }
        PipelineStep::Polish => {
            log::info!("INFO [STEP 9]: Pre-processing for isotools polish started...");
            isotools::polish(step, config, &input_dir, &step_output_dir, executor)
        }
    };

    executor
        .add_jobs(jobs)
        .execute(config, step, global_output_dir.clone());
}
