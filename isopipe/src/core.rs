use crate::{config::*, executor::manager::ParallelExecutor};
use std::path::PathBuf;

pub mod ccs;
pub mod isoseq;
pub mod isotools;
pub mod lima;
pub mod minimap;
pub mod pbindex;
pub mod polya;
pub mod samtools;

pub fn run(
    config: Config,
    global_output_dir: PathBuf,
    mut executor: ParallelExecutor,
) -> Result<(), Box<dyn std::error::Error>> {
    log::info!("SUCCESS: All dependecies are loaded, starting pipeline...");
    // log::info!("INFO: Running with the following config: {:#?}", config);

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
            let input_dir = &global_output_dir.join(input_dir);

            samtools::merge(input_dir, executor, config, prefix);
            lima::lima(step, config, input_dir, &step_output_dir)
        }
        PipelineStep::Refine => {
            log::info!("INFO [STEP 3]: Pre-processing for isoseq::refine started...");
            isoseq::refine(
                step,
                config,
                &global_output_dir.join(input_dir),
                &step_output_dir,
            )
        }
        PipelineStep::Cluster => {
            log::info!("INFO [STEP 4]: Pre-processing for isoseq::cluster started...");
            isoseq::cluster(
                step,
                config,
                &global_output_dir.join(input_dir),
                &step_output_dir,
            )
        }
        PipelineStep::Minimap => {
            log::info!("INFO [STEP 5]: Pre-processing for minimap started...");
            minimap::minimap2(
                step,
                config,
                &global_output_dir.join(input_dir),
                &step_output_dir,
            )
        }
        PipelineStep::Polya => {
            log::info!("INFO [STEP 6]: Pre-processing for polya started...");
            polya::polya(
                step,
                config,
                &global_output_dir.join(input_dir),
                &step_output_dir,
            )
        }
        PipelineStep::Fusion => {
            log::info!("INFO [STEP 7]: Pre-processing for iso-fusion started...");
            isotools::iso_fusion(
                step,
                config,
                &global_output_dir.join(input_dir),
                &step_output_dir,
            )
        }
        PipelineStep::Orf => {
            todo!()
        }
    };

    executor
        .add_jobs(jobs)
        .execute(config, step, global_output_dir.clone());
}
