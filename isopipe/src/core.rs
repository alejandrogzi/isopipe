use crate::{config::*, executor::manager::ParallelExecutor};
use std::path::PathBuf;

pub mod ccs;
pub mod isoseq;
pub mod lima;

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
            log::info!("INFO [STEP 0]: Pre-processing for ccs started...");
            ccs::ccs(step, config, &input_dir, &step_output_dir, prefix)
        }
        PipelineStep::Lima => {
            log::info!("INFO [STEP 1]: Pre-processing for lima started...");
            lima::lima(step, config, &input_dir, &step_output_dir, prefix)
        }
        PipelineStep::Refine => {
            log::info!("INFO [STEP 2]: Pre-processing for isoseq::refine started...");
            isoseq::refine(step, config, &input_dir, &step_output_dir, prefix)
        }
        PipelineStep::Cluster => {
            log::info!("INFO [STEP 3]: Pre-processing for isoseq::cluster started...");
            isoseq::cluster(step, config, &input_dir, &step_output_dir, prefix)
        }
        PipelineStep::Minimap => {
            // let output = std::process::Command::new("minimap3")
            //     .arg("--help")
            //     .output()
            //     .expect("ERROR: Failed to execute process");

            // if output.status.success() {
            //     log::info!("minimap: {}", String::from_utf8_lossy(&output.stdout));
            // } else {
            //     log::error!(
            //         "ERROR: failed to execute minimap\n{}",
            //         String::from_utf8_lossy(&output.stderr)
            //     );
            //     std::process::exit(1);
            // }
            //

            todo!()
        }
        PipelineStep::FilterQuality => {
            todo!()
            // isotools iso-polya filter
            //
            // filterMinimapQuality.perl
            // $input.sam
            // -perID 96
            // -clip3 50
            // -polyAReadSuffix 30

            // let args = config
            //     .params()
            //     .get(&PipelineStep::FilterQuality)
            //     .expect("ERROR: filter-quality not found in config.toml!")
            //     .flat(Some(vec![INPUT_DIR, OUTPUT_DIR, MEMORY, TIME]));

            // let output_dir = format!(
            //     "-outdir {}",
            //     global_output_dir
            //         .join(
            //             config
            //                 .get_param(PipelineStep::FilterQuality, OUTPUT_DIR)
            //                 .expect(
            //                     "ERROR: output_dir not found for filter-quality in config.toml!"
            //                 )
            //                 .to_path_buf(),
            //         )
            //         .display()
            // );

            // let mut file_count = 0;

            // for entry in std::fs::read_dir(input_dir)
            //     .expect("Failed to read assets directory")
            //     .flatten()
            //     .filter(|entry| {
            //         entry
            //             .path()
            //             .extension()
            //             .and_then(|ext| ext.to_str())
            //             .map(|ext| ext.eq_ignore_ascii_case(SAM))
            //             .unwrap_or(false)
            //     })
            // {
            //     file_count += 1;

            //     let sam = format!("--sam {}", entry.path().display());

            //     let job = Job::new()
            //         .task(PipelineStep::FilterQuality)
            //         .arg(&sam)
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
        PipelineStep::LoadGenome => {
            todo!()
        }
    };

    executor
        .add_jobs(jobs)
        .execute(config, step, global_output_dir.clone());
}
