use crate::config::*;
use crate::executor::job::Job;
use crate::executor::manager::ParallelExecutor;

use std::path::PathBuf;

const BAM: &str = "bam";
const SAM: &str = "sam";

const OUTPUT_DIR: &str = "output_dir";
const PRIMERS: &str = "primers";
const INPUT_DIR: &str = "input_dir";
const MEMORY: &str = "memory";
const TIME: &str = "time";
const PREFIX: &str = "prefix";

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
    let mut jobs = Vec::new();

    match step {
        PipelineStep::Ccs => {
            log::info!("INFO [STEP0]: Pre-processing for ccs started...");

            let args = config
                .params()
                .get(&PipelineStep::Ccs)
                .expect("ERROR: ccs not found in config.toml!")
                .flat(Some(vec![
                    INPUT_DIR,
                    PREFIX,
                    OUTPUT_DIR,
                    "chunk",
                    MEMORY,
                    TIME,
                    "report-file",
                ]));

            let prefix = config.get_data_prefix();
            let input_dir = config
                .get_param(PipelineStep::Ccs, INPUT_DIR)
                .expect("ERROR: input_dir not found for ccs in config.toml!")
                .to_path_buf();
            let chunks = config
                .get_param(PipelineStep::Ccs, "chunk")
                .expect("ERROR: chunk not found for ccs in config.toml!")
                .to_string();
            let report = config
                .get_param(PipelineStep::Ccs, "report-file")
                .expect("ERROR: report-file not found for ccs in config.toml!")
                .to_string();

            let step_output_dir = global_output_dir.join(
                config
                    .get_param(PipelineStep::Ccs, OUTPUT_DIR)
                    .expect("ERROR: output_dir not found for ccs in config.toml!")
                    .to_path_buf(),
            );
            std::fs::create_dir_all(&step_output_dir).expect("ERROR: failed to create output dir!");

            // WARN: ignoring prefix + .subreads ending -> forcing to isolate samples
            for (idx, entry) in std::fs::read_dir(input_dir)
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
                .enumerate()
            {
                let idx = idx + 1; // WARN: 1-indexed
                let bam = entry.path();
                let out_bam = step_output_dir.join(format!("{}.ccs.{}.bam", prefix, idx));
                let chunks = format!("--chunk {}/{}", idx, chunks);
                let report = format!("--report-file {}/{}", step_output_dir.display(), report);

                // WARN: need to check if bam has a .pbi file -> if not, run pbindex
                let mut pbi = bam.clone();
                pbi.set_extension("bam.pbi");
                if !pbi.exists() {
                    log::warn!(
                        "WARN: pbi file not found for {}, generating index...",
                        bam.display()
                    );

                    generate_pb_index(&bam);
                }

                let job = Job::new()
                    .task(PipelineStep::Ccs)
                    .arg(bam.to_str().expect("ERROR: failed to convert path to str"))
                    .arg(
                        out_bam
                            .to_str()
                            .expect("ERROR: failed to convert path to str"),
                    )
                    .arg(&chunks)
                    .arg(&args)
                    .arg(&report);

                jobs.push(job)
            }

            log::info!("INFO [STEP0]: Pre-processing completed -> Running...");
            executor
                .add_jobs(jobs)
                .execute(config, step, global_output_dir.clone());
        }
        PipelineStep::Lima => {
            // lima
            // --isoseq
            // --peek-guess
            // ${P_out_CCS}/${dataPrefix}.ccs.${i}.bam
            // ${f_primers}
            // ${P_out_lima}/${dataPrefix}.fl.${i}.bam
            log::info!("INFO [STEP1]: Pre-processing for lima started...");

            let args = config
                .params()
                .get(&PipelineStep::Lima)
                .expect("ERROR: ccs not found in config.toml!")
                .flat(Some(vec![
                    INPUT_DIR, PREFIX, OUTPUT_DIR, MEMORY, TIME, "primers",
                ]));

            let prefix = config.get_data_prefix();
            let input_dir = global_output_dir.join(
                config
                    .get_param(PipelineStep::Lima, INPUT_DIR)
                    .expect("ERROR: input_dir not found for lima in config.toml!")
                    .to_path_buf(),
            );
            let primers = config
                .get_param(PipelineStep::Lima, PRIMERS)
                .expect("ERROR: primers file path not found for lima in config.toml!")
                .to_string();

            let step_output_dir = global_output_dir.join(
                config
                    .get_param(PipelineStep::Lima, OUTPUT_DIR)
                    .expect("ERROR: output_dir not found for lima in config.toml!")
                    .to_path_buf(),
            );
            std::fs::create_dir_all(&step_output_dir).expect("ERROR: failed to create output dir!");

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
                let bam = entry.path();
                let identifier = bam
                    .file_stem()
                    .expect("ERROR: failed to get file stem")
                    .to_str()
                    .expect("ERROR: failed to convert path to str")
                    .split('.')
                    .last()
                    .expect("ERROR: failed to get last element from .bam name!");
                let out_bam = step_output_dir.join(format!("{}.fl.{}.bam", prefix, identifier));

                let job = Job::new()
                    .task(PipelineStep::Lima)
                    .arg(&args)
                    .arg(bam.to_str().expect("ERROR: failed to convert path to str"))
                    .arg(&primers)
                    .arg(
                        out_bam
                            .to_str()
                            .expect("ERROR: failed to convert path to str"),
                    );

                jobs.push(job)
            }

            log::info!("INFO [STEP1]: Pre-processing completed -> Running...");
            executor
                .add_jobs(jobs)
                .execute(config, step, global_output_dir.clone());
        }
        PipelineStep::Refine => {
            // isoseq3
            // refine
            // ${P_out_lima}/${dataPrefix}.fl.${i}.${primer_p5}${primer_p3}.bam
            // ${f_primers}
            // ${P_out_isoR}/${dataPrefix}.flnc.${i}.bam
            // --num-threads $defVars{'nThreadsIsoSeq3'}

            let args = config
                .params()
                .get(&PipelineStep::Refine)
                .expect("ERROR: ccs not found in config.toml!")
                .flat(Some(vec![INPUT_DIR, PREFIX, OUTPUT_DIR, MEMORY, TIME]));

            let prefix = config.get_data_prefix();
            let input_dir = config
                .get_param(PipelineStep::Refine, INPUT_DIR)
                .expect("ERROR: input_dir not found for refine in config.toml!")
                .to_path_buf();
            let primers = config
                .get_param(PipelineStep::Refine, PRIMERS)
                .expect("ERROR: primers file path not found for refine in config.toml!")
                .to_string();

            let step_output_dir = global_output_dir.join(
                config
                    .get_param(PipelineStep::Refine, OUTPUT_DIR)
                    .expect("ERROR: output_dir not found for refine in config.toml!")
                    .to_path_buf(),
            );
            std::fs::create_dir_all(&step_output_dir).expect("ERROR: failed to create output dir!");

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
                let bam = entry.path();

                let identifier = bam
                    .file_stem()
                    .expect("ERROR: failed to get file stem")
                    .to_str()
                    .expect("ERROR: failed to convert path to str")
                    .split('.')
                    .collect::<Vec<&str>>();
                let identifier = identifier.get(2).expect("ERROR: failed to get identifier");
                let out_bam = step_output_dir.join(format!("{}.flnc.{}.bam", prefix, identifier));

                let job = Job::new()
                    .task(PipelineStep::Refine)
                    .arg(&args)
                    .arg(bam.to_str().expect("ERROR: failed to convert path to str"))
                    .arg(&primers)
                    .arg(
                        out_bam
                            .to_str()
                            .expect("ERROR: failed to convert path to str"),
                    );

                jobs.push(job)
            }
        }
        PipelineStep::Cluster => {
            // isoseq3
            // cluster
            // ${P_out_isoC}/ALL.flnc.fofn ${P_out_isoC}/ALL.CuP.bam
            // --singletons
            // -verbose
            // --split-bam $splitBAM
            // --num-threads $defVars{'nThreadsIsoSeq3'}
            // --log-file ${P_out_isoC}/ALL.CuP.log

            let args = config
                .params()
                .get(&PipelineStep::Cluster)
                .expect("ERROR: ccs not found in config.toml!")
                .flat(Some(vec![INPUT_DIR, PREFIX, OUTPUT_DIR, MEMORY, TIME]));

            let prefix = config.get_data_prefix();
            let input_dir = config
                .get_param(PipelineStep::Cluster, INPUT_DIR)
                .expect("ERROR: input_dir not found for cluster in config.toml!")
                .to_path_buf();

            let step_output_dir = global_output_dir.join(
                config
                    .get_param(PipelineStep::Cluster, OUTPUT_DIR)
                    .expect("ERROR: output_dir not found for cluster in config.toml!")
                    .to_path_buf(),
            );
            std::fs::create_dir_all(&step_output_dir).expect("ERROR: failed to create output dir!");

            let mut file_count = 0;

            // for entry in std::fs::read_dir(input_dir)
            //     .expect("Failed to read assets directory")
            //     .flatten()
            //     .filter(|entry| {
            //         entry
            //             .path()
            //             .extension()
            //             .expect("ERROR: no extension found!")
            //             == BAM
            //     })
            // {
            //     file_count += 1;

            //     let job = Job::new()
            //         .task(PipelineStep::Cluster)
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
        }
        PipelineStep::FilterQuality => {
            // isotools iso-polya filter
            //
            // filterMinimapQuality.perl
            // $input.sam
            // -perID 96
            // -clip3 50
            // -polyAReadSuffix 30

            let args = config
                .params()
                .get(&PipelineStep::FilterQuality)
                .expect("ERROR: filter-quality not found in config.toml!")
                .flat(Some(vec![INPUT_DIR, OUTPUT_DIR, MEMORY, TIME]));

            let output_dir = format!(
                "-outdir {}",
                global_output_dir
                    .join(
                        config
                            .get_param(PipelineStep::FilterQuality, OUTPUT_DIR)
                            .expect(
                                "ERROR: output_dir not found for filter-quality in config.toml!"
                            )
                            .to_path_buf(),
                    )
                    .display()
            );
            let input_dir = config
                .get_param(PipelineStep::FilterQuality, INPUT_DIR)
                .expect("ERROR: input_dir not found for filter-quality in config.toml!")
                .to_path_buf();

            let mut file_count = 0;

            for entry in std::fs::read_dir(input_dir)
                .expect("Failed to read assets directory")
                .flatten()
                .filter(|entry| {
                    entry
                        .path()
                        .extension()
                        .and_then(|ext| ext.to_str())
                        .map(|ext| ext.eq_ignore_ascii_case(SAM))
                        .unwrap_or(false)
                })
            {
                file_count += 1;

                let sam = format!("--sam {}", entry.path().display());

                let job = Job::new()
                    .task(PipelineStep::FilterQuality)
                    .arg(&sam)
                    .arg(&args)
                    .arg(&output_dir);

                jobs.push(job);

                if file_count > 1 {
                    log::error!(
                        "ERROR: more than one .sam file found in input_dir. This is a bug!"
                    );
                    std::process::exit(1);
                }
            }

            executor
                .add_jobs(jobs)
                .execute(config, step, global_output_dir.clone());
        }
        PipelineStep::LoadGenome => {
            let output = std::process::Command::new("load_genome")
                .arg("--help")
                .output()
                .expect("ERROR: Failed to execute process");

            if output.status.success() {
                log::info!("load_genome: {}", String::from_utf8_lossy(&output.stdout));
            } else {
                log::error!(
                    "ERROR: failed to execute load_genome\n{}",
                    String::from_utf8_lossy(&output.stderr)
                );
                std::process::exit(1);
            }
        }
    }
}

fn generate_pb_index(bam: &PathBuf) {
    let output = std::process::Command::new("pbindex")
        .arg(bam)
        .output()
        .expect("ERROR: Failed to execute process");

    if output.status.success() {
        log::info!(
            "INFO [PBINDEX]: Successfully generated pbindex for {}",
            bam.display()
        );
    } else {
        log::error!(
            "ERROR: failed to execute pbindex\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        std::process::exit(1);
    }
}
