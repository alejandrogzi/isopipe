use crate::config::*;
use crate::executor::job::Job;

use std::path::PathBuf;

const BAM: &str = "bam";
const SAM: &str = "sam";

pub fn run(config: Config, global_output_dir: PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    log::info!("SUCCESS: All dependecies are loaded, starting pipeline...");
    log::info!("INFO: Running with the following config: {:#?}", config);

    // config.steps().iter().for_each(|step| {
    //     run_step(step, &config, output);
    // });

    Ok(())
}

pub fn run_step(step: &PipelineStep, config: &Config, global_output_dir: PathBuf) {
    let mut jobs = Vec::new();

    match step {
        PipelineStep::Ccs => {
            // ccs
            // ${raWData_Link}/${dataPrefix}.subreads.bam
            // ${P_out_CCS}/${dataPrefix}.ccs.$i.bam
            // --min-rq 0.95
            // --chunk $i/$nChunks
            // -j $defVars{'nThreadsIsoSeq3'}
            // --report-file ${runID}_ccs_report.txt

            let args = config
                .params()
                .get(&PipelineStep::Ccs)
                .expect("ERROR: ccs not found in config.toml!")
                .flat(Some(vec!["input_dir", "prefix", "output_dir", "chunk"]));

            let prefix = config.get_data_prefix();
            let input_dir = config
                .get_param(PipelineStep::Ccs, "input_dir")
                .to_path_buf();
            let chunks = config.get_param(PipelineStep::Ccs, "chunk").to_string();

            let step_output_dir = global_output_dir.join(
                config
                    .get_param(PipelineStep::Ccs, "output_dir")
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
                        .expect("ERROR: no extension found!")
                        == BAM
                })
                .enumerate()
            {
                let bam = entry.path();
                let out_bam = step_output_dir.join(format!("{}.ccs.{}.bam", prefix, idx));
                let chunks = format!("{}/{}", idx, chunks);

                let job = Job::new()
                    .task(PipelineStep::Ccs)
                    .arg(bam.to_str().expect("ERROR: failed to convert path to str"))
                    .arg(
                        out_bam
                            .to_str()
                            .expect("ERROR: failed to convert path to str"),
                    )
                    .arg(&chunks)
                    .arg(&args);

                jobs.push(job)
            }

            // executor.add_jobs(jobs).execute();
        }
        PipelineStep::Lima => {
            // lima
            // --isoseq
            // --peek-guess
            // ${P_out_CCS}/${dataPrefix}.ccs.${i}.bam
            // ${f_primers}
            // ${P_out_lima}/${dataPrefix}.fl.${i}.bam

            let args = config
                .params()
                .get(&PipelineStep::Lima)
                .expect("ERROR: ccs not found in config.toml!")
                .flat(Some(vec!["input_dir", "prefix", "output_dir"]));

            let prefix = config.get_data_prefix();
            let input_dir = config
                .get_param(PipelineStep::Lima, "input_dir")
                .to_path_buf();
            let primers = config.get_param(PipelineStep::Lima, "primers").to_string();

            let step_output_dir = global_output_dir.join(
                config
                    .get_param(PipelineStep::Lima, "output_dir")
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
                        .expect("ERROR: no extension found!")
                        == BAM
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

            // executor.add_jobs(jobs).execute();
        }
        PipelineStep::Refine => {
            // isoseq3
            // refine
            // ${P_out_lima}/${dataPrefix}.fl.${i}.${primer_p5}${primer_p3}.bam
            // ${f_primers}
            // ${P_out_isoR}/${dataPrefix}.flnc.${i}.bam
            // --num-threads $defVars{'nThreadsIsoSeq3'}
            let output = std::process::Command::new("isoseq3")
                .arg("refine")
                .arg("--help")
                .output()
                .expect("ERROR: Failed to execute process");

            if output.status.success() {
                log::info!("refine: {}", String::from_utf8_lossy(&output.stdout));
            } else {
                log::error!(
                    "ERROR: failed to execute refine\n{}",
                    String::from_utf8_lossy(&output.stderr)
                );
                std::process::exit(1);
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
            let output = std::process::Command::new("isoseq3")
                .arg("cluster")
                .arg("--help")
                .output()
                .expect("ERROR: Failed to execute process");

            if output.status.success() {
                log::info!("cluster: {}", String::from_utf8_lossy(&output.stdout));
            } else {
                log::error!(
                    "ERROR: failed to execute cluster\n{}",
                    String::from_utf8_lossy(&output.stderr)
                );
                std::process::exit(1);
            }
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
                .flat(Some(vec!["input_dir", "output_dir"]));

            let output_dir = format!(
                "-outdir {}",
                global_output_dir
                    .join(
                        config
                            .get_param(PipelineStep::FilterQuality, "output_dir")
                            .to_path_buf(),
                    )
                    .display()
            );
            let input_dir = config
                .get_param(PipelineStep::FilterQuality, "input_dir")
                .to_path_buf();

            let mut file_count = 0;

            for entry in std::fs::read_dir(input_dir)
                .expect("Failed to read assets directory")
                .flatten()
                .filter(|entry| {
                    entry
                        .path()
                        .extension()
                        .expect("ERROR: no extension found!")
                        == SAM
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

            // dbg!(jobs);
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
