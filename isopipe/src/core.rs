use crate::config::*;

pub fn run(config: Config) -> Result<(), Box<dyn std::error::Error>> {
    log::info!("SUCCESS: All dependecies are loaded, starting pipeline...");
    log::info!("INFO: Running with the following config: {:#?}", config);

    Ok(())
}

pub fn run_step(step: &PipelineStep) {
    match step {
        PipelineStep::Ccs => {
            let output = std::process::Command::new("ccs")
                .arg("--help")
                .output()
                .expect("ERROR: Failed to execute process");

            if output.status.success() {
                log::info!("ccs: {}", String::from_utf8_lossy(&output.stdout));
            } else {
                log::error!(
                    "ERROR: failed to execute ccs\n{}",
                    String::from_utf8_lossy(&output.stderr)
                );
                std::process::exit(1);
            }
        }
        PipelineStep::Lima => {
            let output = std::process::Command::new("lima")
                .arg("--help")
                .output()
                .expect("ERROR: Failed to execute process");

            if output.status.success() {
                log::info!("lima: {}", String::from_utf8_lossy(&output.stdout));
            } else {
                log::error!(
                    "ERROR: failed to execute lima\n{}",
                    String::from_utf8_lossy(&output.stderr)
                );
                std::process::exit(1);
            }
        }
        PipelineStep::Refine => {
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
            let output = std::process::Command::new("minimap3")
                .arg("--help")
                .output()
                .expect("ERROR: Failed to execute process");

            if output.status.success() {
                log::info!("minimap: {}", String::from_utf8_lossy(&output.stdout));
            } else {
                log::error!(
                    "ERROR: failed to execute minimap\n{}",
                    String::from_utf8_lossy(&output.stderr)
                );
                std::process::exit(1);
            }
        }
        PipelineStep::FilterQuality => {
            let output = std::process::Command::new("isotools")
                .arg("--bin")
                .arg("iso-polya")
                .arg("filter")
                .arg("--help")
                .output()
                .expect("ERROR: Failed to execute process");

            if output.status.success() {
                log::info!(
                    "filter_quality: {}",
                    String::from_utf8_lossy(&output.stdout)
                );
            } else {
                log::error!(
                    "ERROR: failed to execute filter_quality\n{}",
                    String::from_utf8_lossy(&output.stderr)
                );
                std::process::exit(1);
            }
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
