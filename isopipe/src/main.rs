use clap::{self, Parser};
use log::{error, info};
use simple_logger::init_with_level;

use isopipe::{
    cli::{Args, SubArgs},
    core::{__tags, run},
    executor::manager::ParallelManager,
};

#[allow(unused_variables)]
fn main() {
    let start = std::time::Instant::now();
    let args: Args = Args::parse();

    init_with_level(args.level).unwrap();

    let executor = match args.manager {
        ParallelManager::Nextflow | ParallelManager::Para => {
            info!("INFO: Initializing parallel executor...");
            args.manager.init()
        }
        _ => {
            error!("ERROR: Unknown executor or has not been implemented!");
            std::process::exit(1);
        }
    };

    match args.command {
        SubArgs::Run { args } => {
            let config = isopipe::config::Config::read(args.config)
                .unwrap_or_else(|e| panic!("ERROR: Could not read config file -> {e}"))
                .load()
                .check();

            let global_output_dir = config.create_global_output_dir();

            run(config, global_output_dir, executor).unwrap_or_else(|e| {
                error!("{}", e);
                std::process::exit(1);
            });
        }
        SubArgs::Step { args } => {
            let config = isopipe::config::Config::read(args.config.clone())
                .unwrap_or_else(|e| panic!("ERROR: Could not read config file -> {e}"))
                .aware(args.clone())
                .load()
                .check();

            let global_output_dir = config.create_global_output_dir();

            run(config, global_output_dir, executor).unwrap_or_else(|e| {
                error!("{}", e);
                std::process::exit(1);
            });
        }
        SubArgs::Write { args } => todo!(),
        SubArgs::Tag { args } => __tags(args),
    }

    let elapsed = start.elapsed();
    info!("Elapsed time: {:.3?}", elapsed);
}
