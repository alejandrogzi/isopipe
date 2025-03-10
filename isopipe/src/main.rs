use clap::{self, Parser};
use log::{error, info, Level};
use simple_logger::init_with_level;

use isopipe::{
    cli::{Args, SubArgs},
    core::run,
};

fn main() {
    let start = std::time::Instant::now();
    init_with_level(Level::Info).unwrap();

    let args: Args = Args::parse();

    match args.command {
        SubArgs::Run { args } => {
            // args.check().unwrap_or_else(|e| {
            //     error!("{}", e);
            //     std::process::exit(1);
            // });

            let config = isopipe::config::Config::read(args.config)
                .expect("ERROR: Could not read config file");
            config.load().expect("ERROR: Could not load config file");

            run(config).unwrap_or_else(|e| {
                error!("{}", e);
                std::process::exit(1);
            });
        }
        SubArgs::Step { args } => {
            // args.check().unwrap_or_else(|e| {
            //     error!("{}", e);
            //     std::process::exit(1);
            // });

            let mut config = isopipe::config::Config::read(args.config.clone())
                .expect("ERROR: Could not read config file");

            config
                .aware(args)
                .load()
                .expect("ERROR: Could not load config file");

            run(config).unwrap_or_else(|e| {
                error!("{}", e);
                std::process::exit(1);
            });
        }
        SubArgs::Write { args } => {
            todo!()
            // args.check().unwrap_or_else(|e| {
            //     error!("{}", e);
            //     std::process::exit(1);
            // });

            // write(args).unwrap_or_else(|e| {
            //     error!("{}", e);
            //     std::process::exit(1);
            // });
        }
    }

    let elapsed = start.elapsed();
    info!("Elapsed time: {:.3?}", elapsed);
}
