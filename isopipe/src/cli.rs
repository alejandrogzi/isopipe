use clap::{Parser, Subcommand};
use std::path::PathBuf;

use crate::config::*;
use crate::executor::manager::ParallelManager;

pub const MIN_STEP: &str = "0";
pub const MAX_STEP: &str = "8";

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    #[command(subcommand)]
    pub command: SubArgs,

    #[arg(
        short = 'm',
        long = "manager",
        help = "Parallel executor strategy",
        value_name = "MANAGER",
        required = false,
        default_value = "para"
    )]
    pub manager: ParallelManager,
}

impl Args {}

#[derive(Debug, Subcommand)]
pub enum SubArgs {
    #[command(name = "run")]
    Run {
        #[command(flatten)]
        args: RunArgs,
    },
    #[command(name = "run-step")]
    Step {
        #[command(flatten)]
        args: StepArgs,
    },
    #[command(name = "write")]
    Write {
        #[command(flatten)]
        args: WriteArgs,
    },
}

/// Run the pipeline from start to finish
///
/// # Example
///
/// ```bash,no_run
/// isopipe run -c config.toml
/// ```
///
/// # Arguments
///
/// * `config` - Path to the configuration file
///
/// # Note
///
/// * The pipeline will run from start to finish
/// * If not config.toml is provided, it will default to config.toml
#[derive(Debug, Parser)]
pub struct RunArgs {
    #[arg(
        short = 'c',
        long = "config",
        help = "Path to the configuration file",
        value_name = "STEP",
        required = true,
        default_value = "config.toml"
    )]
    pub config: PathBuf,
}

// impl ArgCheck for RunArgs {}

/// Run the pipeline from a specific step to another step
/// or only run a specific step
/// or skip a specific step or multiple steps
///
/// # Example
///
/// ```bash,no_run
/// isopipe run-step -c config.toml -f 0 -t 7
/// isopipe run-step -c config.toml -o 0
/// isopipe run-step -c config.toml -s 0,1,2
/// ```
///
/// # Arguments
///
/// * `config` - Path to the configuration file
/// * `from` - Start from a specific step
/// * `to` - Stop at a specific step
/// * `only` - Only run a specific step
/// * `skip` - Skip a specific step or multiple steps
/// * `dry_run` - Dry run the pipeline
/// * `verbose` - Increase verbosity
/// * `quiet` - Decrease verbosity
///
/// # Note
///
/// * `from` and `to` are mutually exclusive
/// * `only` and `skip` are mutually exclusive
/// * `skip` can take multiple values
///
/// # Warning
///
/// * Assuming 9 steps in the pipeline
///
#[derive(Debug, Parser, Clone)]
pub struct StepArgs {
    #[arg(
        short = 'c',
        long = "config",
        help = "Path to the configuration file",
        value_name = "STEP",
        required = true,
        default_value = "config.toml"
    )]
    pub config: PathBuf,

    #[arg(
        short = 'f',
        long = "from",
        help = "Start from a specific step. Can be a step number or step name.",
        value_name = "STEP",
        default_value = MIN_STEP,
        conflicts_with = "only"
    )]
    pub from: String,

    #[arg(
        short = 't',
        long = "to",
        help = "Stop at a specific step. Can be a step number or step name.",
        value_name = "STEP",
        conflicts_with = "only",
        requires = "from",
        default_value = MAX_STEP
    )]
    pub to: String,

    #[arg(
        short = 'o',
        long = "only",
        help = "Only run a specific step (or steps). Specify the step number or step name.",
        value_name = "STEP",
        value_delimiter = ',',
        conflicts_with = "to",
        conflicts_with = "skip",
        conflicts_with = "from",
        num_args = 1..,
    )]
    pub only: Option<Vec<String>>,

    #[arg(
        short = 's',
        long = "skip",
        help = "Skip a specific step or multiple steps. Specify the step number or step name.",
        value_delimiter = ',',
        value_name = "STEPs",
        conflicts_with = "only",
        num_args = 1..,
    )]
    pub skip: Option<Vec<String>>,

    #[arg(short = 'd', long = "dry-run", help = "Dry run the pipeline")]
    pub dry_run: bool,

    #[arg(short = 'v', long = "verbose", help = "Increase verbosity")]
    pub verbose: bool,

    #[arg(short = 'q', long = "quiet", help = "Decrease verbosity")]
    pub quiet: bool,
}

impl StepArgs {
    /// Build an absolute list of steps to run based on args
    ///
    /// # Returns
    ///
    /// * A list of steps to run
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// let args = StepArgs {
    ///    config: PathBuf::from("config.toml"),
    ///   from: "0".to_string(),
    ///  to: "6".to_string(),
    /// only: None,
    /// skip: None,
    /// dry_run: false,
    /// verbose: false,
    /// quiet: false,
    /// };
    ///
    /// let steps = args.abs_steps().unwrap();
    ///
    /// assert_eq!(steps.len(), 7);
    /// ```
    pub fn abs_steps(&self) -> Result<Vec<PipelineStep>, Box<dyn std::error::Error>> {
        let max_step = MAX_STEP
            .parse::<usize>()
            .expect("ERROR: Could not parse max step!");

        fn parse_step(step: &str) -> Result<usize, Box<dyn std::error::Error>> {
            step.parse::<usize>()
                .or_else(|_| PipelineStep::from_str(step).map(|s| s.to_int()))
                .map_err(|_| format!("ERROR: invalid step '{}'", step).into())
        }

        fn validate_step(
            step: usize,
            max: usize,
            flag: &str,
        ) -> Result<usize, Box<dyn std::error::Error>> {
            if step > max {
                return Err(
                    format!("ERROR: --{} must be less than or equal to {}", flag, max).into(),
                );
            }
            Ok(step)
        }

        let from = validate_step(parse_step(&self.from)?, max_step, "from")?;
        let to = validate_step(parse_step(&self.to)?, max_step, "to")?;

        if from > to {
            return Err("ERROR: --from must be less than --to".into());
        } else if from == to {
            return Err("ERROR: --from must not be equal to --to".into());
        }

        if let Some(only) = &self.only {
            let mut steps: Vec<usize> = only
                .iter()
                .map(|s| validate_step(parse_step(s)?, max_step, "only"))
                .collect::<Result<Vec<_>, _>>()?;
            steps.sort_unstable();

            log::info!("INFO: running step/s {:?} only...", steps);

            return Ok(steps
                .into_iter()
                .map(|s| {
                    PipelineStep::from_int(s).expect("ERROR: Tried to parse an invalid step index!")
                })
                .collect());
        }

        let mut steps = (from..=to).collect::<Vec<_>>();
        steps.sort_unstable();

        let skips = if let Some(skip) = &self.skip {
            skip.iter()
                .map(|s| validate_step(parse_step(s)?, max_step, "skip"))
                .collect::<Result<Vec<_>, _>>()?
        } else {
            Vec::new()
        };

        let result_steps = steps
            .into_iter()
            .filter(|s| !skips.contains(s))
            .map(|s| {
                PipelineStep::from_int(s).expect("ERROR: Tried to parse an invalid step index!")
            })
            .collect::<Vec<_>>();

        if result_steps.is_empty() {
            return Err("ERROR: No steps to run".into());
        } else if result_steps.len() == max_step {
            log::warn!("WARN: Running all steps... Next time use run instead of run-step!");
        }

        Ok(result_steps)
    }
}

// impl ArgCheck for StepArgs {}

/// Write commands to a .sh file
/// with parameters specified in --config
///
/// # Example
///
/// ```bash,no_run
/// isopipe write -c config.toml --cmd ccs,lima,refine
/// ```
///
/// # Arguments
///
/// * `config` - Path to the configuration file
/// * `cmd` - Command(s) to write to a .sh file with paramenters specified in --config
///
/// # Note
///
/// * `cmd` can take multiple values
/// * If not config.toml is provided, it will default to config.toml
///
#[derive(Debug, Parser, Clone)]
pub struct WriteArgs {
    #[arg(
        short = 'c',
        long = "config",
        help = "Path to the configuration file",
        value_name = "STEP",
        required = true,
        default_value = "config.toml"
    )]
    pub config: PathBuf,

    #[arg(
        long = "cmd",
        help = "Command(s) to write to a .sh file with paramenters specified in --config",
        value_delimiter = ',',
        value_name = "CMD",
        num_args = 1..,
    )]
    pub cmd: Vec<String>,
}
