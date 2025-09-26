// collapse run --bed 1.bed,2.bed,3.bed --index [binary or .gz?] --write [writes each collpased file always to /collapsed] --outdir
// collapse read --index <INDEX> --read <READ> [print to stdout] --write [writes utf8]
// collapse run --bed 1.bed --extend [appends what would be --index as an additional column + automatically writes + merges all beds] --outdir

use clap::{self, ArgAction, ArgGroup, Parser};

use std::path::PathBuf;

#[derive(Parser, Debug)]
#[clap(name="collapse", author = env!("CARGO_PKG_AUTHORS"), version = env!("CARGO_PKG_VERSION"), about = "shrink your .beds", long_about = None)]
pub struct Args {
    #[clap(subcommand)]
    pub command: Command,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        short = 'L',
        long = "level",
        help = "Logging level",
        value_name = "LEVEL",
        default_value_t = log::Level::Info,
    )]
    pub level: log::Level,
}

#[derive(Parser, Debug)]
pub enum Command {
    Run(RunArgs),
    Read(ReadArgs),
}

#[derive(Parser, Debug)]
#[command(group(ArgGroup::new("indexing").required(false).multiple(true).args(&["index", "write"])))]
#[command(group(ArgGroup::new("mode").required(true).multiple(true).args(&["index", "extend"])))]
#[command(group(ArgGroup::new("writer").required(false).multiple(false).args(&["write", "extend"])))]
pub struct RunArgs {
    #[arg(
        short = 'b',
        long = "bed",
        required = true,
        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Paths to BED12 files delimited by comma"
    )]
    pub bed: Vec<PathBuf>,

    #[arg(
        short = 'I',
        long = "index",
        help = "Create an index for the given bed file(s)",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub index: bool,

    #[arg(
        short = 'W',
        long = "write",
        help = "Writes the collapse counterpart of each bed file using the index (only if --index is set)",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub write: bool,

    #[arg(
        short = 'E',
        long = "extend",
        help = "Add the index (queue of collapse reads) as an additional column to the bed file(s)",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
        conflicts_with = "write"
    )]
    pub extend: bool,

    #[arg(
        short = 'o',
        long = "outdir",
        help = "Path to output directory (/collapse)",
        value_name = "PATH",
        required = false,
        default_value = "."
    )]
    pub outdir: PathBuf,
}

#[derive(Parser, Debug)]
pub struct ReadArgs {
    #[arg(short, long)]
    pub index: String,

    #[arg(short, long)]
    pub read: String,

    #[arg(short, long)]
    pub write: bool,
}
