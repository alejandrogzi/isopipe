use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "orf", about = "Open reading frame pipeline wrapper", version = env!("CARGO_PKG_VERSION"))]
pub struct Args {
    #[command(subcommand)]
    pub command: Commands,

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

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Run ORF detection using orfipy + diamond
    Blast(BlastArgs),

    /// Run ORF detection using TranslationAi
    Tai(TaiArgs),

    /// Merge diamond and translationAi results
    Merge(MergeArgs),

    /// Read and merge TOGA results
    Toga(TogaArgs),
}

#[derive(Debug, Parser)]
pub struct CommonArgs {
    #[arg(
        short = 'f',
        long = "fasta",
        required = true,
        help = "Path to .fa file"
    )]
    pub fasta: PathBuf,

    #[arg(
        short = 'a',
        long = "alignments",
        required = true,
        help = "Path to .bed file"
    )]
    pub alignments: PathBuf,

    #[arg(
        short = 'o',
        long = "outdir",
        default_value = ".",
        help = "Output directory"
    )]
    pub outdir: PathBuf,
}

#[derive(Debug, Parser)]
pub struct BlastArgs {
    #[clap(flatten)]
    pub common: CommonArgs,

    #[arg(
        short = 'e',
        long = "executable",
        default_value = "orfipy",
        help = "Path to orfipy executable"
    )]
    pub orfipy: PathBuf,

    #[arg(
        short = 'd',
        long = "db",
        required = true,
        help = "Path to BLAST database"
    )]
    pub database: PathBuf,

    #[arg(
        long = "orf-min-len",
        short = 'l',
        default_value = "50",
        help = "Minimum ORF length"
    )]
    pub orf_min_len: usize,

    #[arg(
        short = 'p',
        long = "orf-min-percent",
        default_value = "0.25",
        help = "Minimum ORF percent"
    )]
    pub orf_min_percent: f32,

    #[arg(
        short = 'P',
        long = "pattern",
        default_value = "M",
        help = "Pattern to subsequence [aa/nt]"
    )]
    pub pattern: String,
}

#[derive(Debug, Parser)]
pub struct TaiArgs {
    #[clap(flatten)]
    pub common: CommonArgs,

    #[arg(
        short = 't',
        long = "threshold",
        default_value = "0.01",
        help = "TranslationAI threshold"
    )]
    pub threshold: String,
}

#[derive(Debug, Parser)]
pub struct MergeArgs {
    #[arg(
        short = 'b',
        long = "blast",
        required = true,
        help = "Path to BLAST results file"
    )]
    pub blast: PathBuf,

    #[arg(
        short = 'T',
        long = "tai",
        required = true,
        help = "Path to TranslationAI results file"
    )]
    pub tai: PathBuf,

    #[arg(
        short = 't',
        long = "toga",
        required = true,
        help = "Path to TOGA merged file"
    )]
    pub toga: PathBuf,

    #[arg(
        short = 'a',
        long = "alignments",
        required = true,
        help = "Path to .bed file"
    )]
    pub alignments: PathBuf,

    #[arg(
        short = 'o',
        long = "outdir",
        default_value = ".",
        help = "Output directory"
    )]
    pub outdir: PathBuf,
}

#[derive(Debug, Parser)]
pub struct TogaArgs {
    #[arg(
        short = 'p',
        long = "path",
        required = true,
        help = "Path to TOGA results directory"
    )]
    pub path: PathBuf,

    #[arg(
        short = 'o',
        long = "outdir",
        default_value = ".",
        help = "Output directory"
    )]
    pub outdir: PathBuf,
}
