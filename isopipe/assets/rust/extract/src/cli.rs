//! Core module for extracting sequences using .2bit from a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for efficiently getting whole
//! sequences from .2bit genomic files using plain coordinates. Additionally,
//! it provides the option of deduplicate repeated entries in order to save
//! space.
//!
//! In short, every sequence is extracted from a .2bit and holded in memory
//! as plain bytes for every read in the query set [rev comp for neg strand].
//! The command line interface provides the option of specifying indexing, leading
//! to a one-pass deduplication step and the creation of an index where
//! simple integers map to read identifiers [all of them as plain bytes].
//! The process is heavily parallelized to offer fast performance on large datasets.

use clap::{ArgAction, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "extract", about = "Extract sequences from a .2bit using .bed [and index them]", version = env!("CARGO_PKG_VERSION"))]
pub struct Args {
    #[command(subcommand)]
    pub command: Commands,

    #[arg(
        short = 'T',
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
    /// Run extract using .2bit + .bed
    Base(ExtractArgs),

    /// Reads an index created with the 'base'
    Read(IndexArgs),
}

#[derive(Debug, Parser)]
pub struct ExtractArgs {
    #[arg(
        short = 't',
        long = "twobit",
        required = true,
        help = "Path to fasta file",
        value_name = "PATH"
    )]
    pub twobit: PathBuf,

    #[arg(
        short = 'b',
        required = true,
        long = "bed",
        help = "Path to bed file",
        value_name = "PATH"
    )]
    pub bed: PathBuf,

    #[arg(
        short = 'o',
        long = "outdir",
        help = "Output directory",
        value_name = "PATH",
        default_value(".")
    )]
    pub output_dir: PathBuf,

    #[arg(
        short = 'I',
        long = "index",
        help = "Flag to return indexed sequences [light .fa + index]",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub mode: bool,

    #[arg(
        short = 'P',
        long = "prefix",
        required = false,
        help = "Outdir preffix",
        value_name = "VALUE",
        default_value("seqs")
    )]
    pub dir_prefix: String,

    #[arg(
        short = 'S',
        long = "suffix",
        required = false,
        help = "Outdir suffix",
        value_name = "VALUE",
        default_value("tmp")
    )]
    pub suffix: String,

    #[arg(
        short = 's',
        long = "sequence-mode",
        required = false,
        help = "Sequence extraction mode [whole sequence, exon, intron, cds]",
        value_name = "MODE",
        default_value("exon")
    )]
    pub seq_mode: SeqMode,

    #[arg(
        short = 'c',
        long = "chunk-size",
        required = false,
        help = "Chunk size",
        value_name = "CHUNK_SIZE",
        default_value("30000")
    )]
    pub chunk_size: usize,

    #[arg(
        short = 'J',
        long = "join",
        help = "Flag to join all chunked .fa files",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub join: bool,

    #[arg(
        short = 'X',
        long = "translate",
        help = "Flag to translate the sequences to amino acids",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub translate: bool,

    #[arg(
        short = 'N',
        long = "no-reduced-bed",
        help = "Flag to not create a reduced .bed file",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub no_reduced_bed: bool,

    #[arg(
        short = 'u',
        long = "flank-upstream",
        help = "Flank upstream in bp",
        value_name = "VALUE",
        default_value("0")
    )]
    pub flank_upstream: usize,

    #[arg(
        short = 'd',
        long = "flank-downstream",
        help = "Flank downstream in bp",
        value_name = "VALUE",
        default_value("0")
    )]
    pub flank_downstream: usize,

    #[arg(
        short = 'K',
        long = "split-extraction",
        help = "Split extracted sequences as individual components (each exon or each intron)",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub split_extraction: bool,

    #[arg(
        short = 'Z',
        long = "icc",
        help = "Output intron IC sequences in the format of name-flank-seq-flank, tab-separated",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
        requires = "split-extraction"
    )]
    pub intron_ic_output_fmt: bool,
}

#[derive(ValueEnum, Clone, Debug)]
pub enum SeqMode {
    #[value(name = "genome", alias = "g")]
    Genome,
    #[value(name = "exon", alias = "e")]
    Exon,
    #[value(name = "intron", alias = "i")]
    Intron,
    #[value(name = "cds", alias = "c")]
    CDS,
}

#[derive(Debug, Parser)]
pub struct IndexArgs {
    #[arg(
        short = 'i',
        long = "index",
        required = true,
        help = "Path to .index file produced by extract step",
        value_name = "PATH"
    )]
    pub index: PathBuf,

    #[arg(
        short = 'w',
        long = "write",
        help = "Flag to write the whole index as plain text",
        action = ArgAction::SetTrue
    )]
    pub write: bool,

    #[arg(
        long = "id",
        help = "Index ID for lookup",
        value_name = "VALUE",
        conflicts_with = "write"
    )]
    pub id: Option<u32>,

    #[arg(
        short = 'o',
        long = "outdir",
        help = "Output directory",
        value_name = "PATH",
        default_value("index")
    )]
    pub output_dir: PathBuf,
}
