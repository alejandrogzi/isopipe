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

use clap::{ArgAction, Parser};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "orf", about = "Extract sequences from a .2bit using .bed", version = env!("CARGO_PKG_VERSION"))]
pub struct Args {
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
        short = 'i',
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
        short = 'P',
        long = "preffix",
        required = false,
        help = "Outdir suffix",
        value_name = "VALUE",
        default_value("tmp")
    )]
    pub suffix: String,

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
