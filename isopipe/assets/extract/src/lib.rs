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
//!
//! # Usage
//!
//! The `extract` command-line utility provide a simple invocation pattern:
//!
//! ```bash
//! extract [OPTIONS] --twobit <PATH> --bed <PATH>
//! ```
//!
//! ## Arguments
//!
//! * `-t`, `--twobit <PATH>`: Path to the 2bit genome file.
//! * `-b`, `--bed <PATH>`: Path to the input BED file.
//! * `-o`, `--outdir <PATH>`: Output directory for the extracted sequences. Defaults to the current directory (`.`).
//! * `-i`, `--index[=<FLAG>]`: Flag to return indexed sequences, resulting in a lighter FASTA file
//!     and an accompanying index. Defaults to `false`. Possible values are `true` or `false`.
//! * `-P`, `--prefix <VALUE>`: Prefix for the output directory. Defaults to `seqs`.
//! * `-S`, `--suffix <VALUE>`: Suffix for the output directory. Defaults to `tmp`.
//!     (Note: The original prompt had `-P, --preffix <VALUE>` for suffix, which is likely a
//!     typo and should be a different short flag like `-S` to avoid conflict with `--prefix`).
//! * `-c`, `--chunk-size <CHUNK_SIZE>`: Number of records to process in each chunk. Defaults to `30000`.
//! * `-T`, `--threads <THREADS>`: Specifies the number of threads to utilize
//!     for parallel processing. Defaults to the number of logical CPUs available.
//! * `-L`, `--level <LEVEL>`: Sets the logging level for output messages (e.g., `Info`, `Debug`, `Warn`, `Error`).
//!     Defaults to `INFO`.
//! * `-h`, `--help`: Prints help information.
//!
//! ## Example
//!
//! ```bash
//! extract \
//!   --twobit /path/to/genome.2bit \
//!   --bed /path/to/alignments.bed \
//!   --outdir /path/to/output_dir \
//! ```

#[cfg(feature = "core")]
pub mod cli;

#[cfg(feature = "core")]
pub mod ext;

#[cfg(feature = "read")]
pub mod read;
