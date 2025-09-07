//! Core module for detecting open reading frames in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for finding open reading frames (ORFs)
//! in a set of aligned reads.
//!
//! In short, every possible open reading frame (ORF) is detected for every
//! read in the query set. For every potential ORF, learning models and databases
//! are used to determine whether the ORF is a true ORF, a false positive.
//! All the data from each reliable ORF is collected and subjected to another
//! learning model trained with true ORFs and false positives. The process is
//! heavily parallelized to offer fast performance on large datasets.
//!
//! # Usage
//!
//! The `orf` command-line utility provides several subcommands to orchestrate
//! the ORF detection pipeline. The general invocation pattern is:
//!
//! ```bash
//! orf [OPTIONS] <COMMAND>
//! ```
//!
//! Global options applicable to all commands include:
//!
//! * `-t`, `--threads <THREADS>`: Specifies the number of threads to utilize for parallel processing.
//!   Defaults to the number of logical CPUs available.
//! * `-L`, `--level <LEVEL>`: Sets the logging level for output messages (e.g., `Info`, `Debug`, `Warn`, `Error`).
//!   Defaults to `Info`.
//!
//! ## Commands
//!
//! ### `blast`
//!
//! Executes ORF detection using `orfipy` for candidate ORF identification and `diamond` for homology-based filtering.
//!
//! #### Arguments
//!
//! * `--fasta <PATH>`: Required. Path to the input FASTA file containing query sequences.
//! * `--alignments <PATH>`: Required. Path to the input BED file containing alignment information for the query sequences.
//! * `--outdir <DIR>`: Output directory for results. Defaults to the current directory (`.`).
//! * `-e`, `--executable <PATH>`: Path to the `orfipy` executable. Defaults to `orfipy`.
//! * `-d`, `--db <PATH>`: Required. Path to the DIAMOND BLAST database.
//! * `-l`, `--orf-min-len <LENGTH>`: Minimum length for an ORF. Defaults to `50`.
//! * `-p`, `--orf-min-percent <PERCENT>`: Minimum percentage of an ORF's length relative to the full sequence. Defaults to `0.25`.
//! * `-P`, `--pattern <PATTERN>`: Pattern for subsequence (amino acid `aa` or nucleotide `nt`). Defaults to `M`.
//!
//! #### Example
//!
//! ```bash
//! orf blast \
//!   --fasta /path/to/reads.fasta \
//!   --alignments /path/to/alignments.bed \
//!   --outdir /path/to/blast_results \
//!   --db /path/to/uniprot.dmnd \
//!   --orf-min-len 100
//! ```
//!
//! ### `tai`
//!
//! Performs ORF detection utilizing the Translation AI (TAI) model.
//!
//! #### Arguments
//!
//! * `--fasta <PATH>`: Required. Path to the input FASTA file.
//! * `--alignments <PATH>`: Required. Path to the input BED file.
//! * `--outdir <DIR>`: Output directory for results. Defaults to the current directory (`.`).
//! * `-t`, `--threshold <THRESHOLD>`: Prediction threshold for TranslationAI. Defaults to `0.01`.
//!
//! #### Example
//!
//! ```bash
//! orf tai \
//!   --fasta /path/to/reads.fasta \
//!   --alignments /path/to/alignments.bed \
//!   --outdir /path/to/tai_results \
//!   --threshold 0.05
//! ```
//!
//! ### `merge`
//!
//! Consolidates the results obtained from the `blast` and `tai` commands, optionally incorporating TOGA results.
//!
//! #### Arguments
//!
//! * `-b`, `--blast <PATH>`: Required. Path to the BLAST results file.
//! * `-T`, `--tai <PATH>`: Required. Path to the TranslationAI results file.
//! * `-t`, `--toga <PATH>`: Required. Path to the TOGA merged results file.
//! * `-a`, `--alignments <PATH>`: Required. Path to the original alignment BED file.
//! * `-o`, `--outdir <DIR>`: Output directory for the merged results. Defaults to the current directory (`.`).
//!
//! #### Example
//!
//! ```bash
//! orf merge \
//!   --blast /path/to/blast_results/blast.out \
//!   --tai /path/to/tai_results/tai.out \
//!   --toga /path/to/toga_results/toga_merged.tsv \
//!   --alignments /path/to/alignments.bed \
//!   --outdir /path/to/final_results
//! ```
//!
//! ### `toga`
//!
//! Processes and merges results from the TOGA (Tool to infer Orthologs from Genome Alignments) pipeline.
//!
//! #### Arguments
//!
//! * `-p`, `--path <PATH>`: Required. Path to the TOGA results directory.
//! * `-o`, `--outdir <DIR>`: Output directory for the merged TOGA results. Defaults to the current directory (`.`).
//!
//! #### Example
//!
//! ```bash
//! orf toga \
//!   --path /path/to/toga_raw_output \
//!   --outdir /path/to/processed_toga
//! ```

use clap::{self, Parser};
use log::info;
use simple_logger::init_with_level;

use orf::{
    blast::run_blast,
    cli::{Args, Commands},
    merge::merge,
    tai::run_tai,
    toga::run_toga,
};

fn main() {
    let start = std::time::Instant::now();

    let args = Args::parse();

    init_with_level(args.level).unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();

    match args.command {
        Commands::Blast(args) => run_blast(args),
        Commands::Tai(args) => run_tai(args),
        Commands::Toga(args) => run_toga(args),
        Commands::Merge(args) => merge(args),
    }

    let elapsed = start.elapsed();
    info!("Elapsed time: {:.3?}", elapsed);
}
