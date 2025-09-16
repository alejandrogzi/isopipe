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

use config::bed_to_custom_struct_collection;
use dashmap::DashMap;
use hashbrown::HashMap;
use isopipe::config::depure;
use packbed::{reader as bed_reader, record::GenePred};

use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;
use std::sync::Arc;

use crate::{
    blast::{cannonical::*, core::*, indexed::*},
    cli::BlastArgs,
    utils::*,
};

pub mod cannonical;
pub mod core;
pub mod indexed;

const ORF_PEP: &str = "orfs.pep.fa";
const RESULT: &str = "dmd.result";
const HEADER_REGEX: &str = r"\[(\d+)-(\d+)\]\(([+-])\)";
pub const VENV: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tai/.venv/bin/activate");

/// Runs a complete BLAST analysis pipeline, including ORF prediction, deduplication,
/// and alignment against a DIAMOND database.
///
/// This function orchestrates the entire process:
/// 1. Predicts Open Reading Frames (ORFs) from the input FASTA file using `orfipy`.
/// 2. Deduplicates the predicted ORFs based on length and percentage.
/// 3. Performs a protein-protein BLAST search using DIAMOND against the specified database.
/// 4. Processes and inflates the BLAST results, writing them to an output file.
///
/// # Arguments
///
/// * `args` - A `BlastArgs` struct containing all the necessary arguments for the BLAST run,
///            including paths to input files, output directories, and various
///            parameters for `orfipy` and DIAMOND.
///
/// # Panics
///
/// This function will panic if any of the underlying commands (`orfipy`, `diamond`) fail
/// to execute, or if there are issues with file I/O or data parsing.
///
/// # Example
///
/// ```rust, no_run
/// run_blast(args);
/// ```
pub fn run_blast(args: BlastArgs) {
    let mode = Mode::from(&args.common.index);

    let dir = &args.common.outdir.join("orf");
    std::fs::create_dir_all(&dir)
        .unwrap_or_else(|e| panic!("ERROR: could not create directory -> {e}!"));

    let pep = orfipy(&args.common.fasta, &dir);

    let regex = regex::Regex::new(HEADER_REGEX).unwrap_or_else(|e| {
        panic!("ERROR: failed to compile regex for header parsing -> {e}");
    });

    let records = bed_to_custom_struct_collection::<GenePred>(
        bed_reader(&args.common.alignments)
            .unwrap_or_else(|e| panic!("ERROR: failed to read BED file -> {e}"))
            .into(),
        config::BedColumn::Name,
        config::BedOperation::SplitName(config::BIG_SEP, 0), // INFO: R9834_chr16__FC37#TC0#PA0#PR0#IY887 or 0
    )
    .unwrap_or_else(|e| panic!("ERROR: failed construct BED to GenePred collection -> {e}"));

    let (dedup, mut index, inner_idx_to_idxs) = deduplicate(
        &pep,
        true,
        args.orf_min_len,
        args.orf_min_percent,
        "M".as_bytes(),
        SeqType::AminoAcid,
        &mode,
        &regex,
    );

    match mode {
        Mode::Indexed => {
            index = args
                .common
                .index
                .unwrap_or_else(|| panic!("ERROR: could not unwrap index path, this is a bug!"))
        }
        _ => {}
    }

    diamond(
        &dedup,
        &args.database,
        &index,
        records,
        mode,
        &regex,
        inner_idx_to_idxs,
    );

    isopipe::depure!(dir, "result");
}

/// Predicts Open Reading Frames (ORFs) from a FASTA file using `orfipy`.
///
/// This function creates a temporary directory for `orfipy`'s output, constructs
/// and executes the `orfipy` command, and returns the path to the generated
/// protein FASTA file.
///
/// # Arguments
///
/// * `fasta` - A `PathBuf` representing the path to the input FASTA file.
/// * `outdir` - A `PathBuf` representing the base output directory where `orfipy`'s
///              temporary directory will be created.
/// * `executable` - A `PathBuf` representing the path to the `orfipy` executable.
///
/// # Returns
///
/// A `PathBuf` pointing to the generated protein FASTA file by `orfipy`.
///
/// # Panics
///
/// This function will panic if it cannot create the necessary output directory
/// or if the `orfipy` command fails to execute.
///
/// # Example
///
/// ```rust, no_run
/// let output = orfipy(&fasta_path, &output_dir, &orfipy_executable);
/// ```
fn orfipy(fasta: &PathBuf, dir: &PathBuf) -> PathBuf {
    let cmd = format!(
        "source {} && orfipy {} --pep {} --partial-5 --partial-3 --start ATG --include-stop --min 100 --ignore-case --outdir {}",
        VENV,
        fasta.display(),
        ORF_PEP,
        &dir.display()
    );

    std::process::Command::new("bash")
        .arg("-c")
        .arg(&cmd)
        .status()
        .unwrap_or_else(|e| panic!("ERROR: failed to execute orfipy command -> {e}"));

    // INFO: checking if orfipy output is empty and make it run again!
    let outfile = dir.join(ORF_PEP);
    if !outfile.exists() || outfile.metadata().unwrap().len() == 0 {
        // INFO: forcing orfipy to run again
        std::process::Command::new("bash")
            .arg("-c")
            .arg(cmd)
            .status()
            .unwrap_or_else(|e| panic!("ERROR: failed to execute orfipy command -> {e}"));
    }

    return outfile;
}

/// Performs a protein-protein BLAST search using DIAMOND and processes the results.
///
/// This function executes a DIAMOND `blastp` command, filters and processes the
/// resulting alignments, and inflates the results with additional information
/// from the original gene predictions. It then writes the final, processed
/// BLAST records to an output file. Unused IDs (those without a DIAMOND hit)
/// are also reported with a specific tag [#DM].
///
/// # Arguments
///
/// * `dedup` - A `PathBuf` representing the path to the deduplicated protein sequences
///             (query sequences for DIAMOND).
/// * `database` - A `PathBuf` representing the path to the DIAMOND database.
/// * `index` - A `PathBuf` representing the path to an index file containing
///             information about the original sequences, used for result inflation.
/// * `alignments` - A `PathBuf` representing the path to a BED file containing
///                  gene predictions, used to map absolute CDS coordinates.
///
/// # Panics
///
/// This function will panic if:
/// - The DIAMOND command fails to execute.
/// - It cannot read the BLAST predictions file or the BED alignment file.
/// - It cannot create the output file for BLAST results.
/// - It encounters issues parsing data from the BLAST output or the index.
/// - It cannot find a corresponding query ID or chromosome in the provided indices
///   or records during result inflation.
///
/// # Example
///
/// ```rust, no_run
/// diamond(&dedup_path, &database_path, &index_path, &alignments_path);
/// ```
fn diamond(
    dedup: &PathBuf,
    database: &PathBuf,
    index: &PathBuf,
    records: DashMap<String, HashMap<String, GenePred>>,
    mode: Mode,
    regex: &regex::Regex,
    inner_idx_to_idxs: HashMap<u32, Vec<Arc<[u8]>>>,
) {
    let diamond = dedup.with_extension("diamond");
    let cmd = format!(
        "diamond blastp --query {} --db {} --out {} --outfmt 6 qseqid pident qlen slen length qstart qend sstart send evalue --threads 8 --sensitive -e 1e-10",
        dedup.display(),
        database.display(),
        diamond.display()
    );

    std::process::Command::new("bash")
        .arg("-c")
        .arg(cmd)
        .status()
        .unwrap_or_else(|e| panic!("ERROR: failed to execute diamond command -> {e}"));

    let predictions = parse_predictions(&diamond, &mode, regex, &inner_idx_to_idxs);

    let writer = BufWriter::new(
        File::create(diamond.with_extension(RESULT)).unwrap_or_else(|e| {
            panic!("ERROR: failed to create output file for blast results -> {e}");
        }),
    );

    match mode {
        Mode::Raw => cannonical(index, predictions, records, writer),
        Mode::Indexed => indexed(index, predictions, records, writer),
    }
}
