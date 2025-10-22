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

use hashbrown::{HashMap, hash_map::Entry};

use std::path::PathBuf;
use std::sync::Arc;

use crate::utils::*;

/// Deduplicates sequences in a FASTA file based on exact sequence matches,
/// and optionally performs subsequence nesting by splitting records at a
/// specified start signal (e.g., 'M' for amino acids, 'ATG' for nucleotides).
///
/// This function reads sequences from the input FASTA, stores them in a HashMap
/// for deduplication, and writes the unique (and potentially nested) sequences
/// to a new FASTA file. It also generates an index file mapping original
/// sequence IDs to their deduplicated counterparts.
///
/// # Arguments
///
/// * `fasta` - A `PathBuf` representing the path to the input FASTA file.
/// * `do_nesting` - A boolean flag indicating whether to perform subsequence nesting.
/// * `min_len` - The minimum length for a sequence or subsequence to be considered.
/// * `min_percent` - The minimum percentage of the original sequence length
///                   that a subsequence must represent to be considered.
/// * `pattern` - A byte slice representing the start signal for splitting sequences
///               during nesting (e.g., `b"M"` for methionine, `b"ATG"` for a start codon).
/// * `seq_type` - A `SeqType` enum indicating whether the sequences are nucleotides
///                or amino acids, which affects the splitting logic.
///
/// # Returns
///
/// A tuple containing two `PathBuf` instances:
/// - The path to the deduplicated FASTA file.
/// - The path to the generated index file.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to parse the input FASTA file.
/// - It cannot create the output deduplicated FASTA file or the index file.
/// - The regular expression for header parsing fails to compile.
/// - Errors occur during the indexing process.
pub fn deduplicate(
    fasta: &PathBuf,
    do_nesting: bool,
    min_len: usize,
    min_percent: f32,
    pattern: &[u8],
) -> HashMap<HashHead, Vec<OrfRecord>> {
    let seqs =
        parse_fa(fasta).unwrap_or_else(|e| panic!("ERROR: failed to parse FASTA file -> {e}"));

    // INFO  HH<idx  seq>         Vec<OrfRecord>
    // INFO: { (1, MADAD*) : [0.ORF.1, 1.ORF.1@1] }
    let mut mapper: HashMap<HashHead, Vec<OrfRecord>> = HashMap::new(); // INFO: seq -> name
    let mut idx = 0; // INFO: stores deduplicated index of sequence

    // INFO: loops through sequences and populates mapper
    for (mut header, seq) in seqs.into_iter() {
        // INFO: sequence used as key in mapper -> hhead
        let len = seq.len();
        let mut key = Vec::with_capacity(len);
        for b in &seq {
            if *b != b'\n' {
                key.push(*b);
            }
        }

        let hhead = HashHead::new(Arc::from(key));

        // INFO: if hhead exists -> header inherits seq_idx
        match mapper.entry(hhead) {
            Entry::Occupied(mut entry) => {
                // INFO: seq exists -> push header to vec
                // INFO: seq_idx is inherited from first header
                entry.get_mut().push(header.clone());
            }
            Entry::Vacant(entry) => {
                header.seq_idx = idx;
                entry.insert(vec![header.clone()]);
                idx += 1; // INFO: increment index count
            }
        }

        if do_nesting {
            __split_record(
                &header,
                &seq,
                len,
                min_len,
                min_percent,
                &mut mapper,
                pattern,
                &mut idx,
            )
        }
    }

    mapper
}
