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

use dashmap::DashMap;
use hashbrown::HashMap;
use log::warn;
use packbed::record::GenePred;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::str::from_utf8;
use std::sync::Arc;

use crate::utils::*;

/// Inflates and processes indexed Blast prediction records, and writes them to a file.
///
/// This function is part of the `indexed` mode processing. It takes a map of
/// Blast predictions (where the keys are the sequential IDs from the indexing step)
/// and a binary index file. It "inflates" the results by mapping each unique
/// prediction back to all the original read IDs that shared that same sequence.
/// It also processes any indexed sequences that did not have a Blast hit, but
/// are associated with Translation AI (TAI) predictions (via the `helper` map),
/// creating a record for them with a "DM" tag (Denoted Missing). Finally, it
/// writes all of these inflated records to a specified output file.
///
/// # Arguments
///
/// * `index` - A `PathBuf` to the binary index file for the chunk.
/// * `predictions` - A `DashMap` where keys are the sequential IDs from the index
///                   and values are `BlastRecord`s.
/// * `records` - A `DashMap` containing the original `GenePred` records, keyed
///               by chromosome and then by read ID. Used to get CDS coordinates.
/// * `writer` - A `BufWriter<File>` to which the processed and inflated records are written.
/// * `helper` - A `HashMap` where keys are sequential IDs and values are `Vec<String>`
///              of TAI prediction targets (used for inflating TAI-only hits).
///
/// # Panics
///
/// This function will panic if:
/// - It cannot extract the chromosome name from the `index` file path.
/// - It fails to read the index file.
/// - A sequential ID from `predictions` is not found in the `index` map.
/// - It fails to write to the output `writer`.
/// - `get_cds_coords` or any parsing of TAI prediction strings fails.
pub fn indexed(
    index: &PathBuf,
    predictions: DashMap<String, Arc<BlastRecord>>,
    records: DashMap<String, HashMap<String, GenePred>>,
    mut writer: BufWriter<File>,
    idx_to_name: HashMap<usize, Vec<String>>,
    inner_idx_to_idxs: HashMap<u32, Vec<Arc<[u8]>>>,
) {
    let chr = get_chr_from_path(index);
    let mut index = extract::read::read_index(index, &chr);

    // INFO: inflate results!
    predictions.iter_mut().for_each(|mut record| {
        // INFO: id type -> 0_ORF.1_[1-10](+)@2
        let (id, data) = record.pair_mut();
        let parts = id.split('_').collect::<Vec<&str>>();
        let seq_id = parts[0]
            .parse::<u32>()
            .unwrap_or_else(|e| panic!("ERROR: {e}; could not parse id to u32 -> {id:?}"));
        let n_orf = parts[1]
            .strip_prefix("ORF.")
            .unwrap_or_else(|| panic!("ERROR: could not get ORF number from {parts:?}"));
        let n_nested = parts[2].split("@").last().unwrap_or("");

        let (start, end, _) = if let Some(coords) = data.coords {
            coords
        } else {
            (0, 0, '+')
        };

        // INFO: unpacking index reference -> queries
        // INFO: for each query all blast records [ mapping 0 -> R12,R15,R20]
        let queries = index.get(&seq_id).unwrap_or_else(|| {
            panic!("ERROR: no queries found for ID: {}", id);
        });

        for query in queries.into_iter() {
            let (orf_start, orf_end) = get_cds_coords(
                &records,
                &chr,
                &format!("{}", seq_id),
                start as u64,
                end as u64,
            );

            // WARN: skipping unreliable ORFs for the current alignment
            // INFO: none of these will match any other prediction because
            // INFO: the fall off any exonic boundary
            if orf_start == 0 && orf_end == 0 {
                log::debug!("WARN: ORF start and end are zero for ID: {}, skipping!", id);
                continue;
            }

            // WARN: query -> R{}_chr{} but
            let mut cannonical_id = format!("R{}_{}__OR{}", query, chr, n_orf);
            if !n_nested.is_empty() {
                cannonical_id += "#NE";
                cannonical_id += n_nested;
            }

            writer
                .write_all(
                    format!(
                        "{}\t{}\t{:.2}\t{:e}\t{}\t{}\t{:2}\t{}\t{}\t{}\t{}\n",
                        cannonical_id,
                        data.blast_idx_id,
                        data.blast_pid,
                        data.blast_e_value,
                        data.blast_offset,
                        data.blast_alignment_len,
                        data.percent_aligned,
                        start,
                        end,
                        orf_start,
                        orf_end
                    )
                    .as_bytes(),
                )
                .unwrap_or_else(|e| {
                    panic!("ERROR: failed to write blast record to file -> {e}");
                });
        }

        // INFO: removing the ID from the index to remain with unused IDs
        index.remove(&seq_id);
    });

    // INFO: repeating the process for unused ids
    // INFO: add tag DM to the ID -> identify unused ids
    // INFO: 0 -> {R}12,{R}15
    for (id, queries) in index.iter() {
        // INFO: each query mapping to that id needs to inflate lines with each target!
        // INFO: 0_ORF.87 [4632-4770](-) -> { 0 : [ 0_ORF.87_[4632-4770](-) ] }
        // INFO: 0_ORF.95 [4632-4770](-) -> { 0 : [ 0_ORF.87_[4632-4770](-), 0_ORF.95_[4632-4770](-) ] }
        // WARN: 0 -> { 0 : [ 0_ORF.87_[4632-4770](-), 0_ORF.95_[4632-4770](-), 0 ] } -> from double index!
        let Some(targets) = idx_to_name.get(&(*id as usize)) else {
            warn!("WARN: index number {id} with queries {queries:?} does not have any targets!");
            continue;
        };

        // INFO: looping on R12,R15
        for original in queries {
            // INFO: { 0 : [ 0_ORF.87_[4632-4770](-), 0_ORF.95@1_[4632-4770](-) ] }
            for target in targets {
                // INFO: checking if the current id was a secondary idx -> if can be parsed as u32
                let u_header = target.parse::<u32>();
                if let Ok(u_header) = u_header {
                    let handle = inner_idx_to_idxs.get(&u_header).unwrap_or_else(|| {
                        panic!("ERROR: could not find {u_header} in secondary index -> {target:?}!")
                    });

                    for record in handle {
                        let rc = from_utf8(record).unwrap().to_owned(); // INFO: safe to unwrap
                        map_target(original, &rc, &records, &chr, &mut writer);
                    }
                } else {
                    map_target(original, target, &records, &chr, &mut writer);
                };
            }
        }
    }
}

/// Parses a custom target string, retrieves corresponding gene data, and writes a
/// formatted record to a file buffer.
///
/// The `target` string is expected to be in a specific format, e.g.,
/// `0_ORF.87_[4632-4770](-)`. This function extracts the ID (`0`), ORF number (`87`),
/// genomic coordinates (`4632-4770`), and an optional nested ORF tag. It uses a
/// collection of gene records to get additional CDS coordinates before writing the
/// final record to the specified writer.
///
/// # Arguments
///
/// * `original` - The original record ID (e.g., from a query sequence).
/// * `target` - The formatted string containing the target gene information.
/// * `records` - A `DashMap` containing gene predictions, structured as
///   `DashMap<chromosome, HashMap<gene_id, GenePred>>`.
/// * `chr` - The chromosome name, used as a key for `records`.
/// * `writer` - A mutable reference to a `BufWriter` for writing the output.
///
/// # Panics
///
/// This function will panic if the `target` string does not conform to the expected
/// format, or if any of the following parsing or I/O operations fail:
/// - Splitting the string on `_`.
/// - Stripping the "ORF." prefix.
/// - Extracting the coordinate part within `[]`.
/// - Parsing coordinates into `u64`.
/// - Writing the final formatted string to the `writer`.
fn map_target(
    original: &String,
    target: &String,
    records: &DashMap<String, HashMap<String, GenePred>>,
    chr: &str,
    writer: &mut BufWriter<File>,
) {
    // INFO: target ->  0_ORF.87_[4632-4770](-)
    let parts: Vec<&str> = target.split('_').collect();

    let id = parts[0];
    let n_orf = parts[1]
        .strip_prefix("ORF.")
        .unwrap_or_else(|| panic!("ERROR: could not get ORF number from {parts:?}"));
    let n_nested = parts[2].split("@").last().unwrap_or("");

    let coords = parts[2]
        .strip_prefix('[')
        .unwrap_or_else(|| {
            panic!("ERROR: failed to get coords part from -> {target}");
        })
        .split(']')
        .next()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to get coords part from -> {target}");
        })
        .split('-')
        .collect::<Vec<&str>>(); // INFO: [4632, 4770]

    let start = coords[0].parse::<u64>().unwrap_or_else(|e| {
        panic!("ERROR: failed to parse start from {parts:?} -> {e}");
    });
    let end = coords[1].parse::<u64>().unwrap_or_else(|e| {
        panic!("ERROR: failed to parse end from {parts:?} -> {e}");
    });

    let (orf_start, orf_end) = get_cds_coords(records, chr, id, start, end);

    let mut cannonical_id = format!("R{}_{}__OR{}", original, chr, n_orf);
    if !n_nested.is_empty() {
        cannonical_id += "#NE";
        cannonical_id += n_nested;
    }

    writer
        .write_all(
            format!(
                "{}\t{}\t{:.2}\t{:e}\t{}\t{}\t{:2}\t{}\t{}\t{}\t{}\n",
                cannonical_id, id, 0.0, 1.0, 0, 0, 0.0, start, end, orf_start, orf_end
            )
            .as_bytes(),
        )
        .unwrap_or_else(|e| {
            panic!("ERROR: failed to write unused blast record to file -> {e}");
        });
}
