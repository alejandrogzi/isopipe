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
    predictions: DashMap<u32, BlastRecord>,
    records: DashMap<String, HashMap<String, GenePred>>,
    mut writer: BufWriter<File>,
    helper: HashMap<usize, Vec<String>>,
) {
    let chr = get_chr_from_path(index);
    let mut index = extract::read::read_index(index, &chr);

    // INFO: inflate results!
    predictions.iter_mut().for_each(|mut record| {
        let (id, data) = record.pair_mut();
        let (start, end, _) = if let Some(coords) = data.coords {
            coords
        } else {
            (0, 0, '+')
        };

        // INFO: unpacking index reference -> queries
        // INFO: for each query all blast records
        // INFO: { index_id : [(read_id: u16, chr_bytes: [u8; chr_len], orf: u16, subseq_orf: u16, seq_len: usize)] }
        // INFO: { 0 : [(read_id: 5903, chr_bytes: [16, 32], orf: 1, subseq_orf: 3, seq_len: 350)] }
        let queries = index.get(id).unwrap_or_else(|| {
            panic!("ERROR: no queries found for ID: {}", id);
        });

        for query in queries.into_iter() {
            let (orf_start, orf_end) =
                get_cds_coords(&records, &chr, query, start as u64, end as u64);

            // WARN: skipping unreliable ORFs for the current alignment
            // INFO: none of these will match any other prediction because
            // INFO: the fall off any exonic boundary
            if orf_start == 0 && orf_end == 0 {
                // warn!(
                //     "WARN: ORF start and end are zero for ID: {}, skipping!",
                //     query_id
                // );
                continue;
            }

            // WARN: query -> R{}_chr{} but
            data.set_id(query.clone());

            writer
                .write_all(
                    format!(
                        "{}\t{}\t{:.2}\t{:e}\t{}\t{}\t{:2}\t{}\t{}\t{}\t{}\n",
                        data.blast_id,
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
        index.remove(id);
    });

    // INFO: repeating the process for unused ids
    // INFO: add tag DM to the ID -> identify unused ids
    for (id, queries) in index.iter() {
        let Some(targets) = helper.get(&(*id as usize)) else {
            warn!("WARN: index number {id} with queries {queries:?} does not have any targets!");
            continue;
        };

        // INFO: each query mapping to that id needs to inflate lines with each target!
        for query in queries {
            // INFO: { 0 : [ 0_ORF.87_[4632-4770](-), 0_ORF.95@1_[4632-4770](-) ] }
            for target in targets {
                let parts = target
                    .split('[')
                    .last()
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

                let start = parts[0].parse::<u64>().unwrap_or_else(|e| {
                    panic!("ERROR: failed to parse start from {parts:?} -> {e}");
                });
                let end = parts[1].parse::<u64>().unwrap_or_else(|e| {
                    panic!("ERROR: failed to parse end from {parts:?} -> {e}");
                });

                let (orf_start, orf_end) = get_cds_coords(&records, &chr, query, start, end);

                // WARN: skipping unreliable ORFs for the current alignment
                // INFO: none of these will match any other prediction because
                // INFO: the fall off any exonic boundary
                if orf_start == 0 && orf_end == 0 {
                    warn!(
                        "WARN: ORF start and end are zero for ID: {}, skipping!",
                        query
                    );
                    continue;
                }

                writer
                    .write_all(
                        format!(
                            "{}\t{}\t{:.2}\t{:e}\t{}\t{}\t{:2}\t{}\t{}\t{}\t{}\n",
                            query, id, 0.0, 1.0, 0, 0, 0.0, start, end, orf_start, orf_end
                        )
                        .as_bytes(),
                    )
                    .unwrap_or_else(|e| {
                        panic!("ERROR: failed to write unused blast record to file -> {e}");
                    });
            }
        }
    }
}
