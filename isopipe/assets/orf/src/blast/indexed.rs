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
    predictions: DashMap<String, BlastRecord>,
    records: DashMap<String, HashMap<String, GenePred>>,
    mut writer: BufWriter<File>,
) {
    let chr = get_chr_from_path(index);
    let index = extract::read::read_index(index, &chr);

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
        let n_nested = parts[2].split('@').nth(1).unwrap_or("");

        let (start, end, _) = if let Some(coords) = data.coords {
            coords
        } else {
            log::error!("ERROR: coords are 0! {seq_id} in {chr} -> {data:?}; {parts:?}");
            std::process::exit(1);
        };

        assert!(
            end > 0,
            "ERROR: orf_end coords are 0! {seq_id} in {chr} -> {data:?}; {parts:?}"
        );

        // INFO: unpacking index reference -> queries
        // INFO: for each query all blast records [ mapping 0 -> R12,R15,R20]
        let queries = index.get(&seq_id).unwrap_or_else(|| {
            panic!("ERROR: no queries found for ID: {seq_id} in {chr} -> {data:?}; {parts:?}\n{index:?}");
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
            let mut cannonical_id = format!("{}__OR{}", query, n_orf);
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
    });
}
