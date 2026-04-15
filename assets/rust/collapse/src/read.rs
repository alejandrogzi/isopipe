// Copyright (c) 2026 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for collapsing BED files with deduplication and indexing
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for efficiently collapsing BED files
//! by identifying and deduplicating identical genomic intervals. The module provides
//! flexible output modes including collapsed BED files with queue annotations or
//! separate binary index files for space-efficient storage, extending the original
//! BED file with additional columns for read name and queue information and fast
//! lookups fo read names.
//!
//! In short, every BED entry is fingerprinted using byte-level hashing accounting
//! for specific columns from the original BED file to ensure that only truly identical
//! rows are grouped together. The deduplicated entries are held in memory alongside
//! their corresponding read identifier queues [maintaining original order for reconstruction].

use std::fmt::Debug;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;

use crate::utils::read_and_swap_index;

/// Looks up a read name in the index file and prints the corresponding queue.
///
/// If no read name is provided, the function will panic.
///
/// # Arguments
///
/// * `index` - A file path that implements `AsRef<Path>` pointing to the index file
/// * `read` - An optional read name to lookup in the index
///
/// # Panics
///
/// Panics if the read name is not provided or if it cannot be found in the index.
///
/// # Example
///
/// ```rust, ignore
/// let reconstructed_tracks = read_and_swap_index("genome.idx");
/// __lookup(&reconstructed_tracks, Some("read1"));
/// ```
pub fn __lookup<P: AsRef<Path> + Debug>(index: P, read: Option<String>) {
    let read = read.unwrap_or_else(|| panic!("ERROR: Read name not provided"));
    let map = read_and_swap_index(&index);

    let lookup = map
        .get(read.as_bytes())
        .unwrap_or_else(|| panic!("ERROR: Could not find read for {read} in {index:?}"));

    println!("{}", lookup);
}

/// Expands a read name in the index file and writes the corresponding queue to a file.
///
/// If no read name is provided, the function will panic.
///
/// # Arguments
///
/// * `index` - A file path that implements `AsRef<Path>` pointing to the index file
/// * `output` - A file path that implements `AsRef<Path>` where the expanded file will be written
///
/// # Panics
///
/// Panics if the read name is not provided or if it cannot be found in the index.
///
/// # Example
///
/// ```rust, ignore
/// let reconstructed_tracks = read_and_swap_index("genome.idx");
/// __expand(&reconstructed_tracks, "expanded.bed");
/// ```
pub fn __expand<P: AsRef<Path> + Debug>(index: P, output: P) {
    let mut writer = BufWriter::new(
        File::create(&output)
            .unwrap_or_else(|e| panic!("ERROR: Could not create file {:?} -> {e}", output)),
    );

    let map = read_and_swap_index(&index);

    for (_, queue) in map {
        let line = std::str::from_utf8(&queue.rep_line)
            .unwrap_or_else(|e| panic!("ERROR: Could not convert rep_line to utf8 -> {e}"));

        writeln!(writer, "{}", line)
            .unwrap_or_else(|e| panic!("ERROR: Could not write to expanded file -> {e}"));

        for read in &queue.reads {
            let read = std::str::from_utf8(read)
                .unwrap_or_else(|e| panic!("ERROR: Could not convert read to utf8 -> {e}"));

            let mut parts = line.split('\t').collect::<Vec<&str>>();
            parts[3] = read;

            writeln!(writer, "{}", parts.join("\t"))
                .unwrap_or_else(|e| panic!("ERROR: Could not write to expanded file -> {e}"));
        }
    }

    log::info!(
        "INFO: Wrote {} lines to expanded file",
        writer.buffer().lines().count()
    );
}
