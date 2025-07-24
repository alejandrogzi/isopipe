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

use config::{BedColumn, BedColumnValue, bed_to_nested_collection};
use dashmap::DashMap;
use hashbrown::HashMap;
use log::warn;
use memchr::memchr;
use memmap2::Mmap;
use packbed::reader as bed_reader;
use smol_str::SmolStr;

use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::str::from_utf8;
use std::sync::Arc;

/// Parses a FASTA file into a `HashMap` where keys are sequence headers
/// and values are the sequences themselves (as `Vec<u8>`).
///
/// This function memory-maps the input FASTA file for efficient reading.
/// It iterates through the file, identifying each record by the `>` character,
/// extracts the header and sequence, and stores them in the `HashMap`.
/// Newline and carriage return characters are removed from the sequences.
///
/// # Type Parameters
///
/// * `F` - A type that can be converted into a `Path`, typically `&Path` or `PathBuf`.
///
/// # Arguments
///
/// * `f` - The path to the FASTA file.
///
/// # Returns
///
/// A `Result` which is:
/// - `Ok(HashMap<SmolStr, Vec<u8>>)`: A `HashMap` containing the parsed FASTA data.
///   `SmolStr` is used for efficient string storage of headers.
/// - `Err(Box<dyn std::error::Error>)`: If any error occurs during file opening,
///   memory mapping, or parsing.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to open the file.
/// - It fails to memory-map the file.
/// - It encounters a malformed FASTA entry where a newline character is not found
///   after a header.
/// - It fails to convert a header byte slice to a UTF-8 string.
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
/// use std::fs::File;
/// use std::io::Write;
/// use std::collections::HashMap;
/// use tempfile::NamedTempFile;
///
/// // Create a dummy FASTA file for the example
/// let mut temp_file = NamedTempFile::new().unwrap();
/// temp_file.write_all(b">seq1\nATGC\n>seq2\nTGCA\n").unwrap();
/// let path = temp_file.path().to_path_buf();
///
/// let result = parse_fa(&path);
/// assert!(result.is_ok());
/// let data = result.unwrap();
/// assert_eq!(data.len(), 2);
/// assert_eq!(data["seq1"], b"ATGC".to_vec());
/// assert_eq!(data["seq2"], b"TGCA".to_vec());
/// ```
pub fn parse_fa<F: AsRef<Path>>(
    f: F,
) -> Result<HashMap<SmolStr, Vec<u8>>, Box<dyn std::error::Error>> {
    let file = File::open(f).unwrap();
    let mmap = unsafe { Mmap::map(&file).unwrap() };
    let data = mmap.as_ref();

    let mut acc = HashMap::new();
    let mut pos = 0;

    while let Some(start) = memchr(b'>', &data[pos..]) {
        let start = pos + start;
        let end = memchr(b'>', &data[start + 1..]).map_or(data.len(), |e| start + 1 + e);
        let entry = &data[start + 1..end];
        let header_end = memchr(b'\n', entry).unwrap();
        let header = from_utf8(&entry[..header_end]).unwrap().trim();
        let record = &entry[header_end + 1..];
        let seq = record
            .iter()
            .filter(|&&b| b != b'\n' && b != b'\r') // Remove newlines and carriage returns
            .cloned()
            .collect::<Vec<u8>>();

        acc.insert(SmolStr::new(header), seq);
        pos = end;
    }

    Ok(acc)
}

/// Creates an empty index file with a `.dedup.index` extension for a given FASTA file.
///
/// This function constructs the output path for the index file by appending
/// `.dedup.index` to the base name of the input FASTA file. It then attempts
/// to create and open this file for writing, returning a `BufWriter` for buffered I/O.
///
/// # Arguments
///
/// * `fasta` - A `PathBuf` representing the path to the original FASTA file.
///
/// # Returns
///
/// A `BufWriter<File>` ready for writing to the newly created index file.
///
/// # Panics
///
/// This function will panic if it fails to create the index file.
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
/// use std::fs;
///
/// let fasta_path = PathBuf::from("my_sequences.fa");
/// // In a real scenario, you might create a dummy file first if `fasta_path` doesn't exist
/// // File::create(&fasta_path).unwrap();
///
/// let index_writer = create_index(&fasta_path);
/// // The index file "my_sequences.dedup.index" is now created and ready for writing.
///
/// // Clean up the dummy file for the example
/// let index_file_path = fasta_path.with_extension("dedup.index");
/// fs::remove_file(&index_file_path).unwrap();
/// ```
pub fn create_index(fasta: &PathBuf) -> BufWriter<File> {
    let index = fasta.with_extension("dedup.index");
    let writer = BufWriter::new(
        File::create(&index)
            .unwrap_or_else(|e| panic!("ERROR: cannot create index from sequences -> {e}")),
    );

    writer
}

/// Parses a BED file into a nested `DashMap` structure.
///
/// This function reads the specified BED file, converts its records into a
/// nested collection based on a primary grouping column (e.g., `BedColumn::Name`),
/// and extracts specified columns into `BedColumnValue` enums.
/// `DashMap` is used for concurrent access, suggesting this function might be
/// used in a parallel processing context.
///
/// # Type Parameters
///
/// * `K` - A type that implements the `config::BedParser` trait, which defines
///         how BED records are parsed into a specific structure.
///
/// # Arguments
///
/// * `bed` - A `PathBuf` representing the path to the input BED file.
/// * `columns` - A `Vec<BedColumn>` specifying which columns from the BED file
///               should be parsed and stored.
///
/// # Returns
///
/// A `DashMap<String, HashMap<String, Vec<BedColumnValue>>>` where:
/// - The outer `DashMap` keys are typically chromosome names or other primary identifiers.
/// - The inner `HashMap` keys are usually record names (from `BedColumn::Name`).
/// - The innermost `Vec<BedColumnValue>` contains the parsed values for the
///   specified `columns` for that record.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to read the BED file.
/// - It fails to convert the BED data into the specified nested collection structure.
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
/// use std::fs::File;
/// use std::io::Write;
/// use dashmap::DashMap;
/// use tempfile::NamedTempFile;
/// use crate::config::{self, BedColumn, BedColumnValue, BedParser};
///
/// // Assume Bed6 is a struct that implements BedParser
/// // and reader() and bed_to_nested_collection() are available.
///
/// // Create a dummy BED file for the example
/// let mut temp_file = NamedTempFile::new().unwrap();
/// temp_file.write_all(b"chr1\t10\t20\tgeneA\t0\t+\nchr1\t30\t40\tgeneB\t0\t-\n").unwrap();
/// let bed_path = temp_file.path().to_path_buf();
///
/// // Define the columns to parse
/// let columns_to_parse = vec![BedColumn::Chrom, BedColumn::Start, BedColumn::End, BedColumn::Name, BedColumn::Strand];
///
/// // Call parse_bed (assuming Bed6 is a valid BedParser implementation)
/// // let records = parse_bed::<Bed6>(&bed_path, columns_to_parse);
/// // assert!(!records.is_empty());
/// // assert!(records.contains_key("chr1"));
/// ```
pub fn parse_bed<K: config::BedParser>(
    bed: &PathBuf,
    columns: Vec<BedColumn>,
) -> DashMap<String, HashMap<String, Vec<BedColumnValue>>> {
    let rows = Arc::from(
        bed_reader(bed).unwrap_or_else(|e| panic!("ERROR: failed to read BED file -> {e}")),
    );
    let records = bed_to_nested_collection::<K>(rows, BedColumn::Name, columns)
        .unwrap_or_else(|e| panic!("ERROR: failed to convert BED to nested collection -> {e}"));

    records
}

/// Creates a new FASTA file with a specified extension.
///
/// This function constructs the output path by replacing the original
/// FASTA file's extension with the provided `extension`. If a file
/// with the new path does not already exist, it creates the file and
/// returns a `BufWriter` for buffered writing. If the file already exists,
/// it logs a warning and returns `None`.
///
/// # Arguments
///
/// * `fasta` - A `PathBuf` representing the path to the original FASTA file.
/// * `extension` - A string slice representing the new file extension
///                 (e.g., "dedup.fa", "output.fasta").
///
/// # Returns
///
/// An `Option<BufWriter<File>>`:
/// - `Some(BufWriter<File>)`: If the file was successfully created.
/// - `None`: If a file with the specified path and extension already exists.
///
/// # Panics
///
/// This function will panic if it attempts to create the file but fails
/// to do so due to underlying I/O errors (e.g., permissions issues, invalid path).
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
/// use std::fs;
/// use std::io::Write;
///
/// let original_fasta = PathBuf::from("sequences.fa");
/// // Create a dummy original file if it doesn't exist for the example
/// // File::create(&original_fasta).unwrap();
///
/// // Example 1: Creating a new file
/// let new_fasta_writer = create_fasta(&original_fasta, "dedup.fa");
/// assert!(new_fasta_writer.is_some());
/// let new_file_path = original_fasta.with_extension("dedup.fa");
/// assert!(new_file_path.exists());
/// fs::remove_file(&new_file_path).unwrap(); // Clean up
///
/// // Example 2: Attempting to create an existing file
/// File::create(&new_file_path).unwrap(); // Create it first
/// let existing_fasta_writer = create_fasta(&original_fasta, "dedup.fa");
/// assert!(existing_fasta_writer.is_none()); // Should return None because it exists
/// fs::remove_file(&new_file_path).unwrap(); // Clean up
/// ```
pub fn create_fasta(fasta: &PathBuf, extension: &str) -> Option<BufWriter<File>> {
    let file = fasta.with_extension(extension);
    if !file.exists() {
        let writer = BufWriter::new(
            File::create(&file).unwrap_or_else(|e| panic!("ERROR: cannot create file -> {e}")),
        );

        Some(writer)
    } else {
        warn!("WARN: file already exists -> {:?}!", file.display());
        None
    }
}
