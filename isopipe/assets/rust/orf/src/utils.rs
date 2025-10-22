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
use hashbrown::{HashMap, hash_map::Entry};
use log::warn;
use memchr::memchr;
use memmap2::Mmap;
use packbed::{reader as bed_reader, record::GenePred};

use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::str::{FromStr, from_utf8};
use std::sync::Arc;

pub const CODON_TABLE: [(&[u8], &str); 64] = [
    // Phenylalanine (F)
    (b"TTT", "F"),
    (b"TTC", "F"),
    // Leucine (L)
    (b"TTA", "L"),
    (b"TTG", "L"),
    (b"CTA", "L"),
    (b"CTC", "L"),
    (b"CTG", "L"),
    (b"CTT", "L"),
    // Isoleucine (I)
    (b"ATT", "I"),
    (b"ATC", "I"),
    (b"ATA", "I"),
    // Methionine (M) - Start codon
    (b"ATG", "M"),
    // Valine (V)
    (b"GTA", "V"),
    (b"GTC", "V"),
    (b"GTG", "V"),
    (b"GTT", "V"),
    // Serine (S)
    (b"TCA", "S"),
    (b"TCC", "S"),
    (b"TCG", "S"),
    (b"TCT", "S"),
    (b"AGT", "S"),
    (b"AGC", "S"),
    // Proline (P)
    (b"CCA", "P"),
    (b"CCC", "P"),
    (b"CCG", "P"),
    (b"CCT", "P"),
    // Threonine (T)
    (b"ACA", "T"),
    (b"ACC", "T"),
    (b"ACG", "T"),
    (b"ACT", "T"),
    // Alanine (A)
    (b"GCA", "A"),
    (b"GCC", "A"),
    (b"GCG", "A"),
    (b"GCT", "A"),
    // Tyrosine (Y)
    (b"TAT", "Y"),
    (b"TAC", "Y"),
    // Stop codons (*)
    (b"TAA", "*"),
    (b"TAG", "*"),
    (b"TGA", "*"),
    // Histidine (H)
    (b"CAT", "H"),
    (b"CAC", "H"),
    // Glutamine (Q)
    (b"CAA", "Q"),
    (b"CAG", "Q"),
    // Asparagine (N)
    (b"AAT", "N"),
    (b"AAC", "N"),
    // Lysine (K)
    (b"AAA", "K"),
    (b"AAG", "K"),
    // Aspartic acid (D)
    (b"GAT", "D"),
    (b"GAC", "D"),
    // Glutamic acid (E)
    (b"GAA", "E"),
    (b"GAG", "E"),
    // Cysteine (C)
    (b"TGT", "C"),
    (b"TGC", "C"),
    // Tryptophan (W)
    (b"TGG", "W"),
    // Arginine (R)
    (b"CGA", "R"),
    (b"CGC", "R"),
    (b"CGG", "R"),
    (b"CGT", "R"),
    (b"AGA", "R"),
    (b"AGG", "R"),
    // Glycine (G)
    (b"GGA", "G"),
    (b"GGC", "G"),
    (b"GGG", "G"),
    (b"GGT", "G"),
];

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
) -> Result<HashMap<OrfRecord, Vec<u8>>, Box<dyn std::error::Error>> {
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
            .filter(|&&b| b != b'\n' && b != b'\r') // INFO: remove newlines and carriage returns
            .cloned()
            .collect::<Vec<u8>>();

        acc.insert(OrfRecord::parse(header), seq);
        pos = end;
    }

    Ok(acc)
}

/// An enum representing the type of an ORF
#[derive(Debug, Clone, Copy, Eq, Hash, PartialEq)]
pub enum OrfType {
    Complete,
    ThreePartial,
    FivePartial,
    FivePartialNested,
    ThreePartialNested,
    CompleteNested,
}

impl std::str::FromStr for OrfType {
    type Err = Box<dyn std::error::Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "type:complete" => Ok(OrfType::Complete),
            "type:3-prime-partial" => Ok(OrfType::ThreePartial),
            "type:5-prime-partial" => Ok(OrfType::FivePartial),
            "type:5-prime-partial-nested" => Ok(OrfType::FivePartialNested),
            "type:3-prime-partial-nested" => Ok(OrfType::ThreePartialNested),
            "type:complete-nested" => Ok(OrfType::CompleteNested),
            _ => Err(Box::new(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Invalid ORF type: {}", s),
            ))),
        }
    }
}

impl std::fmt::Display for OrfType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OrfType::Complete => write!(f, "CO"),
            OrfType::ThreePartial => write!(f, "TP"),
            OrfType::FivePartial => write!(f, "FP"),
            OrfType::FivePartialNested => write!(f, "FN"),
            OrfType::ThreePartialNested => write!(f, "TN"),
            OrfType::CompleteNested => write!(f, "CN"),
        }
    }
}

/// Represents a single ORF record from orfipy
///
/// # Fields
///
/// * `id` - The ID of the ORF
/// * `start` - The start position of the ORF
/// * `end` - The end position of the ORF
/// * `orf_type` - The type of the ORF
/// * `frame` - The frame of the ORF
/// * `start_codon` - The start codon of the ORF
/// * `stop_codon` - The stop codon of the ORF
/// * `strand` - The strand of the ORF
#[derive(Debug, Clone, Eq, Hash, PartialEq)]
pub struct OrfRecord {
    pub id: String,
    pub start: usize,
    pub end: usize,
    pub orf_type: OrfType,
    pub frame: i8,
    pub start_codon: String,
    pub stop_codon: String,
    pub strand: char,
    pub seq_idx: usize,
}

/// Splits a header string into its components
///
/// # Arguments
///
/// * `header` - The header string
///
/// # Returns
///
/// * `Option<(usize, usize, char)>` - The start, end, and strand of the ORF
///
/// # Example
///
/// ```rust
/// let (start, end, strand) = split_coords(coords).unwrap_or_else(|| {
///     panic!("ERROR: failed to parse coords from line: {}", line);
/// });
/// ```
pub fn split_coords(header: &str) -> Option<(usize, usize, char)> {
    let mut parts = header.split('-');

    let start = parts
        .next()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse ORF start from header: {}", header);
        })
        .strip_prefix('[')
        .unwrap_or_else(|| {
            panic!(
                "ERROR: failed to strip prefix from ORF start from header: {}",
                header
            );
        })
        .parse::<usize>()
        .unwrap_or_else(|e| {
            panic!(
                "ERROR: failed to parse ORF start from header: {} -> {e}",
                header
            );
        });

    parts = parts
        .next()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse ORF end from header: {}", header);
        })
        .split(']');

    let end = parts
        .next()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse ORF end from header: {}", header);
        })
        .parse::<usize>()
        .unwrap_or_else(|e| {
            panic!(
                "ERROR: failed to parse ORF end from header: {} -> {e}",
                header
            );
        });

    let strand = parts
        .next()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse ORF strand from header: {}", header);
        })
        .char_indices()
        .nth(1)
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse ORF strand from header: {}", header);
        })
        .1;

    Some((start, end, strand))
}

impl OrfRecord {
    /// Parses an ORF record from a string
    ///
    /// # Arguments
    ///
    /// * `line` - The string to parse
    ///
    /// # Returns
    ///
    /// * `Self` - The parsed ORF record
    ///
    /// # Example
    ///
    /// ```rust
    /// let record = OrfRecord::parse(line);
    /// ```
    pub fn parse(line: &str) -> Self {
        let mut parts = line.split(' ');

        // WARN: if headers have spaces splitting will fail!
        // INFO: >ENSMUST00000167914.2_ORF.1 [69-330](+) type:complete length:258 frame:1 start:ATG stop:TGA
        let (id, coords, orf_type, _, frame, start_codon, stop_codon) = (
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse ORF ID from line: {}", line);
                })
                .to_string(),
            parts.next().unwrap_or_else(|| {
                panic!(
                    "ERROR: failed to parse ORF coords from line: {} -> {:?}",
                    line, parts
                );
            }),
            OrfType::from_str(parts.next().unwrap_or_else(|| {
                panic!("ERROR: failed to parse ORF type from line: {}", line);
            }))
            .unwrap_or_else(|e| {
                panic!("ERROR: failed to parse ORF type from line: {} -> {e}", line);
            }),
            parts.next(), // INFO: length
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse ORF frame from line: {}", line);
                })
                .strip_prefix("frame:")
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to strip prefix from ORF frame from line: {}",
                        line
                    );
                })
                .parse::<i8>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse ORF frame from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse ORF start codon from line: {}", line);
                })
                .strip_prefix("start:")
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to strip prefix from ORF start codon from line: {}",
                        line
                    );
                })
                .to_string(),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse ORF stop codon from line: {}", line);
                })
                .strip_prefix("stop:")
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to strip prefix from ORF stop codon from line: {}",
                        line
                    );
                })
                .to_string(),
        );

        let (start, end, strand) = split_coords(coords).unwrap_or_else(|| {
            panic!("ERROR: failed to parse coords from line: {}", line);
        });

        Self {
            id,
            start,
            end,
            orf_type,
            frame,
            start_codon,
            stop_codon,
            strand,
            seq_idx: 0, // INFO: set on the fly
        }
    }
}

/// Splits an ORF record into nested ORF records
///
/// # Arguments
///
/// * `header` - The header of the ORF record
/// * `seq` - The sequence of the ORF record
/// * `seq_length` - The length of the sequence
/// * `min_len` - The minimum length of the ORF
/// * `min_percent` - The minimum percentage of the ORF
/// * `mapper` - The mapper to store the nested ORF records
/// * `needle` - The needle to match for the nested ORF
/// * `idx_count` - The index count to increment
///
/// # Returns
///
/// None
///
/// # Example
///
/// ```rust
/// let mapper = HashMap::new();
/// let needle = b"M";
/// let idx_count = 0;
/// __split_record(&header, &seq, seq_length, min_len, min_percent, &mut mapper, &needle, &mut idx_count);
/// ```
#[allow(clippy::too_many_arguments)]
pub fn __split_record(
    header: &OrfRecord,
    seq: &[u8],
    seq_length: usize,
    min_len: usize,
    min_percent: f32,
    mapper: &mut HashMap<HashHead, Vec<OrfRecord>>,
    needle: &[u8], // b"ATG" or b"M"
    idx_count: &mut usize,
) {
    let mut orf_count = 0;

    for (pos, &aa) in seq.iter().enumerate().skip(1) {
        // INFO: default needle is 'M' -> start signal
        if aa == needle[0] {
            let len_remaining = seq_length - pos;
            let percent = len_remaining as f32 / seq_length as f32;

            // INFO: we have found a nested ORF -> create a new record
            // WARN: updating start codons only in non-ATG cases because we are only looking for inner ATGs
            if len_remaining >= min_len && percent >= min_percent {
                orf_count += 1;
                *idx_count += 1;

                let sub_seq = &seq[pos..];
                let mut nested_header = header.clone();

                nested_header.start += pos * 3;
                nested_header.id += format!("@{}", orf_count).as_str();

                // INFO: only updating start codon if it is not ATG
                if header.start_codon != "ATG" {
                    nested_header.start_codon = "ATG".to_string();
                }

                nested_header.orf_type = match header.orf_type {
                    OrfType::Complete => OrfType::CompleteNested,
                    OrfType::FivePartial => OrfType::FivePartialNested,
                    OrfType::ThreePartial => OrfType::ThreePartialNested,
                    _ => unreachable!(),
                };

                // INFO: subsequence is now a new key in the mapper
                let mut inner_key = Vec::with_capacity(sub_seq.len());
                for &b in sub_seq {
                    if b != b'\n' {
                        inner_key.push(b);
                    }
                }

                let hhead = HashHead::new(inner_key.into());
                match mapper.entry(hhead) {
                    Entry::Occupied(mut entry) => {
                        // INFO: seq exists -> push header to vec
                        // INFO: seq_idx is inherited from first header
                        entry.get_mut().push(nested_header.clone());
                    }
                    Entry::Vacant(entry) => {
                        nested_header.seq_idx = *idx_count;
                        entry.insert(vec![nested_header.clone()]);
                        *idx_count += 1; // INFO: increment index count
                    }
                }
            }
        }
    }
}

/// A HashHead is a tuple of (idx, seq) where idx is the index of the hash
#[derive(Debug, Clone, Eq, Hash, PartialEq)]
pub struct HashHead {
    pub seq: Arc<Vec<u8>>,
}

impl HashHead {
    /// Creates a new `HashHead` from an index and a sequence
    ///
    /// # Arguments
    ///
    /// * `idx` - The index of the hash
    /// * `seq` - The sequence of the hash
    ///
    /// # Returns
    ///
    /// A new `HashHead` instance
    ///
    /// # Example
    ///
    /// ```rust
    /// let hash_head = HashHead::new(0, Arc::from(vec![b'A', b'T', b'G']));
    /// ```
    pub fn new(seq: Arc<Vec<u8>>) -> Self {
        Self { seq }
    }
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
pub fn create_fasta(fasta: &Path, extension: &str) -> Option<BufWriter<File>> {
    let file = fasta.with_extension(extension);
    if !file.exists() {
        let writer = BufWriter::new(
            File::create(&file).unwrap_or_else(|e| panic!("ERROR: cannot create file -> {e}")),
        );

        Some(writer)
    } else {
        warn!("WARN: file already exists -> {:?}!", file.display());

        // INFO: overwrite existing file
        std::fs::remove_file(&file).unwrap_or_else(|e| {
            panic!("ERROR: cannot remove existing file -> {e}");
        });

        Some(BufWriter::new(File::create(&file).unwrap_or_else(|e| {
            panic!("ERROR: cannot create file -> {e}")
        })))
    }
}

/// Creates a new file with the given extension if it doesn't exist
///
/// # Arguments
///
/// * `file` - The path to the file
/// * `extension` - The extension of the file
///
/// # Returns
///
/// An `Option` containing a `BufWriter` for the file
///
/// # Example
///
/// ```rust
/// let writer = create_file("path/to/file.txt", "txt").unwrap();
/// ```
pub fn create_file(file: &Path, extension: &str) -> Option<BufWriter<File>> {
    let file = file.with_extension(extension);
    if !file.exists() {
        let writer = BufWriter::new(
            File::create(&file).unwrap_or_else(|e| panic!("ERROR: cannot create file -> {e}")),
        );

        Some(writer)
    } else {
        // INFO: overwrite existing fil
        dbg!("WARN: overwriting existing file -> {:?}", &file);
        std::fs::remove_file(&file).unwrap_or_else(|e| {
            panic!("ERROR: cannot remove existing file -> {e}");
        });

        Some(BufWriter::new(File::create(&file).unwrap_or_else(|e| {
            panic!("ERROR: cannot create file -> {e}")
        })))
    }
}

/// An enum representing the mode for sequence extraction.
///
/// This enum determines whether sequences are extracted directly ("Raw")
/// or if an indexing approach is used (though `Indexed`)
pub enum Mode {
    Raw,
    Indexed,
}

impl Mode {
    /// Creates an `Mode` from a boolean value.
    ///
    /// # Arguments
    ///
    /// * `mode` - A boolean value. `true` maps to `Mode::Indexed`,
    ///            `false` maps to `Mode::Raw`.
    ///
    /// # Returns
    ///
    /// An `Mode` variant.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// let raw_mode = Mode::from(false);
    /// assert!(matches!(raw_mode, Mode::Raw));
    ///
    /// let indexed_mode = Mode::from(true);
    /// assert!(matches!(indexed_mode, Mode::Indexed));
    /// ```
    pub fn from(mode: &Option<PathBuf>) -> Self {
        match mode {
            Some(_) => Self::Indexed,
            None => Self::Raw,
        }
    }
}

/// Extracts the chromosome name from a file path generated by the pipeline.
///
/// The file path is expected to follow a specific pattern, such as
/// `.../tmp_chunk_chr16:{n}.bed`. This function splits the basename to
/// extract "chr16".
///
/// # Arguments
///
/// * `path` - A `PathBuf` representing the path to the file.
///
/// # Returns
///
/// A `String` containing the chromosome name.
///
/// # Panics
///
/// This function will panic if the file path does not conform to the expected
/// naming convention and its components cannot be extracted successfully.
///
/// # Example
///
/// ```rust
/// use std::path::PathBuf;
///
/// let path1 = PathBuf::from("/path/to/tmp_chunk_chr16:0.bed");
/// let chr1 = get_chr_from_path(&path1);
/// assert_eq!(chr1, "chr16");
///
/// let path2 = PathBuf::from("./outputs/tmp_chunk_chr1:123.index");
/// let chr2 = get_chr_from_path(&path2);
/// assert_eq!(chr2, "chr1");
///
/// // Panicking example
/// // let path3 = PathBuf::from("invalid_name.txt");
/// // let _ = get_chr_from_path(&path3); // This would panic
/// ```
pub fn get_chr_from_path(path: &PathBuf) -> String {
    // INFO: {...}/step8_orf/seqs_free/chr10_KZ289072_fix:0/{filename}
    let basename = path
        .parent()
        .unwrap_or_else(|| panic!("ERROR: could not get basename from {:?}", path))
        .file_stem()
        .unwrap_or_else(|| panic!("ERROR: could not get file stem from {:?}", path))
        .to_string_lossy();

    let chr = basename.split(':').collect::<Vec<&str>>()[0].to_string();

    chr
}

/// Retrieves the absolute CDS (Coding Sequence) coordinates for a given transcript.
///
/// This function looks up the `GenePred` record corresponding to the provided
/// chromosome and ID within the `records` map. It then uses the `map_absolute_cds`
/// method of the `GenePred` struct to convert the alignment start and end coordinates
/// into absolute CDS coordinates.
///
/// # Arguments
///
/// * `records` - A reference to a `DashMap` containing `GenePred` records,
///               keyed by chromosome and then by transcript ID.
/// * `chr` - A string slice representing the chromosome name.
/// * `id` - A string slice representing the transcript ID (e.g., "R1_chr1").
/// * `start` - The start coordinate of the alignment.
/// * `end` - The end coordinate of the alignment.
///
/// # Returns
///
/// A tuple `(u64, u64)` representing the absolute start and end coordinates of the CDS.
/// If the mapping fails or results in default values, it returns `(0, 0)`.
///
/// # Panics
///
/// This function will panic if:
/// - The specified `chr` is not found as a key in the `records` map.
/// - The specified `id` is not found as a key within the chromosome's `HashMap` in `records`.
///
pub fn get_cds_coords(
    records: &DashMap<String, HashMap<String, GenePred>>,
    chr: &str,
    id: &str,
    start: u64,
    end: u64,
) -> (u64, u64) {
    // INFO: retrieving the reference gene prediction record
    let (orf_start, orf_end) = records
        .get_mut(chr)
        .unwrap_or_else(|| {
            panic!(
                "ERROR: chromosome from {} not found in sequences -> {}!",
                id, chr
            );
        })
        .get_mut(id)
        .unwrap_or_else(|| {
            panic!("ERROR: id not found in BED, this is a bug -> {}!", id);
        })
        .map_absolute_cds(start, end);

    (orf_start, orf_end)
}

/// Represents a single BLAST alignment record.
#[derive(Debug, Clone, PartialEq)]
pub struct BlastRecord {
    pub blast_pid: f32,           // Percentage of identical matches
    pub blast_e_value: f64,       // E-value of the match
    pub blast_offset: i32,        // Offset in the query sequence where the match starts
    pub blast_alignment_len: u32, // Length of the alignment
    pub percent_aligned: f32,     // Percentage of the query sequence that is aligned
}

impl BlastRecord {
    /// Creates a new `BlastRecord` from a slice of string parts, typically
    /// obtained by splitting a line from a DIAMOND BLAST output.
    ///
    /// The expected format of `parts` corresponds to `diamond blastp --outfmt 6`
    /// output: `qseqid pident qlen slen length qstart qend sstart send evalue`.
    ///
    /// # Arguments
    ///
    /// * `parts` - A slice of string slices, where each element represents a column from the BLAST output.
    ///
    /// # Returns
    ///
    /// A `BlastRecord` instance populated with the parsed data. The `blast_id`
    /// field is initially empty and is expected to be set later.
    ///
    /// # Panics
    ///
    /// This function will panic if:
    /// - `parts` does not contain at least 10 elements.
    /// - Any of the numeric fields (`blast_idx_id`, `blast_pid`, `blast_e_value`,
    ///   `blast_offset` components, `blast_alignment_len`, query length for `percent_aligned`)
    ///   cannot be successfully parsed into their respective types.
    ///
    /// # Example
    ///
    /// Follows this format:
    ///
    /// qseqid pident  qlen    slen   length qstart    qend   sstart   send     evalue
    ///  17      97.2    142     357     141     1       141     217     357     5.09e-93
    ///
    /// ```rust
    /// let parts = ["1", "99.0", "500", "0", "100", "1", "100", "1", "100", "1e-10"];
    /// let record = BlastRecord::from_parts(&parts);
    /// ```
    pub fn from_parts(parts: &[&str]) -> Self {
        if parts.len() < 10 {
            panic!("ERROR: not enough parts to create BlastRecord -> {parts:?}");
        }

        let blast_pid = parts[1].parse::<f32>().unwrap_or_else(|_| {
            panic!("ERROR: failed to parse blast PID from parts: {:?}", parts);
        });

        let blast_e_value = parts[9].parse::<f64>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast E-value from parts: {:?}",
                parts
            );
        });
        // INFO: if parsed to zero, but string was not "0.0", it's subnormal
        let blast_e_value = if blast_e_value == 0.0 {
            // INFO: represent it with the minimum positive value
            f64::MIN_POSITIVE // INFO: ~2.225074e-308
        } else {
            blast_e_value
        };

        let qstart = parts[5].parse::<i32>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast offset from parts: {:?}",
                parts
            );
        });
        let qend = parts[6].parse::<i32>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast offset from parts: {:?}",
                parts
            );
        });

        let blast_offset = parts[7].parse::<i32>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast offset from parts: {:?}",
                parts
            )
        }) - qstart;

        let blast_alignment_len = parts[4].parse::<u32>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast offset from parts: {:?}",
                parts
            )
        });

        let percent_aligned = (qend - qstart + 1) as f32
            / parts[2].parse::<u32>().unwrap_or_else(|_| {
                panic!(
                    "ERROR: failed to parse blast length from parts: {:?}",
                    parts
                );
            }) as f32
            * 100.0;

        Self {
            blast_pid,
            blast_e_value,
            blast_offset,
            blast_alignment_len,
            percent_aligned,
        }
    }
}

/// Parses a FASTA header string to extract start, end coordinates, and strand information.
///
/// This function uses a provided regular expression to capture specific groups
/// from the header string, typically in the format `[start-end](strand)`.
///
/// # Arguments
///
/// * `header` - A string slice representing the FASTA header.
/// * `capture` - A reference to a compiled `regex::Regex` with capture groups
///               for start, end, and strand.
///
/// # Returns
///
/// An `Option` containing a tuple `(usize, usize, char)` representing
/// (start coordinate, end coordinate, strand character) if parsing is successful.
/// Returns `None` if the header does not match the regex or if parsing of
/// coordinates or strand fails.
pub fn split_header<'a>(header: &'a str, capture: &regex::Regex) -> Option<(usize, usize, char)> {
    let caps = capture.captures(header)?;

    let start = caps[1].parse().ok()?;
    let end = caps[2].parse().ok()?;
    let strand = caps[3].chars().next()?;
    Some((start, end, strand))
}

/// Reads the translationAi output and converts it to a hashmap of `TaiRecord`s
///
/// # Arguments
///
/// * `tai` - The path to the translationAi output
///
/// # Returns
///
/// A `HashMap` of `String` and their corresponding `TaiRecord`
///
/// # Example
///
/// ```rust
/// let tai = read_tai_table(tai).unwrap();
/// ```
pub fn read_tai_table(tai: PathBuf) -> HashMap<String, TaiRecord> {
    let predictions = bed_reader(tai)
        .unwrap_or_else(|e| panic!("ERROR: failed to read blast predictions file -> {e}"));

    let mut accumulator = HashMap::new();

    // INFO: chr17 26054102 26060105 0.p3 - 0.011934959 0.9905367 242 2126 ATG TAG 8
    for line in predictions.lines() {
        let record = TaiRecord::parse(line);
        let cannonical_id =
            record.id.split(".p").next().unwrap_or_else(|| {
                panic!("ERRROR: failed to parse cannonical ID from line: {}", line)
            }); // INFO: ENSMUST00000176751.2.p3 -> ENSMUST00000176751.2
        let key = format!("{}:{}-{}", cannonical_id, record.start, record.end); // INFO: ENSMUST00000176751.2:26054102-26060105
        accumulator.insert(key, record);
    }

    accumulator
}

/// A struct representing a translationAi record
#[derive(Debug, Clone)]
pub struct TaiRecord {
    pub chr: String,
    pub start: usize,
    pub end: usize,
    pub id: String,
    pub strand: char,
    pub start_score: f32,
    pub stop_score: f32,
    pub relative_orf_start: usize,
    pub relative_orf_end: usize,
    pub start_codon: String,
    pub stop_codon: String,
    pub inner_stops: usize,
    pub aa: String,
}

impl TaiRecord {
    pub fn parse(line: &str) -> Self {
        let mut parts = line.split('\t');

        let (
            chr,
            start,
            end,
            id,
            strand,
            start_score,
            stop_score,
            relative_orf_start,
            relative_orf_end,
            start_codon,
            stop_codon,
            inner_stops,
            aa,
        ) = (
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai chr from line: {}", line);
                })
                .to_string(),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai start from line: {}", line);
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse tai start from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai end from line: {}", line);
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!("ERROR: failed to parse tai end from line: {} -> {e}", line);
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai id from line: {}", line);
                })
                .to_string(),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai strand from line: {}", line);
                })
                .chars()
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai strand from line: {}", line);
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai start score from line: {}", line);
                })
                .parse::<f32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse tai start score from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse tai stop score from line: {}", line);
                })
                .parse::<f32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse tai stop score from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to parse relative orf start from line: {}",
                        line
                    );
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse relative orf start from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to parse relative orf end from line: {}",
                        line
                    );
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse relative orf end from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse start codon from line: {}", line);
                })
                .to_string(),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse stop codon from line: {}", line);
                })
                .to_string(),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse inner stops from line: {}", line);
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse inner stops from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse aa from line: {}", line);
                })
                .to_string(),
        );

        Self {
            chr,
            start,
            end,
            id,
            strand,
            start_score,
            stop_score,
            relative_orf_start,
            relative_orf_end,
            start_codon,
            stop_codon,
            inner_stops,
            aa,
        }
    }
}
