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

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, prelude::*, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::RwLock;
use std::{fmt::Debug, io::BufReader};

use log;
use once_cell::sync::Lazy;
use rayon::prelude::*;

use crate::cli::CollapseMode;
use crate::{
    cli::RunMode,
    record::{BinKey, Queue, QueueState, Record},
};

pub const MAGIC: &[u8; 5] = b"C0IDX";
pub const IDX_VERSION: u8 = 1;

static CHROM_INTERNER: Lazy<RwLock<HashMap<String, u32>>> =
    Lazy::new(|| RwLock::new(HashMap::with_capacity(1024)));
static NEXT_CHROM_ID: AtomicU32 = AtomicU32::new(1);

/// Unpacks and processes multiple BED12 files in parallel.
///
/// This function reads multiple BED12 files concurrently, parses their contents,
/// and groups records by their genomic coordinates into queues for efficient processing.
///
/// # Arguments
///
/// * `files` - A vector of file paths that implement `AsRef<Path> + Debug + Sync + Send`
///
/// # Returns
///
/// A `Result` containing a `HashMap<BinKey, Queue>` where:
/// - Keys are `BinKey` instances representing unique genomic intervals
/// - Values are `Queue` instances containing grouped reads for each interval
///
/// # Errors
///
/// Returns an error if:
/// - Any file cannot be read
/// - File contents cannot be parsed as valid BED12 format
/// - I/O operations fail during file reading
///
/// # Example
///
/// ```rust, ignore
/// let files = vec!["sample1.bed", "sample2.bed", "sample3.bed"];
/// let tracks = unpack(files)?;
/// println!("Processed {} unique genomic intervals", tracks.len());
/// ```
pub fn unpack<P: AsRef<Path> + Debug + Sync + Send>(
    files: Vec<P>,
    mode: CollapseMode,
    merge: bool,
) -> Result<HashMap<BinKey, Queue>, Box<dyn std::error::Error>> {
    let contents = par_reader(files)?;
    let tracks = par_parse_tracks(&contents, mode, merge)?;

    Ok(tracks)
}

/// Parses BED12 content in parallel and groups records by genomic coordinates.
///
/// This function processes multi-line BED12 content using parallel processing,
/// filtering out comments and invalid lines, then grouping valid records by
/// their `BinKey` (genomic coordinates).
///
/// # Arguments
///
/// * `contents` - A string slice containing BED12 formatted data
///
/// # Returns
///
/// A `Result` containing a `HashMap<BinKey, Queue>` with grouped genomic records.
///
/// # Errors
///
/// Returns an error if the parsing operations fail catastrophically.
///
/// # Processing Pipeline
///
/// 1. Split content into parallel lines
/// 2. Filter out comment lines (starting with '#')
/// 3. Parse each line into a `Record` (invalid lines are silently dropped)
/// 4. Group records by `BinKey` using parallel fold/reduce
/// 5. Merge duplicate keys by combining their queues
///
/// # Implementation Notes
///
/// - Uses Rayon's `par_lines()` for parallel line processing
/// - Invalid records are filtered out silently with `filter_map()`
/// - The fold operation creates per-thread `HashMap` instances
/// - The reduce operation merges thread-local maps efficiently
/// - Read counts account for both individual reads and headers
///
/// # Example
///
/// ```rust, ignore
/// let bed_content = "chr1\t1000\t2000\tread1\t100\t+\t1100\t1900\t0\t1\t1000\t0\n";
/// let tracks = par_parse_tracks(bed_content)?;
/// assert!(tracks.len() > 0);
/// ```
fn par_parse_tracks(
    contents: &str,
    mode: CollapseMode,
    merge: bool,
) -> Result<HashMap<BinKey, Queue>, Box<dyn std::error::Error>> {
    let tracks = contents
        .par_lines()
        .filter(|row| !row.starts_with("#"))
        .filter_map(|line| Record::parse(line, &mode).ok())
        .fold(HashMap::new, |mut acc, record| {
            // INFO: if record not in tracks, create a new queue
            acc.entry(record.key)
                .and_modify(|queue: &mut Queue| {
                    queue.count += 1;
                    queue.reads.push(record.read.clone());

                    // INFO: the next block will only happend when collapse_mode is
                    // equal to 'gapped-cds'. The 'exon' and 'transcript' mode will
                    // require an exact BinKey match. The 'cds' mode sets both bounds
                    // to 0 and thus will always fall in the third case.

                    // INFO: update queue bounds, three cases:
                    // 1. if one bound is greater but the other less -> QueueState::Perturbed
                    // and we need to merge
                    //
                    // Example:
                    //
                    // read1: xxxxXXX----XXX--XXXXXXxxxx
                    //        ^^                        ^^^^
                    // read2:   xxXXX----XXX--XXXXXXxxxxxxxx
                    //
                    // or its reverse.
                    //
                    // Here, if merge flag is on, we need to update bounds to min and max,
                    // respectively. After that, we use the updated bounds to create an
                    // artificial read that joins both, replacing the line and updating the
                    // name in that artificial line with the tag ::AR (artificial). Note that
                    // we also need to update queue.reads with the replaced header.
                    if (record.bounds.0 > queue.bounds.0 && record.bounds.1 > queue.bounds.1)
                        || (record.bounds.0 < queue.bounds.0 && record.bounds.1 < queue.bounds.1)
                    {
                        queue.state = QueueState::Perturbed;
                        queue.bounds = (
                            record.bounds.0.min(queue.bounds.0),
                            record.bounds.1.max(queue.bounds.1),
                        );

                        // INFO: not necessary to update header or rep_line, updating reads is enough
                        queue.reads.push(record.read.clone());

                        queue.count += 1;
                        if merge {
                            // INFO: accounting for the current header that will also be part of queue
                            queue.count += 1;
                            queue.reads.push(queue.header.clone());
                        }
                    }
                    // 2. if one of both bounds is greater and the other is either equal or
                    // also greater -> QueueState::Unperturbed, replace current line
                    //
                    // Example:
                    //
                    // read1: xxxxXXXXX----XXXX----XXXXXxxxx
                    //        ^^                           |
                    // read2:   xxXXXXX----XXXX----XXXXXxxxx
                    //
                    // or its reverse
                    else if record.bounds.0 <= queue.bounds.0 && record.bounds.1 >= queue.bounds.1
                    {
                        // INFO: in the case the read was already perturbed, makes
                        // sense to change the state because the new read bounds are
                        // outwards the already updated perturbed bounds
                        queue.state = QueueState::Unperturbed;

                        queue.bounds = record.bounds;
                        queue.rep_line = record.line.clone();
                        queue.count += 1;

                        // INFO: update queue.reads with current read being replaced
                        queue.reads.push(record.read.clone());

                        // INFO: update header
                        queue.header = record.read.clone();
                    }
                    // 3. if both bounds are equal or less -> state remains the same
                    else {
                    }
                })
                .or_insert(Queue {
                    reads: vec![],
                    count: 0,
                    rep_line: record.line,
                    header: record.read,
                    bounds: record.bounds,
                    state: QueueState::Unperturbed,
                });

            acc
        })
        .reduce(HashMap::new, |mut left, right| {
            // INFO: if record not in tracks, create a new queue
            for (key, right_queue) in right {
                left.entry(key)
                    .and_modify(|left_queue| {
                        __compare_queue(left_queue, &right_queue, merge);
                    })
                    .or_insert(right_queue);
            }

            left
        });

    Ok(tracks)
}

/// Compares two `Queue` instances and updates the first one with the second one.
///
/// This function compares two `Queue` instances and updates the first one with the second one.
/// It accounts for the right `Queue` instance's header and updates the left `Queue` instance
/// with perturbed bounds, header, and read information.
///
/// # Arguments
///
/// * `left` - A mutable reference to the left `Queue` instance
/// * `right` - A reference to the right `Queue` instance
/// * `merge` - A boolean flag indicating whether to merge the right `Queue` instance
///
/// # Example
///
/// ```rust, ignore
/// let mut left = Queue::default();
/// let right = Queue::default();
/// __compare_queue(&mut left, &right, false);
/// ```
///
/// See [`par_parse_tracks`] for descriptions of each branch
fn __compare_queue(left: &mut Queue, right: &Queue, merge: bool) {
    // INFO: +1 accounts either for left or right
    left.count += right.count + 1;

    // INFO: accounting for right header
    left.reads.extend(right.reads.clone());

    if (right.bounds.0 > left.bounds.0 && right.bounds.1 > left.bounds.1)
        || (right.bounds.0 < left.bounds.0 && right.bounds.1 < left.bounds.1)
    {
        left.state = QueueState::Perturbed;
        left.bounds = (
            right.bounds.0.min(left.bounds.0),
            right.bounds.1.max(left.bounds.1),
        );
        left.reads.push(left.header.clone()); // INFO: sending left header to queue
        left.header = right.header.clone();

        if merge {
            // INFO: accounting for the current header that will also be part of queue
            left.count += 1;
            left.reads.push(left.header.clone()); // INFO: sending right header to queue
        }
    } else if right.bounds.0 <= left.bounds.0 && right.bounds.1 >= left.bounds.1 {
        // INFO: in the case the read was already perturbed, makes
        // sense to change the state because the new read bounds are
        // outwards the already updated perturbed bounds
        left.state = QueueState::Unperturbed;

        // WARN: not adding +1 because we already do it at the beginning
        left.bounds = right.bounds;
        left.rep_line = right.rep_line.clone();

        // INFO: update queue.reads with current read being replaced
        left.reads.push(left.header.clone());
        left.header = right.header.clone();
    } else {
        left.reads.push(right.header.clone()); // INFO: accounts for right header
    }
}

/// Reads the entire contents of a file into a string.
///
/// # Arguments
///
/// * `file` - A file path that implements `AsRef<Path> + Debug`
///
/// # Returns
///
/// A `Result` containing the file contents as a `String`, or an error if reading fails.
///
/// # Errors
///
/// Returns an error if:
/// - The file cannot be opened (doesn't exist, permission denied, etc.)
/// - I/O errors occur during reading
/// - The file contains invalid UTF-8 sequences
///
/// # Example
///
/// ```rust, ignore
/// let content = reader("data.bed")?;
/// println!("File contains {} characters", content.len());
/// ```
fn reader<P: AsRef<Path> + Debug>(file: P) -> Result<String, Box<dyn std::error::Error>> {
    let mut file = File::open(file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}

/// Reads multiple files in parallel and concatenates their contents.
///
/// This function leverages Rayon's parallel iterator to read multiple files
/// concurrently, then concatenates all file contents into a single string.
///
/// # Arguments
///
/// * `files` - A vector of file paths that implement `AsRef<Path> + Debug + Sync + Send`
///
/// # Returns
///
/// A `Result` containing the concatenated contents of all files as a single `String`.
///
/// # Errors
///
/// Returns an error if any file cannot be read successfully.
///
/// # Panics
///
/// Panics if any individual file reading operation fails, with a descriptive error message
/// indicating which file caused the failure.
///
/// # Performance Notes
///
/// - Files are read concurrently using Rayon's parallel iterators
/// - Memory usage scales with the total size of all input files
/// - No attempt is made to add separators between file contents
///
/// # Example
///
/// ```rust, ignore
/// let files = vec!["file1.bed", "file2.bed", "file3.bed"];
/// let combined_content = par_reader(files)?;
/// println!("Combined {} files into {} characters", files.len(), combined_content.len());
/// ```
fn par_reader<P: AsRef<Path> + Debug + Sync + Send>(
    files: Vec<P>,
) -> Result<String, Box<dyn std::error::Error>> {
    let contents: Vec<String> = files
        .par_iter()
        .map(|path| reader(path).unwrap_or_else(|e| panic!("Error reading file: {:?}", e)))
        .collect();

    Ok(contents.concat())
}

/// Interns a chromosome name to a small u32 identifier in a thread-safe manner.
///
/// This function provides a thread-safe string interning system for chromosome names,
/// converting them to compact u32 identifiers for memory efficiency. Uses a two-phase
/// locking strategy for optimal read performance.
///
/// # Arguments
///
/// * `chrom` - A string slice containing the chromosome name (e.g., "chr1", "chrX", "scaffold_1")
///
/// # Returns
///
/// A `u32` identifier that uniquely represents the chromosome name.
///
/// # Thread Safety
///
/// This function is fully thread-safe and uses the following locking strategy:
/// 1. **Fast path**: Attempts to read with a shared read lock
/// 2. **Slow path**: Acquires exclusive write lock only if chromosome is not found
/// 3. **Double-checked locking**: Rechecks existence after acquiring write lock
///
/// # Panics
///
/// Panics if the internal `RwLock` is poisoned due to a panic in another thread.
///
/// # Example
///
/// ```rust, ignore
/// let chr1_id = intern_chromosome("chr1");
/// let chr2_id = intern_chromosome("chr2");
/// let chr1_again = intern_chromosome("chr1");
///
/// assert_eq!(chr1_id, chr1_again); // Same chromosome returns same ID
/// assert_ne!(chr1_id, chr2_id);    // Different chromosomes get different IDs
/// ```
///
/// # Global State
///
/// This function modifies global static variables:
/// - `CHROM_INTERNER`: Thread-safe map of chromosome names to IDs
/// - `NEXT_CHROM_ID`: Atomic counter for assigning new IDs
pub fn intern_chromosome(chrom: &str) -> u32 {
    // INFO: Fast path (read lock)
    if let Ok(read) = CHROM_INTERNER.read() {
        if let Some(&id) = read.get(chrom) {
            return id;
        }
    }

    // INFO: Slow path (write lock)
    let mut write = CHROM_INTERNER
        .write()
        .unwrap_or_else(|e| panic!("ERROR: poisoned lock: {e} -> {chrom} - {CHROM_INTERNER:?}"));

    if let Some(&id) = write.get(chrom) {
        return id;
    }

    let id = NEXT_CHROM_ID.fetch_add(1, Ordering::Relaxed);
    write.insert(chrom.to_string(), id);
    id
}

/// Writes collapsed genomic tracks to a file in different output modes.
///
/// This function outputs the processed genomic tracks to a file, with different
/// formatting options based on the specified run mode.
///
/// # Arguments
///
/// * `tracks` - A `HashMap<BinKey, Queue>` containing the processed genomic data
/// * `mode` - A `RunMode` enum specifying the output format:
///   - `RunMode::Extend`: Includes read names as an additional tab-separated column
///   - `RunMode::Index`: Outputs only the representative line without extra columns
/// * `output` - A file path that implements `AsRef<Path> + Debug` where output will be written
///
/// # Output Formats
///
/// ## RunMode::Extend
/// ```text
/// chr1\t1000\t2000\titem1\t100\t+\t1100\t1900\t0\t1\t1000\t0\tread1,read2,read3
/// ```
///
/// ## RunMode::Index
/// ```text
/// chr1\t1000\t2000\titem1\t100\t+\t1100\t1900\t0\t1\t1000\t0
/// ```
///
/// # Panics
///
/// Panics if:
/// - The output file cannot be created
/// - Write operations fail
/// - Byte data cannot be converted to valid UTF-8
///
/// # Example
///
/// ```rust, ignore
/// let tracks = unpack(vec!["input.bed"])?;
/// __write_collapsed(tracks, RunMode::Extend, "output.bed");
/// ```
pub fn __write_collapsed<P: AsRef<Path> + Debug>(
    tracks: HashMap<BinKey, Queue>,
    mode: RunMode,
    output: P,
    merge: bool,
) {
    log::info!("INFO: Writing collapsed file");
    let mut writer = BufWriter::new(
        File::create(&output)
            .unwrap_or_else(|e| panic!("ERROR: Could not create file {:?} -> {e}", output)),
    );

    let lines = tracks.len();
    let mut artificial = 0;

    if merge {
        log::warn!("WARN: Merge mode activated, will create artificial reads!");
    }

    for (_, queue) in tracks {
        match mode {
            RunMode::Extend => {
                // INFO: appending queue.reads as extra column
                let tail = queue
                    .reads
                    .iter()
                    .map(|r| {
                        std::str::from_utf8(r).unwrap_or_else(|e| {
                            panic!("ERROR: Could not convert read to utf8 -> {e}")
                        })
                    })
                    .collect::<Vec<_>>()
                    .join(",");

                let line = std::str::from_utf8(&queue.rep_line)
                    .unwrap_or_else(|e| panic!("ERROR: Could not convert rep_line to utf8 -> {e}"));

                match queue.state {
                    QueueState::Unperturbed => {
                        let collapsed = format!("{}\t{}", line, tail);
                        __write_line(&mut writer, &collapsed);
                    }
                    QueueState::Perturbed => {
                        if merge {
                            log::debug!(
                                "DEBUG: Perturbed state, writing artificial read -> {queue:?}"
                            );

                            artificial += 1;
                            let mut parts = line.split("\t").collect::<Vec<_>>();

                            let old_start = parts[1]
                                .parse::<u32>()
                                .unwrap_or_else(|e| panic!("ERROR: Could not parse start -> {e}"));
                            let old_end = parts[2]
                                .parse::<u32>()
                                .unwrap_or_else(|e| panic!("ERROR: Could not parse end -> {e}"));

                            let new_start = format!("{}", queue.bounds.0);
                            let new_end = format!("{}", queue.bounds.1);
                            let new_name = format!("{}#AR", parts[3]);

                            parts[1] = new_start.as_str();
                            parts[2] = new_end.as_str();
                            parts[3] = new_name.as_str();

                            let (new_sizes, new_starts) = update_blocks(
                                &mut parts,
                                old_start,
                                old_end,
                                queue.bounds.0,
                                queue.bounds.1,
                            );

                            parts[10] = new_sizes.as_str();
                            parts[11] = new_starts.as_str();

                            let collapsed = parts.join("\t");
                            __write_line(&mut writer, &collapsed);
                        } else {
                            // INFO: just write the line, no merge
                            let collapsed = format!("{}\t{}", line, tail);
                            __write_line(&mut writer, &collapsed);
                        }
                    }
                }
            }
            RunMode::Index => {
                // INFO: no extra column
                let line = std::str::from_utf8(&queue.rep_line)
                    .unwrap_or_else(|e| panic!("ERROR: Could not convert rep_line to utf8 -> {e}"));

                match queue.state {
                    QueueState::Unperturbed => {
                        __write_line(&mut writer, line);
                    }
                    QueueState::Perturbed => {
                        if merge {
                            log::debug!(
                                "DEBUG: Perturbed state, writing artificial read -> {queue:?}"
                            );

                            artificial += 1;
                            let mut parts = line.split("\t").collect::<Vec<_>>();

                            let old_start = parts[1]
                                .parse::<u32>()
                                .unwrap_or_else(|e| panic!("ERROR: Could not parse start -> {e}"));
                            let old_end = parts[2]
                                .parse::<u32>()
                                .unwrap_or_else(|e| panic!("ERROR: Could not parse end -> {e}"));

                            let new_start = format!("{}", queue.bounds.0);
                            let new_end = format!("{}", queue.bounds.1);
                            let new_name = format!("{}#AR", parts[3]);

                            parts[1] = new_start.as_str();
                            parts[2] = new_end.as_str();
                            parts[3] = new_name.as_str();

                            let (new_sizes, new_starts) = update_blocks(
                                &mut parts,
                                old_start,
                                old_end,
                                queue.bounds.0,
                                queue.bounds.1,
                            );

                            parts[10] = new_sizes.as_str();
                            parts[11] = new_starts.as_str();

                            let collapsed = parts.join("\t");
                            __write_line(&mut writer, &collapsed);
                        } else {
                            // INFO: just write the line, no merge
                            __write_line(&mut writer, line);
                        }
                    }
                }
            }
        }
    }

    log::info!("INFO: Wrote {} lines to collapsed file", lines);

    if merge {
        log::warn!(
            "WARN: Merge mode activated, wrote {} artificial reads",
            artificial
        );
    }
}

/// Writes a line to a file.
///
/// This function writes a line to a file using a buffered writer.
///
/// # Arguments
///
/// * `writer` - A reference to a `BufWriter` instance
/// * `line` - A string slice containing the line to write
///
/// # Panics
///
/// Panics if there is an error writing to the file.
fn __write_line(writer: &mut BufWriter<File>, line: &str) {
    writeln!(writer, "{}", line)
        .unwrap_or_else(|e| panic!("ERROR: Could not write to collapsed file -> {e}"));
}

/// Updates block sizes and starts for a BED12 line.
///
/// This function updates the block sizes and starts for a BED12 line
/// based on the new bounds and the old bounds. It accounts for the
/// following cases:
///
/// 1. Extended at start: The first block size is increased by the difference
///    between the new start and the old start.
/// 2. Trimmed at start: The first block size is decreased by the difference
///    between the new start and the old start.
/// 3. Extended at end: The last block size is increased by the difference
///    between the new end and the old end.
/// 4. Trimmed at end: The last block size is decreased by the difference
///    between the new end and the old end.
/// 5. Adjusted exon starts: The exon starts are adjusted based on the
///    difference between the new start and the old start.
///
/// # Arguments
///
/// * `parts` - A mutable reference to a vector of strings representing the BED12 line
/// * `old_start` - The old start position
/// * `old_end` - The old end position
/// * `new_start` - The new start position
/// * `new_end` - The new end position
///
/// # Returns
///
/// A tuple containing the updated block sizes and starts as strings
///
/// # Example
///
/// ```rust, ignore
/// let mut parts = line.split("\t").collect::<Vec<_>>();
/// let old_start = parts[1]
///     .parse::<u32>()
///     .unwrap_or_else(|e| panic!("ERROR: Could not parse start -> {e}"));
/// let old_end = parts[2]
///     .parse::<u32>()
///     .unwrap_or_else(|e| panic!("ERROR: Could not parse end -> {e}"));
///
/// let new_start = format!("{}", queue.bounds.0);
/// let new_end = format!("{}", queue.bounds.1);
/// let new_name = format!("{}#AR", parts[3]);
///
/// let (new_sizes, new_starts) = update_blocks(
///     &mut parts,
///     old_start,
///     old_end,
///     queue.bounds.0,
///     queue.bounds.1,
/// );
///
/// parts[10] = new_sizes.as_str();
/// parts[11] = new_starts.as_str();
/// ```
fn update_blocks<'a>(
    parts: &mut Vec<&'a str>,
    old_start: u32,
    old_end: u32,
    new_start: u32,
    new_end: u32,
) -> (String, String) {
    // Parse block sizes and starts
    let block_sizes: Vec<u32> = parts[10]
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect();
    let block_starts: Vec<u32> = parts[11]
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect();

    let mut new_block_sizes = block_sizes.clone();
    let mut new_block_starts = block_starts.clone();

    let start_diff = (new_start as i64 - old_start as i64).abs() as u32;

    if new_start < old_start {
        // INFO: Extended at start
        new_block_sizes[0] += start_diff;
        new_block_starts[0] = 0;
    } else {
        // INFO: Trimmed at start
        new_block_sizes[0] -= start_diff;
        new_block_starts[0] = 0;
    }

    // INFO: Adjust all exon starts after the first one
    let first_exon_delta = new_start as i64 - old_start as i64;
    for i in 1..new_block_starts.len() {
        new_block_starts[i] = (block_starts[i] as i64 - first_exon_delta) as u32;
    }

    // INFO: Adjust last exon
    let last_idx = new_block_sizes.len() - 1;
    let end_diff = (new_end as i64 - old_end as i64).abs() as u32;
    if new_end > old_end {
        // INFO: Extended at end
        new_block_sizes[last_idx] += end_diff;
    } else {
        // INFO: Trimmed at end
        new_block_sizes[last_idx] -= end_diff;
    }

    let new_sizes_str = new_block_sizes
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let new_starts_str = new_block_starts
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>()
        .join(",");

    (new_sizes_str, new_starts_str)
}

/// Writes a binary index file containing genomic tracks data.
///
/// This function creates a custom binary index format that stores all genomic
/// track information in a compact, efficiently readable format suitable for
/// fast random access and reconstruction.
///
/// # Arguments
///
/// * `tracks` - A reference to `HashMap<BinKey, Queue>` containing the genomic data
/// * `output` - A file path that implements `AsRef<Path> + Debug` for the output index file
///
/// # Binary Format Structure
///
/// The index file uses the following binary layout:
/// ```text
/// Header:
///   - Magic bytes (5 bytes): File format identifier
///   - Version (1 byte): Index format version
///   - Number of records (8 bytes, little-endian u64)
///
/// For each record:
///   - Key length (4 bytes, little-endian u32)
///   - Key data (variable length): Canonical BinKey bytes
///   - Header length (4 bytes, little-endian u32)
///   - Header data (variable length bytes)
///   - Count (4 bytes, little-endian u32)
///   - Rep line length (4 bytes, little-endian u32)
///   - Rep line data (variable length bytes)
///   - Number of reads (4 bytes, little-endian u32)
///   - For each read:
///     - Read length (4 bytes, little-endian u32)
///     - Read data (variable length bytes)
/// ```
///
/// # Example
///
/// ```rust, ignore
/// let tracks = unpack(vec!["input.bed"])?;
/// __write_index(&tracks, "genome.idx");
///
/// // Later, read the index back
/// let reconstructed = read_index("genome.idx")?;
/// assert_eq!(tracks.len(), reconstructed.len());
/// ```
///
/// # See Also
///
/// - [`read_index`] for reading index files created by this function
/// - [`BinKey::to_bytes_canonical`] for the key serialization format
pub fn __write_index<P: AsRef<Path> + Debug>(tracks: &HashMap<BinKey, Queue>, output: P) {
    log::info!("INFO: Writing index file");

    let mut writer = BufWriter::new(
        File::create(&output)
            .unwrap_or_else(|e| panic!("ERROR: Could not create file {:?} -> {e}", output)),
    );

    writer.write_all(MAGIC).unwrap_or_else(|e| {
        panic!(
            "ERROR: Could not write magic to index file {:?} -> {e}",
            output
        )
    });

    writer.write_all(&[IDX_VERSION]).unwrap_or_else(|e| {
        panic!(
            "ERROR: Could not write version to index file {:?} -> {e}",
            output
        )
    });

    let num_records = tracks.len() as u64;
    writer
        .write_all(&num_records.to_le_bytes())
        .unwrap_or_else(|e| {
            panic!(
                "ERROR: Could not write number of records to index file {:?} -> {e}",
                output
            )
        });

    for (key, queue) in tracks {
        // INFO: write key
        let key_bytes = key.to_bytes_canonical();
        let key_len = key_bytes.len() as u32;
        writer
            .write_all(&key_len.to_le_bytes())
            .unwrap_or_else(|e| {
                panic!(
                    "ERROR: Could not write key length to index file {:?} -> {e}",
                    output
                )
            });
        writer.write_all(&key_bytes).unwrap_or_else(|e| {
            panic!(
                "ERROR: Could not write key to index file {:?} -> {e}",
                output
            )
        });

        // INFO: write header
        let header_len = queue.header.len() as u32;
        writer
            .write_all(&header_len.to_le_bytes())
            .unwrap_or_else(|e| {
                panic!(
                    "ERROR: Could not write header length to index file {:?} -> {e}",
                    output
                )
            });
        writer.write_all(&queue.header).unwrap_or_else(|e| {
            panic!(
                "ERROR: Could not write header to index file {:?} -> {e}",
                output
            )
        });

        // INFO: write count
        writer
            .write_all(&queue.count.to_le_bytes())
            .unwrap_or_else(|e| {
                panic!(
                    "ERROR: Could not write count to index file {:?} -> {e}",
                    output
                )
            });

        // INFO: write rep_line
        let rep_line_len = queue.rep_line.len() as u32;
        writer
            .write_all(&rep_line_len.to_le_bytes())
            .unwrap_or_else(|e| {
                panic!(
                    "ERROR: Could not write rep_line length to index file {:?} -> {e}",
                    output
                )
            });
        writer.write_all(&queue.rep_line).unwrap_or_else(|e| {
            panic!(
                "ERROR: Could not write rep_line to index file {:?} -> {e}",
                output
            )
        });

        // INFO: write reads
        let n_reads = queue.reads.len() as u32;
        writer
            .write_all(&n_reads.to_le_bytes())
            .unwrap_or_else(|e| {
                panic!(
                    "ERROR: Could not write number of reads to index file {:?} -> {e}",
                    output
                )
            });
        for read in &queue.reads {
            let rl = read.len() as u32;
            writer.write_all(&rl.to_le_bytes()).unwrap_or_else(|e| {
                panic!(
                    "ERROR: Could not write read length to index file {:?} -> {e}",
                    output
                )
            });
            writer.write_all(read).unwrap_or_else(|e| {
                panic!(
                    "ERROR: Could not write read to index file {:?} -> {e}",
                    output
                )
            });
        }
    }

    log::info!("INFO: Wrote {} lines to index file", num_records);
}

/// Reads and reconstructs genomic tracks data from a binary index file.
///
/// This function reads a binary index file created by [`__write_index`] and
/// reconstructs the original `HashMap<BinKey, Queue>` data structure with
/// full fidelity to the original data.
///
/// # Arguments
///
/// * `input` - A file path that implements `AsRef<Path>` pointing to the index file
///
/// # Returns
///
/// An `io::Result<HashMap<BinKey, Queue>>` containing:
/// - `Ok(HashMap)`: Successfully reconstructed genomic tracks data
/// - `Err(io::Error)`: File I/O errors or format validation failures
///
/// # Errors
///
/// Returns `io::Error` for the following conditions:
///
/// ## File Access Errors
/// - File does not exist or cannot be opened
/// - Insufficient permissions to read the file
/// - I/O errors during reading operations
///
/// ## Format Validation Errors
/// - **Invalid Magic Bytes**: File is not a valid index format
/// - **Unsupported Version**: Index was created with an incompatible version
/// - **Truncated Data**: File ends unexpectedly during reading
/// - **Corrupted Length Fields**: Invalid length values that exceed remaining data
///
/// # Binary Format Compatibility
///
/// The function expects the exact binary format produced by [`__write_index`]:
/// - Magic bytes must match `MAGIC` constant exactly
/// - Version byte must match `IDX_VERSION` exactly
/// - All multi-byte values must be in little-endian format
/// - Length-prefixed data must be consistent
///
/// # Memory Usage
///
/// - Pre-allocates HashMap with capacity based on record count
/// - Memory usage scales with the total size of all genomic data
/// - Each `BinKey` is reconstructed from canonical bytes with fingerprint computation
///
/// # Example
///
/// ```rust, ignore
/// // Write an index file
/// let original_tracks = unpack(vec!["input.bed"])?;
/// __write_index(&original_tracks, "genome.idx");
///
/// // Read it back
/// match read_index("genome.idx") {
///     Ok(reconstructed_tracks) => {
///         assert_eq!(original_tracks.len(), reconstructed_tracks.len());
///         println!("Successfully loaded {} genomic intervals", reconstructed_tracks.len());
///     }
///     Err(e) => eprintln!("Failed to read index: {}", e),
/// }
/// ```
///
/// # Data Integrity
///
/// - BinKey fingerprints are recomputed during reconstruction
/// - All byte data is preserved exactly as originally stored
/// - HashMap iteration order may differ from original due to hash recomputation
///
/// # See Also
///
/// - [`__write_index`] for creating compatible index files
/// - [`BinKey::from_bytes_canonical`] for key reconstruction details
pub fn read_index<P: AsRef<Path>>(input: P) -> io::Result<HashMap<BinKey, Queue>> {
    let file = File::open(input)?;
    let mut reader = BufReader::new(file);

    // INFO: read magic and version
    let mut magic = [0u8; 5];
    reader.read_exact(&mut magic)?;
    if &magic != MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid index magic: expected {:?}, got {:?}", MAGIC, magic),
        ));
    }
    let mut version_buf = [0u8; 1];
    reader.read_exact(&mut version_buf)?;
    let version = version_buf[0];
    if version != IDX_VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Unsupported index version: {}", version),
        ));
    }

    // INFO: number of records
    let mut num_records_buf = [0u8; 8];
    reader.read_exact(&mut num_records_buf)?;
    let num_records = u64::from_le_bytes(num_records_buf);

    let mut map = HashMap::with_capacity(num_records as usize);

    for _ in 0..num_records {
        // INFO: key
        let mut u32_buf = [0u8; 4];

        reader.read_exact(&mut u32_buf)?;
        let key_len = u32::from_le_bytes(u32_buf) as usize;
        let mut key_bytes = vec![0u8; key_len];
        reader.read_exact(&mut key_bytes)?;

        // INFO: header
        reader.read_exact(&mut u32_buf)?;
        let header_len = u32::from_le_bytes(u32_buf) as usize;
        let mut header = vec![0u8; header_len];
        reader.read_exact(&mut header)?;

        // INFO: count
        reader.read_exact(&mut u32_buf)?;
        let count = u32::from_le_bytes(u32_buf);

        // INFO: rep_line
        reader.read_exact(&mut u32_buf)?;
        let rep_line_len = u32::from_le_bytes(u32_buf) as usize;
        let mut rep_line = vec![0u8; rep_line_len];
        reader.read_exact(&mut rep_line)?;

        // INFO: reads
        reader.read_exact(&mut u32_buf)?;
        let n_reads = u32::from_le_bytes(u32_buf) as usize;
        let mut reads: Vec<Vec<u8>> = Vec::with_capacity(n_reads);
        for _ in 0..n_reads {
            reader.read_exact(&mut u32_buf)?;
            let rl = u32::from_le_bytes(u32_buf) as usize;
            let mut r = vec![0u8; rl];
            reader.read_exact(&mut r)?;
            reads.push(r);
        }

        // INFO: reconstruct BinKey from canonical bytes
        let key = BinKey::from_bytes_canonical(&key_bytes);

        let queue = Queue {
            header,
            reads,
            count,
            rep_line,
            state: QueueState::Unperturbed,
            bounds: (0, 0), // INFO: placeholder bounds for now
        };

        map.insert(key, queue);
    }

    Ok(map)
}

/// Reads and swaps the index file for a `HashMap<Vec<u8>, Queue>`.
///
/// This function reads the index file and swaps the keys and values
/// in the `HashMap` to maintain the original order of the genomic intervals.
///
/// # Arguments
///
/// * `index` - A file path that implements `AsRef<Path>` pointing to the index file
///
/// # Returns
///
/// A `HashMap<Vec<u8>, Queue>` containing the swapped index data.
///
/// # Errors
///
/// Returns an error if the index file cannot be read.
///
/// # Example
///
/// ```rust, ignore
/// let original_tracks = unpack(vec!["input.bed"])?;
/// __write_index(&original_tracks, "genome.idx");
///
/// let reconstructed_tracks = read_and_swap_index("genome.idx");
/// assert_eq!(original_tracks.len(), reconstructed_tracks.len());
/// ```
pub fn read_and_swap_index<P: AsRef<Path> + Debug>(index: P) -> HashMap<Vec<u8>, Queue> {
    log::info!("INFO: Reading index file from {:?}", index.as_ref());
    let mut map = read_index(&index).unwrap_or_else(|e| {
        panic!(
            "ERROR: Could not read index file {:?} -> {e}",
            index.as_ref()
        )
    });

    log::info!("INFO: Swapping index file");
    let mut swapped = HashMap::new();
    for (_, queue) in map.drain() {
        swapped.insert(queue.header.clone(), queue);
    }

    log::info!("INFO: Swapped {} lines", swapped.len());
    swapped
}
