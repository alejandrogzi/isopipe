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

use crate::{
    cli::RunMode,
    record::{BinKey, Queue, Record},
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
) -> Result<HashMap<BinKey, Queue>, Box<dyn std::error::Error>> {
    let contents = par_reader(files)?;
    let tracks = par_parse_tracks(&contents)?;

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
fn par_parse_tracks(contents: &str) -> Result<HashMap<BinKey, Queue>, Box<dyn std::error::Error>> {
    let tracks = contents
        .par_lines()
        .filter(|row| !row.starts_with("#"))
        .filter_map(|line| Record::parse(line).ok())
        .fold(HashMap::new, |mut acc, record| {
            // INFO: if record not in tracks, create a new queue
            acc.entry(record.key)
                .and_modify(|queue: &mut Queue| {
                    queue.count += 1;
                    queue.reads.push(record.read.clone());
                })
                .or_insert(Queue {
                    reads: vec![],
                    count: 0,
                    rep_line: record.line,
                    header: record.read,
                });

            acc
        })
        .reduce(HashMap::new, |mut left, right| {
            // INFO: if record not in tracks, create a new queue
            for (key, right_queue) in right {
                left.entry(key)
                    .and_modify(|left_queue| {
                        left_queue.count += right_queue.count + 1; // INFO: accounts for right header

                        left_queue.reads.extend(right_queue.reads.clone());
                        left_queue.reads.push(right_queue.header.clone()); // INFO: accounts for right header
                    })
                    .or_insert(right_queue);
            }

            left
        });

    Ok(tracks)
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
) {
    log::info!("INFO: Writing collapsed file");
    let mut writer = BufWriter::new(
        File::create(&output)
            .unwrap_or_else(|e| panic!("ERROR: Could not create file {:?} -> {e}", output)),
    );

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

                let collapsed = format!("{}\t{}", line, tail);

                writeln!(writer, "{}", collapsed)
                    .unwrap_or_else(|e| panic!("ERROR: Could not write to collapsed file -> {e}"));
            }
            RunMode::Index => {
                // INFO: no extra column
                let line = std::str::from_utf8(&queue.rep_line)
                    .unwrap_or_else(|e| panic!("ERROR: Could not convert rep_line to utf8 -> {e}"));

                writeln!(writer, "{}", line)
                    .unwrap_or_else(|e| panic!("ERROR: Could not write to collapsed file -> {e}"));
            }
        }
    }

    log::info!(
        "INFO: Wrote {} lines to collapsed file",
        writer.buffer().lines().count()
    );
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

    log::info!(
        "INFO: Wrote {} lines to index file",
        writer.buffer().lines().count()
    );
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
