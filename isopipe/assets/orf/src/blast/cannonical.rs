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
use log::{error, warn};
use packbed::record::GenePred;

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::str::from_utf8;
use std::sync::Arc;

use crate::utils::*;

/// Processes and inflates "canonical" Blast prediction records, and writes them to a file.
///
/// This function is used when the input sequences were *not* indexed in the `extract` step
/// (i.e., `ExtractMode::Raw`). It takes a map of Blast predictions (where keys are the
/// original sequence IDs) and a binary index file that maps these original IDs to their
/// full genomic context.
///
/// For each Blast prediction:
/// 1. It retrieves the corresponding original read information (read ID, chromosome, ORF,
///    subsequence ORF, start, end) from the index.
/// 2. It constructs a full `query_id` string combining these details.
/// 3. It calculates the absolute CDS (Coding Sequence) coordinates using the `GenePred` records.
/// 4. If the CDS coordinates are valid (non-zero), it writes the inflated Blast record
///    to the output file.
///
/// After processing all predictions with hits, it iterates through any remaining IDs in the
/// index (those that did not have a Blast hit) and writes "dummy" records for them to the
/// output file, tagged with "#DM" (Denoted Missing), if their CDS coordinates are valid.
///
/// # Arguments
///
/// * `index` - A `PathBuf` to the binary index file for the chunk.
/// * `predictions` - A `DashMap` where keys are the original sequence IDs (from the raw FASTA)
///                   and values are `BlastRecord`s.
/// * `records` - A `DashMap` containing the original `GenePred` records, keyed
///               by chromosome and then by read ID. Used to get CDS coordinates.
/// * `writer` - A `BufWriter<File>` to which the processed and inflated records are written.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to read the index file.
/// - A sequential ID from `predictions` is not found in the `index` map.
/// - It fails to convert chromosome bytes to a UTF-8 string.
/// - It fails to write to the output `writer`.
/// - `get_cds_coords` or any parsing of coordinates fails.
pub fn cannonical(
    index: &PathBuf,
    predictions: DashMap<String, Arc<BlastRecord>>,
    records: DashMap<String, HashMap<String, GenePred>>,
    mut writer: BufWriter<File>,
) {
    let mut index = read_index(index);

    // INFO: inflate results!
    predictions.iter_mut().for_each(|mut record| {
        let (r_id, data) = record.pair_mut();

        let parts = r_id.split('_').collect::<Vec<&str>>();
        let id = parts[0]
            .parse::<u32>()
            .unwrap_or_else(|e| panic!("ERROR: {e}; could not parse id to u32 -> {r_id:?}"));

        // INFO: unpacking index reference -> queries
        // INFO: for each query all blast records
        // INFO: { index_id : [(read_id: u16, chr_bytes: [u8; chr_len], orf: u16, subseq_orf: u16, seq_len: usize)] }
        // INFO: { 0 : [(read_id: 5903, chr_bytes: [16, 32], orf: 1, subseq_orf: 3, seq_len: 350)] }
        let queries = index.get(&id).unwrap_or_else(|| {
            panic!("ERROR: no queries found for ID: {}", id);
        });

        for query in queries.into_iter() {
            let (read_id, chr, orf, subseq_orf, start, end) = query;

            let chr = format!("chr{}", from_utf8(&chr).unwrap());
            let cannonical_id = format!("R{}_{}", read_id, chr);
            let query_id = format!("{}.p{}@{}", cannonical_id, orf, subseq_orf);

            let (orf_start, orf_end) =
                get_cds_coords(&records, &chr, &cannonical_id, *start as u64, *end as u64);

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

            writer
                .write_all(
                    format!(
                        "{}\t{}\t{:.2}\t{:e}\t{}\t{}\t{:2}\t{}\t{}\t{}\t{}\n",
                        query_id,
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
        index.remove(&id);
    });

    // INFO: repeating the process for unused ids
    // INFO: add tag DM to the ID -> identify unused ids
    index.iter().for_each(|(id, queries)| {
        for query in queries {
            let (read_id, chr, orf, subseq_orf, start, end) = query;

            let chr = format!("chr{}", from_utf8(&chr).unwrap());
            let cannonical_id = format!("R{}_{}", read_id, chr);
            let query_id = format!("{}.p{}@{}#DM", cannonical_id, orf, subseq_orf);

            // INFO: retrieving the reference gene prediction record
            let (orf_start, orf_end) =
                get_cds_coords(&records, &chr, &cannonical_id, *start as u64, *end as u64);

            // WARN: skipping unreliable ORFs for the current alignment
            // INFO: none of these will match any other prediction because
            // INFO: the fall off any exonic boundary
            if orf_start == 0 && orf_end == 0 {
                warn!(
                    "WARN: ORF start and end are zero for ID: {}, skipping!",
                    query_id
                );
                continue;
            }

            writer
                .write_all(
                    format!(
                        "{}\t{}\t{:.2}\t{:e}\t{}\t{}\t{:2}\t{}\t{}\t{}\t{}\n",
                        query_id, id, 0.0, 1.0, 0, 0, 0.0, start, end, orf_start, orf_end
                    )
                    .as_bytes(),
                )
                .unwrap_or_else(|e| {
                    panic!("ERROR: failed to write unused blast record to file -> {e}");
                });
        }
    });
}

/// Creates an index file and a deduplicated FASTA file from a `HashMap` of sequences.
///
/// This function iterates through the provided `mapper` (which contains unique sequences
/// and their associated original record headers). For each unique sequence, it writes
/// the sequence to the deduplicated FASTA file and creates a corresponding entry
/// in the index file. The index entry contains metadata about the original records
/// that map to this unique sequence, including read ID, chromosome, ORF information,
/// and genomic coordinates.
///
/// The index file format is as follows:
/// - `group_id`: `u32` (the unique ID assigned to the deduplicated sequence)
/// - `n_headers`: `u32` (number of original records mapping to this sequence)
/// - For the first original record in the group:
///     - `chr_len`: `u8` (length of the chromosome name in bytes)
///     - `chr_bytes`: `[u8; chr_len]` (chromosome name as bytes)
/// - For each original record in the group:
///     - `read_id`: `u16`
///     - `orf`: `u16`
///     - `subseq_orf`: `u16`
///     - `start`: `u32`
///     - `end`: `u32`
///
/// # Arguments
///
/// * `mapper` - A `HashMap` where keys are unique sequences (`Vec<u8>`) and values
///              are vectors of `Arc<[u8]>` representing the original headers
///              that correspond to that unique sequence.
/// * `index_writer` - A mutable `BufWriter<File>` for writing the index data.
/// * `dedup_writer` - A mutable `BufWriter<File>` for writing the deduplicated FASTA sequences.
///
/// # Panics
///
/// This function will panic if:
/// - There are no sequences in the `mapper`.
/// - It fails to write to the `index_writer` or `dedup_writer`.
/// - It encounters a header that cannot be parsed by `get_read_encoding`.
/// - It fails to convert a byte slice back to a UTF-8 string for writing to FASTA.
///
/// # Example
///
/// ```rust, no_run
/// use std::collections::HashMap;
/// use std::io::{BufWriter, Cursor};
/// use std::fs::File;
/// use std::sync::Arc;
///
/// let mut mapper: HashMap<Vec<u8>, Vec<Arc<[u8]>>> = HashMap::new();
/// mapper.insert(b"ATGC".to_vec(), vec![Arc::from(b"R1_chr1_ORF.1 [1-4](+)_@0".to_vec())]);
///
/// // Create dummy writers for the example
/// let mut index_file = BufWriter::new(File::create("dummy.index").unwrap());
/// let mut dedup_file = BufWriter::new(File::create("dummy.dedup.fa").unwrap());
///
/// make_index(mapper, &mut index_file, &mut dedup_file);
/// ```
pub fn make_index(
    mapper: HashMap<Vec<u8>, Vec<Arc<[u8]>>>,
    index_writer: &mut BufWriter<File>,
    dedup_writer: &mut BufWriter<File>,
) {
    // INFO: ensuring index is filled up -> will write directly in bytes
    if !mapper.is_empty() {
        let mut count: u32 = 0;

        // INFO: grabs each group and writes 1 sequence per group + an index of all records
        // INFO: pointing to that sequence with the following format:
        // INFO: id_len id seq_len n_headers read_chr [read_ids]
        // INFO: where each [read_id] follows the format: id orf subseq_orf start end
        for (seq, records) in mapper {
            index_writer.write_all(&count.to_be_bytes()).unwrap();

            let n_headers = records.len() as u32;
            index_writer.write_all(&n_headers.to_be_bytes()).unwrap();

            let _ = writeln!(dedup_writer, ">{}\n{}", count, from_utf8(&seq).unwrap());

            // INFO: we do not need to write chr for each record -> assuming same chr for all records
            // INFO: consider the following read names:
            //  >R10589_chr7__FC28#TC0#PA0#PR0#IY896_ORF.6_[4443-4650](+)_type:complete_length:204_frame:1_start:CTG_stop:TAG
            //  >R10589_chr7__FC28#TC0#PA0#PR0#IY896_ORF.6_[4443-4650](+)_type:complete_length:204_frame:1_start:CTG_stop:TAG@1
            // INFO: their indexed names will be -> 10589 7 6 0 4443 4650 and 10589 7 6 1 4443 4650
            let mut chr_count = 0;
            for record in records {
                if let Some((read_id, read_chr, read_orf, read_subseq_orf, start, end)) =
                    get_read_encoding(&record)
                {
                    // INFO: only writing chr for 1st record
                    if chr_count < 1 {
                        index_writer.write_all(&[read_chr.len() as u8]).unwrap();
                        index_writer.write_all(&read_chr).unwrap();
                    }
                    index_writer.write_all(&read_id.to_be_bytes()).unwrap();
                    index_writer.write_all(&read_orf.to_be_bytes()).unwrap();
                    index_writer
                        .write_all(&read_subseq_orf.to_be_bytes())
                        .unwrap();
                    index_writer.write_all(&start.to_be_bytes()).unwrap();
                    index_writer.write_all(&end.to_be_bytes()).unwrap();

                    chr_count += 1;
                } else {
                    panic!("Could not parse header: {:?}", std::str::from_utf8(&record));
                }
            }

            count += 1;
        }
    } else {
        error!("ERROR: No sequences found in the FASTA file!");
        std::process::exit(1);
    }
}

/// Parses a FASTA header byte slice to extract encoded read information.
///
/// The expected header format is designed to contain specific delimited information
/// about the read, chromosome, ORF, subsequence ORF, and genomic coordinates.
///
/// Example Header Format:
/// `R<read_id>_chr<chr_num>__FC...#TC...#PA...#PR...#IY..._ORF.<orf_num> [<start>-<end>](<strand>)_{...}@<subseq_orf>`
///
/// # Arguments
///
/// * `header` - A byte slice representing the FASTA header.
///
/// # Returns
///
/// An `Option` containing a tuple `(u16, Vec<u8>, u16, u16, u32, u32)` if parsing is successful.
/// The tuple elements are:
/// - `read_id`: The parsed read ID.
/// - `chr`: The chromosome name as a `Vec<u8>`.
/// - `orf`: The parsed ORF number.
/// - `subseq_orf`: The parsed subsequence ORF number (defaults to 0 if not present).
/// - `start`: The parsed start coordinate.
/// - `end`: The parsed end coordinate.
///
/// Returns `None` if the initial conversion of the header to a UTF-8 string fails,
/// or if any critical part of the header parsing (e.g., splitting by '_',
/// parsing numeric values, stripping prefixes) encounters an unrecoverable error.
///
/// # Panics
///
/// This function will panic if:
/// - The header cannot be converted to a UTF-8 string.
/// - The header does not have enough parts to be parsed.
/// - Any of the required numeric components (ID, ORF, start, end) fail to parse.
/// - The chromosome prefix "chr" is missing.
/// - The coordinate string parts (e.g., `[start-end](+)`) cannot be extracted or split.
///
/// # Example
///
/// ```rust, no_run
/// let header = b"R123_chrX__FC#TC#PA#PR#IY_ORF.5 [100-200](+)_type_length_frame_start_stop@1".to_vec();
/// let encoding = get_read_encoding(&header).unwrap();
/// assert_eq!(encoding.0, 123);
/// assert_eq!(encoding.1, b"X".to_vec());
/// assert_eq!(encoding.2, 5);
/// assert_eq!(encoding.3, 1);
/// assert_eq!(encoding.4, 100);
/// assert_eq!(encoding.5, 200);
/// ```
pub fn get_read_encoding(header: &[u8]) -> Option<(u16, Vec<u8>, u16, u16, u32, u32)> {
    // WARN: cannonical -> R6713_chr16__FC48#TC40#PA0#PR0#IY876
    let header = std::str::from_utf8(header).ok().unwrap_or_else(|| {
        panic!("ERROR: failed to convert header to string");
    });

    // R6456_chr16__FC42#TC48#PA0#PR0#IY907_ORF.247 [44369-44519](+) type:complete length:147 frame:3 start:CTG stop:TAA
    let parts: Vec<&str> = header.split('_').collect();

    if parts.len() < 2 {
        panic!(
            "ERROR: header does not have enough parts to parse: {}",
            &header
        );
    }

    let id = parts[0]
        .strip_prefix('R')?
        .parse::<u16>()
        .ok()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse ID from header: {}", header);
        });
    let chr = parts[1]
        .strip_prefix("chr")
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse chromosome from header: {}", header);
        })
        .as_bytes()
        .to_vec();
    let subseq = parts
        .last()?
        .split('@')
        .nth(1)
        .and_then(|s| s.parse::<u16>().ok())
        .unwrap_or(0);
    let orf = parts
        .get(4)
        .unwrap_or(&"ORF.0")
        .split(" ")
        .next()
        .unwrap_or(&"ORF.0")
        .strip_prefix("ORF.")
        .unwrap_or("0") // WARN: enforcing a default value
        .parse::<u16>()
        .ok()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse ORF from header: {}", header);
        });

    // INFO: ORF.247 [44369-44519](+) type:complete length:147 frame:3 start:CTG stop:TAA
    let coords = parts
        .get(4)
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse coordinates from header: {}", header);
        })
        .split(' ')
        .nth(1)
        .and_then(|s| s.strip_prefix('[')) // 44369-44519](+)
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse coordinates from header: {}", header);
        })
        .split(']')
        .next()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse coordinates from header: {}", header);
        })
        .split('-')
        .collect::<Vec<&str>>();

    let start = coords
        .get(0)
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or_else(|| {
            panic!(
                "ERROR: failed to parse start coordinate from header: {}",
                header
            );
        });
    let end = coords
        .get(1)
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or_else(|| {
            panic!(
                "ERROR: failed to parse end coordinate from header: {}",
                header
            );
        });

    return Some((id, chr, orf, subseq, start, end));
}

/// Reads an index file generated by `make_index` into a `HashMap`.
///
/// The index file is expected to contain a serialized representation of
/// deduplicated sequence IDs and associated original record metadata.
///
/// The structure of the binary index file is described in the `make_index` documentation.
///
/// # Arguments
///
/// * `index` - A `PathBuf` representing the path to the index file.
///
/// # Returns
///
/// A `HashMap<u32, Vec<(u16, Vec<u8>, u16, u16, u32, u32)>>` where:
/// - Keys are the unique `group_id`s (from the deduplicated sequences).
/// - Values are vectors of tuples, each tuple containing:
///     - `read_id`: The ID of the original read.
///     - `chr_bytes`: The chromosome name as a `Vec<u8>`.
///     - `orf`: The ORF number.
///     - `subseq`: The subsequence ORF number.
///     - `start`: The start coordinate from the original header.
///     - `end`: The end coordinate from the original header.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to open the index file.
/// - It encounters an `io::Error` while reading from the file,
///   unless it's an EOF error indicating the end of the file.
/// - The byte arrays read from the file cannot be converted back
///   to their respective integer types.
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
/// use std::fs::File;
/// use std::io::{BufWriter, Cursor};
/// use std::collections::HashMap;
/// use std::sync::Arc;
///
/// // Simulate creating an index file
/// let index_path = PathBuf::from("temp_read.index");
/// {
///     let mut writer = BufWriter::new(File::create(&index_path).unwrap());
///     // Write a dummy entry for group_id 0
///     writer.write_all(&0u32.to_be_bytes()).unwrap(); // group_id
///     writer.write_all(&1u32.to_be_bytes()).unwrap(); // n_headers
///     writer.write_all(&1u8.to_be_bytes()).unwrap();  // chr_len ('X')
///     writer.write_all(b"X").unwrap();                 // chr_bytes
///     writer.write_all(&123u16.to_be_bytes()).unwrap(); // read_id
///     writer.write_all(&5u16.to_be_bytes()).unwrap();  // orf
///     writer.write_all(&1u16.to_be_bytes()).unwrap();  // subseq
///     writer.write_all(&100u32.to_be_bytes()).unwrap(); // start
///     writer.write_all(&200u32.to_be_bytes()).unwrap(); // end
/// }
///
/// let index_data = read_index(&index_path);
/// assert!(index_data.contains_key(&0));
/// let records = index_data.get(&0).unwrap();
/// assert_eq!(records.len(), 1);
/// assert_eq!(records[0].0, 123);
/// assert_eq!(records[0].1, b"X".to_vec());
/// assert_eq!(records[0].2, 5);
/// assert_eq!(records[0].3, 1);
/// assert_eq!(records[0].4, 100);
/// assert_eq!(records[0].5, 200);
///
/// // Clean up the dummy file
/// std::fs::remove_file(&index_path).unwrap();
/// ```
///
pub fn read_index(index: &PathBuf) -> HashMap<u32, Vec<(u16, Vec<u8>, u16, u16, u32, u32)>> {
    let mut reader = BufReader::new(
        File::open(index).unwrap_or_else(|e| panic!("ERROR: failed to open index -> {e}")),
    );

    let mut result = HashMap::new();

    loop {
        let mut group_id_buf = [0u8; 4];
        if reader.read_exact(&mut group_id_buf).is_err() {
            break;
        }
        let group_id = u32::from_be_bytes(group_id_buf);

        // 4. Read n_headers
        let mut n_headers_buf = [0u8; 4];
        reader.read_exact(&mut n_headers_buf).unwrap();
        let n_headers = u32::from_be_bytes(n_headers_buf);

        // 5. Read chromosome
        let mut chr_len_buf = [0u8; 1];
        reader.read_exact(&mut chr_len_buf).unwrap();
        let chr_len = chr_len_buf[0] as usize;

        let mut chr_buf = vec![0u8; chr_len];
        reader.read_exact(&mut chr_buf).unwrap();
        let chr = chr_buf;

        // 6. Read n_headers records
        let mut records = Vec::with_capacity(n_headers as usize);
        for _ in 0..n_headers {
            let mut read_id_buf = [0u8; 2];
            reader.read_exact(&mut read_id_buf).unwrap();
            let read_id = u16::from_be_bytes(read_id_buf);

            let mut orf_buf = [0u8; 2];
            reader.read_exact(&mut orf_buf).unwrap();
            let orf = u16::from_be_bytes(orf_buf);

            let mut subseq_buf = [0u8; 2];
            reader.read_exact(&mut subseq_buf).unwrap();
            let subseq = u16::from_be_bytes(subseq_buf);

            let mut start_buf = [0u8; 4];
            reader.read_exact(&mut start_buf).unwrap();
            let start = u32::from_be_bytes(start_buf);

            let mut end_buf = [0u8; 4];
            reader.read_exact(&mut end_buf).unwrap();
            let end = u32::from_be_bytes(end_buf);

            records.push((read_id, chr.clone(), orf, subseq, start, end));
        }

        // INFO: { index_id : [(read_id: u16, chr_bytes: [u8; chr_len], orf: u16, subseq_orf: u16)] }
        // INFO: { 0 : [(read_id: 5903, chr_bytes: [16, 32], orf: 1, subseq_orf: 3)] }
        result.insert(group_id, records);
    }

    result
}
