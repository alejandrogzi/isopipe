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
use log::error;
use packbed::reader as bed_reader;
use rayon::prelude::*;

use std::io::Write;
use std::path::PathBuf;
use std::str::from_utf8;
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
    seq_type: SeqType,
    mode: &Mode,
    regex: &regex::Regex,
) -> (PathBuf, PathBuf, HashMap<usize, Vec<String>>) {
    let seqs =
        parse_fa(fasta).unwrap_or_else(|e| panic!("ERROR: failed to parse FASTA file -> {e}"));

    let mut mapper = HashMap::new();
    let mut helper = HashMap::new(); // INFO: idx -> name

    // INFO: loops through sequences and populates mapper
    for (header, seq) in seqs.iter() {
        let len = seq.len();
        let mut key = Vec::with_capacity(len);

        // INFO: >R10589_chr7__FC28#TC0#PA0#PR0#IY896_ORF.87 [4632-4770](-) [...]
        // INFO: >0_ORF.87 [4632-4770](-) [...] -> [if indexed during extract]
        // INFO: ends up being -> >0_ORF.87 [4632-4770](-)
        let hdr = header.to_string();
        let parts = hdr.split(' ').take(2).collect::<Vec<&str>>();
        let header = parts.join(" "); // INFO: orfipy headers!

        // INFO: getting id and strand from header!
        let _idx = parts[0].split('_').take(1).next().unwrap_or_else(|| {
            panic!(
                "ERROR: could not get index number from name -> {:?}",
                header
            )
        });

        let strand = parts[1]
            .split('(')
            .nth(1)
            .unwrap_or_else(|| panic!("ERROR: could not get strand from parts -> {:?}", parts))
            .strip_suffix(')')
            .unwrap();

        // WARN: skipping on reverse stranded preditions
        // WARN: orfipy will take the reference strand as '+'!
        if strand != "+" {
            continue;
        };

        // INFO: 0_ORF.87 [4632-4770](-) -> { 0 : [ 0_ORF.87_[4632-4770](-) ] }
        // INFO: 0_ORF.95 [4632-4770](-) -> { 0 : [ 0_ORF.87_[4632-4770](-), 0_ORF.95_[4632-4770](-) ] }
        match mode {
            Mode::Indexed => {
                helper
                    .entry(_idx.parse::<usize>().unwrap_or_else(|e| {
                        panic!(
                            "ERROR: could not parse number from name -> {:?}. {e}",
                            header
                        )
                    }))
                    .or_insert_with(Vec::new)
                    .push(header.clone().replace(" ", "_"));
            }
            _ => {}
        }

        for &b in seq {
            if b != b'\n' {
                key.push(b);
            }
        }

        if do_nesting {
            split_record(
                &header,
                seq,
                len,
                min_len,
                min_percent,
                &mut mapper,
                pattern,
                seq_type,
                &regex,
                &mut helper,
                mode,
            )
        }

        let record = Arc::<[u8]>::from(header.replace(" ", "_").into_bytes());
        mapper.entry(key).or_insert(Vec::new()).push(record);
    }

    let mut dedup = create_fasta(fasta, "dedup.fa")
        .unwrap_or_else(|| panic!("ERROR: could not create file {:?}", fasta));

    match mode {
        Mode::Indexed => {
            if !mapper.is_empty() {
                for (seq, records) in mapper {
                    // INFO: we expect to only have records of len 1
                    for rc in records {
                        let _ = writeln!(
                            dedup,
                            ">{}\n{}",
                            from_utf8(&rc).unwrap(),
                            from_utf8(&seq).unwrap()
                        );
                    }
                }
            } else {
                error!("ERROR: mapper is empty -> no records were read!");
                std::process::exit(1);
            }
        }
        Mode::Raw => {
            let mut index = create_index(fasta);
            let _ = crate::blast::cannonical::make_index(mapper, &mut index, &mut dedup);
        }
    }

    return (
        fasta.with_extension("dedup.fa"),
        fasta.with_extension("dedup.index"),
        helper,
    );
}

/// An enum representing the type of biological sequence.
///
/// This is used to determine the appropriate logic for sequence processing,
/// such as searching for start codons in nucleotides or start residues in amino acids.
#[derive(Debug, Clone, Copy)]
pub enum SeqType {
    Nucleotide, // search codons like ATG
    AminoAcid,  // search residues like M
}

/// Splits a FASTA record into potentially multiple sub-sequences based on a
/// specified start signal (needle), respecting codon logic for nucleotide sequences.
///
/// This function iterates through a sequence and, if `do_nesting` is enabled,
/// identifies subsequences starting with the `needle` (e.g., 'M' for proteins
/// or 'ATG' for nucleotides). These subsequences are then added to the `mapper`
/// if they meet the minimum length and percentage criteria.
///
/// # Arguments
///
/// * `header` - A `String` reference to the header of the original FASTA record.
/// * `seq` - A `Vec<u8>` reference to the byte representation of the sequence.
/// * `seq_length` - The total length of the original sequence.
/// * `min_len` - The minimum length for a subsequence to be considered valid.
/// * `min_percent` - The minimum percentage of the original sequence length
///                   that a subsequence must represent to be considered valid.
/// * `mapper` - A mutable `HashMap` that stores the deduplicated sequences.
///              Keys are sequence bytes, and values are vectors of ARC-wrapped
///              headers of the sequences that map to that key.
/// * `needle` - A byte slice representing the start signal to search for
///              (e.g., `b"ATG"` or `b"M"`).
/// * `seq_type` - A `SeqType` enum indicating whether the sequences are nucleotides
///                or amino acids, influencing the splitting and length calculation logic.
/// * `regex` - A reference to a compiled `regex::Regex` used for parsing information
///             (like original start, end, and strand) from the header for amino acid sequences.
///
/// # Panics
///
/// This function will panic if:
/// - It encounters an unknown strand character during header parsing for amino acid sequences.
/// - It fails to parse necessary information (e.g., start, end) from the header
///   when processing amino acid sequences.
pub fn split_record(
    header: &String,
    seq: &Vec<u8>,
    seq_length: usize,
    min_len: usize,
    min_percent: f32,
    mapper: &mut HashMap<Vec<u8>, Vec<Arc<[u8]>>>,
    needle: &[u8],        // b"ATG" or b"M"
    seq_type: SeqType,    // determines scanning logic
    regex: &regex::Regex, // regex for parsing header
    helper: &mut HashMap<usize, Vec<String>>,
    mode: &Mode,
) {
    // INFO: always write the original full sequence
    let mut orf_count = 0;

    match seq_type {
        SeqType::Nucleotide => {
            let codon_len = 3;
            let mut pos = codon_len; // INFO: skip first codon
            while pos + codon_len <= seq_length {
                if &seq[pos..pos + codon_len] == needle {
                    let len_remaining = seq_length - pos;
                    let percent = (len_remaining as f32 / 3_f32) / (seq_length as f32 / 3_f32);

                    if (len_remaining as f32 / 3_f32) >= min_len as f32 && percent >= min_percent {
                        orf_count += 1;
                        let sub_seq = &seq[pos..];
                        let sub_id = format!("{}@{}", header, orf_count);

                        let mut inner_key = Vec::with_capacity(sub_seq.len());

                        for &b in sub_seq {
                            if b != b'\n' {
                                inner_key.push(b);
                            }
                        }

                        let record = Arc::<[u8]>::from(sub_id.clone().into_bytes());
                        mapper.entry(inner_key).or_insert(Vec::new()).push(record);
                    }
                }
                pos += codon_len;
            }
        }
        SeqType::AminoAcid => {
            for (pos, &aa) in seq.iter().enumerate().skip(1) {
                // INFO: default needle is 'M' -> start signal
                if aa == needle[0] {
                    let len_remaining = seq_length - pos;
                    let percent = len_remaining as f32 / seq_length as f32;

                    if len_remaining >= min_len && percent >= min_percent {
                        orf_count += 1;
                        let sub_seq = &seq[pos..];
                        let (orig_start, orig_end, strand) =
                            split_header(&header, regex) // >0_ORF.87 [4632-4770](-)
                                .unwrap_or_else(|| {
                                    panic!("ERROR: failed to parse header: {}", header);
                                });

                        let (nested_start, nested_end) = match strand {
                            '+' => {
                                let start = orig_start + pos * 3;
                                (start, orig_end)
                            }
                            '-' => {
                                let start = orig_end - pos * 3;
                                (orig_start, start)
                            }
                            _ => panic!("ERROR: unknown strand -> {strand} in header: {header}"),
                        };

                        // INFO: >R10589_chr7__FC28#TC0#PA0#PR0#IY896_ORF.87 [4632-4770](-) [...]
                        // INFO: >0_ORF.87 [4632-4770](-) [...] -> [if indexed during extract]
                        let cannonical_id = header.split(' ').next().unwrap_or_else(|| {
                            panic!(
                                "ERROR: failed to parse cannonical ID from header: {}",
                                header
                            );
                        });
                        let _idx = cannonical_id.split('_').next().unwrap_or_else(|| {
                            panic!(
                                "ERROR: failed to parse idx ID from header: {}",
                                &cannonical_id
                            );
                        });
                        let sub_id = format!(
                            "{}_[{}-{}]({})@{}",
                            cannonical_id, nested_start, nested_end, strand, orf_count
                        );

                        match mode {
                            Mode::Indexed => {
                                helper
                                    .entry(_idx.parse::<usize>().unwrap_or_else(|e| {
                                        panic!(
                                            "ERROR: could not parse number from name -> {:?}. {e}",
                                            header
                                        )
                                    }))
                                    .or_insert_with(Vec::new)
                                    .push(sub_id.clone());
                            }
                            _ => {}
                        }

                        let mut inner_key = Vec::with_capacity(sub_seq.len());

                        for &b in sub_seq {
                            if b != b'\n' {
                                inner_key.push(b);
                            }
                        }

                        let record = Arc::<[u8]>::from(sub_id.clone().into_bytes());
                        mapper.entry(inner_key).or_insert(Vec::new()).push(record);
                    }
                }
            }
        }
    }
}

/// Parses a DIAMOND BLAST output file and filters for the best hits.
///
/// This function reads a DIAMOND output file line by line, parallelizing the
/// parsing process. For each line (representing a hit), it creates a `BlastRecord`.
/// If multiple hits exist for the same query ID, it keeps only the one with the
/// highest percentage identity (`blast_pid`). The results are accumulated into a
/// `DashMap` for concurrent access.
///
/// # Arguments
///
/// * `diamond` - A `PathBuf` representing the path to the DIAMOND output file.
/// * `mode` - A reference to an `Mode` enum, indicating whether the IDs
///            should be parsed as `Indexed` (e.g., "0_ORF...") or `Raw` (e.g., "17").
/// * `regex` - A reference to a `regex::Regex` object used for parsing
///             specific fields within the `BlastRecord`.
///
/// # Returns
///
/// A `DashMap<u32, BlastRecord>` containing the best Blast hit for each unique query ID.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to read the `diamond` output file.
/// - It fails to parse a query ID or any other field from a line.
pub fn parse_predictions(
    diamond: &PathBuf,
    mode: &Mode,
    regex: &regex::Regex,
) -> DashMap<u32, BlastRecord> {
    let predictions = bed_reader(diamond)
        .unwrap_or_else(|e| panic!("ERROR: failed to read blast predictions file -> {e}"));

    let accumulator = DashMap::new();

    // INFO: filtering repeated blast hits by percent_identity -> preserving best
    predictions
        .par_lines()
        .filter(|line| !line.starts_with('#'))
        .for_each(|line| {
            let parts: Vec<&str> = line.split('\t').collect();

            // qseqid pident  qlen    slen   length qstart    qend   sstart   send     evalue
            //  17      97.2    142     357     141     1       141     217     357     5.09e-93
            // WARN: with extract indexing qseqid -> 0_ORF.1_[1-10](+)
            let id = match mode {
                Mode::Indexed => parts[0].split('_').collect::<Vec<&str>>()[0]
                    .parse::<u32>()
                    .unwrap_or_else(|_| {
                        panic!("ERROR: failed to parse ID from line: {}", line);
                    }),
                Mode::Raw => parts[0].parse::<u32>().unwrap_or_else(|_| {
                    panic!("ERROR: failed to parse ID from line: {}", line);
                }),
            };

            let data = BlastRecord::from_parts(&parts, mode, regex);

            // WARN: using a transition collection to retain the best blast record based on % aligned
            accumulator
                .entry(id)
                .and_modify(|existing_data: &mut BlastRecord| {
                    if data.blast_pid > existing_data.blast_pid {
                        *existing_data = data.clone();
                    }
                })
                .or_insert(data);
        });

    accumulator
}
