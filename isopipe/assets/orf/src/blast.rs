use config::bed_to_custom_struct_collection;
use dashmap::DashMap;
use hashbrown::HashMap;
use log::{error, warn};
use packbed::{reader as bed_reader, record::GenePred};
use rayon::prelude::*;

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::str::from_utf8;
use std::sync::Arc;

use crate::{cli::BlastArgs, utils::*};

const ORF_PEP: &str = "orfs.pep.fa";
const ORF_BED: &str = "orfs.bed";
const RESULT: &str = "dmd.result";
const HEADER_REGEX: &str = r"\[(\d+)-(\d+)\]\(([+-])\)";

/// Runs a complete BLAST analysis pipeline, including ORF prediction, deduplication,
/// and alignment against a DIAMOND database.
///
/// This function orchestrates the entire process:
/// 1. Predicts Open Reading Frames (ORFs) from the input FASTA file using `orfipy`.
/// 2. Deduplicates the predicted ORFs based on length and percentage.
/// 3. Performs a protein-protein BLAST search using DIAMOND against the specified database.
/// 4. Processes and inflates the BLAST results, writing them to an output file.
///
/// # Arguments
///
/// * `args` - A `BlastArgs` struct containing all the necessary arguments for the BLAST run,
///            including paths to input files, output directories, and various
///            parameters for `orfipy` and DIAMOND.
///
/// # Panics
///
/// This function will panic if any of the underlying commands (`orfipy`, `diamond`) fail
/// to execute, or if there are issues with file I/O or data parsing.
///
/// # Example
///
/// ```rust, no_run
/// run_blast(args);
/// ```
pub fn run_blast(args: BlastArgs) {
    let pep = orfipy(&args.common.fasta, &args.common.outdir, &args.orfipy);

    let (dedup, index) = deduplicate(
        &pep,
        true,
        args.orf_min_len,
        args.orf_min_percent,
        "M".as_bytes(),
        SeqType::AminoAcid,
    );

    diamond(&dedup, &args.database, &index, &args.common.alignments);
}

/// Predicts Open Reading Frames (ORFs) from a FASTA file using `orfipy`.
///
/// This function creates a temporary directory for `orfipy`'s output, constructs
/// and executes the `orfipy` command, and returns the path to the generated
/// protein FASTA file.
///
/// # Arguments
///
/// * `fasta` - A `PathBuf` representing the path to the input FASTA file.
/// * `outdir` - A `PathBuf` representing the base output directory where `orfipy`'s
///              temporary directory will be created.
/// * `executable` - A `PathBuf` representing the path to the `orfipy` executable.
///
/// # Returns
///
/// A `PathBuf` pointing to the generated protein FASTA file by `orfipy`.
///
/// # Panics
///
/// This function will panic if it cannot create the necessary output directory
/// or if the `orfipy` command fails to execute.
///
/// # Example
///
/// ```rust, no_run
/// let output = orfipy(&fasta_path, &output_dir, &orfipy_executable);
/// ```
fn orfipy(fasta: &PathBuf, outdir: &PathBuf, executable: &PathBuf) -> PathBuf {
    let dir = outdir.join("orf");

    std::fs::create_dir_all(&dir)
        .unwrap_or_else(|e| panic!("ERROR: could not create directory -> {e}!"));

    let cmd = format!(
        "{} {} --pep {} --bed {} --partial-5 --partial-3 --include-stop --min 100 --ignore-case --outdir {}",
        executable.display(),
        fasta.display(),
        ORF_PEP,
        ORF_BED,
        &dir.display()
    );

    std::process::Command::new("bash")
        .arg("-c")
        .arg(cmd)
        .status()
        .unwrap_or_else(|e| panic!("ERROR: failed to execute orfipy command -> {e}"));

    dir.join(ORF_PEP)
}

/// Performs a protein-protein BLAST search using DIAMOND and processes the results.
///
/// This function executes a DIAMOND `blastp` command, filters and processes the
/// resulting alignments, and inflates the results with additional information
/// from the original gene predictions. It then writes the final, processed
/// BLAST records to an output file. Unused IDs (those without a DIAMOND hit)
/// are also reported with a specific tag [#DM].
///
/// # Arguments
///
/// * `dedup` - A `PathBuf` representing the path to the deduplicated protein sequences
///             (query sequences for DIAMOND).
/// * `database` - A `PathBuf` representing the path to the DIAMOND database.
/// * `index` - A `PathBuf` representing the path to an index file containing
///             information about the original sequences, used for result inflation.
/// * `alignments` - A `PathBuf` representing the path to a BED file containing
///                  gene predictions, used to map absolute CDS coordinates.
///
/// # Panics
///
/// This function will panic if:
/// - The DIAMOND command fails to execute.
/// - It cannot read the BLAST predictions file or the BED alignment file.
/// - It cannot create the output file for BLAST results.
/// - It encounters issues parsing data from the BLAST output or the index.
/// - It cannot find a corresponding query ID or chromosome in the provided indices
///   or records during result inflation.
///
/// # Example
///
/// ```rust, no_run
/// diamond(&dedup_path, &database_path, &index_path, &alignments_path);
/// ```
fn diamond(dedup: &PathBuf, database: &PathBuf, index: &PathBuf, alignments: &PathBuf) {
    let dmd = dedup.with_extension("diamond");
    let cmd = format!(
        "diamond blastp --query {} --db {} --out {} --outfmt 6 qseqid pident qlen slen length qstart qend sstart send evalue --threads 8 --sensitive -e 1e-10",
        dedup.display(),
        database.display(),
        dmd.display()
    );

    std::process::Command::new("bash")
        .arg("-c")
        .arg(cmd)
        .status()
        .unwrap_or_else(|e| panic!("ERROR: failed to execute diamond command -> {e}"));

    let mut index = read_index(&index);
    let predictions = bed_reader(&dmd)
        .unwrap_or_else(|e| panic!("ERROR: failed to read blast predictions file -> {e}"));

    let mut writer = BufWriter::new(
        File::create(dmd.with_extension(RESULT)).unwrap_or_else(|e| {
            panic!("ERROR: failed to create output file for blast results -> {e}");
        }),
    );

    let accumulator = DashMap::new();

    let records = bed_to_custom_struct_collection::<GenePred>(
        bed_reader(alignments)
            .unwrap_or_else(|e| panic!("ERROR: failed to read BED file -> {e}"))
            .into(),
        config::BedColumn::Name,
        config::BedOperation::SplitName("__", 0), // INFO: R9834_chr16__FC37#TC0#PA0#PR0#IY887
    )
    .unwrap_or_else(|e| panic!("ERROR: failed construct BED to GenePred collection -> {e}"));

    // INFO: filtering repeated blast hits by percent_identity -> preserving best
    predictions
        .par_lines()
        .filter(|line| !line.starts_with('#'))
        .for_each(|line| {
            let parts: Vec<&str> = line.split('\t').collect();

            // INFO: 16 sp|Q9QX47|SON_MOUSE 100 497 0 0 42 538 1089 1585 1.20e-163 515
            let id = parts[0].parse::<u32>().unwrap_or_else(|_| {
                panic!("ERROR: failed to parse ID from line: {}", line);
            });

            let data = BlastRecord::from_parts(&parts);

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

    // INFO: inflate results!
    accumulator.iter_mut().for_each(|mut record| {
        let (id, data) = record.pair_mut();

        // INFO: unpacking index reference -> queries
        // INFO: for each query all blast records
        // INFO: { index_id : [(read_id: u16, chr_bytes: [u8; chr_len], orf: u16, subseq_orf: u16, seq_len: usize)] }
        // INFO: { 0 : [(read_id: 5903, chr_bytes: [16, 32], orf: 1, subseq_orf: 3, seq_len: 350)] }
        let queries = index.get(id).unwrap_or_else(|| {
            panic!("ERROR: no queries found for ID: {}", id);
        });

        for query in queries.into_iter() {
            let (read_id, chr, orf, subseq_orf, start, end) = query;

            let chr = format!("chr{}", from_utf8(&chr).unwrap());
            let cannonical_id = format!("R{}_{}", read_id, chr);
            let query_id = format!("{}.p{}@{}", cannonical_id, orf, subseq_orf);

            // INFO: retrieving the reference gene prediction record
            let (orf_start, orf_end) = records
                .get_mut(&chr)
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: chromosome from {} not found in sequences -> {}!",
                        cannonical_id, chr
                    );
                })
                .get_mut(&cannonical_id)
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: id not found in BED, this is a bug -> {}!",
                        cannonical_id
                    );
                })
                .map_absolute_cds(*start as u64, *end as u64)
                .unwrap_or_default();

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

            data.set_id(query_id);

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
    index.iter().for_each(|(id, queries)| {
        for query in queries {
            let (read_id, chr, orf, subseq_orf, start, end) = query;

            let chr = format!("chr{}", from_utf8(&chr).unwrap());
            let cannonical_id = format!("R{}_{}", read_id, chr);
            let query_id = format!("{}.p{}@{}#DM", cannonical_id, orf, subseq_orf);

            // INFO: retrieving the reference gene prediction record
            let (orf_start, orf_end) = records
                .get_mut(&chr)
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: chromosome from {} not found in sequences -> {}!",
                        cannonical_id, chr
                    );
                })
                .get_mut(&cannonical_id)
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: id not found in BED, this is a bug -> {}!",
                        cannonical_id
                    );
                })
                .map_absolute_cds(*start as u64, *end as u64)
                .unwrap_or_default();

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

/// Represents a single BLAST alignment record.
#[derive(Debug, Clone, PartialEq)]
pub struct BlastRecord {
    pub blast_id: String,         // ID of the blast record
    pub blast_idx_id: u32,        // Indexed ID of the blast record
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
    /// * `parts` - A slice of string slices, where each element represents
    ///             a column from the BLAST output.
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

        let blast_idx_id = parts[0].parse::<u32>().unwrap_or_else(|_| {
            panic!("ERROR: failed to parse blast ID from parts: {:?}", parts);
        });
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
        let blast_e_value = if blast_e_value == 0.0 && parts[9] != "0.0" {
            // INFO: represent it with the minimum positive value
            f64::MIN_POSITIVE // INFO: ~2.225074e-308
        } else {
            blast_e_value
        };

        let blast_offset = parts[7].parse::<i32>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast offset from parts: {:?}",
                parts
            )
        }) - parts[5].parse::<i32>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast offset from parts: {:?}",
                parts
            );
        });
        let blast_alignment_len = parts[4].parse::<u32>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast offset from parts: {:?}",
                parts
            )
        });

        let percent_aligned = blast_alignment_len as f32
            / parts[2].parse::<u32>().unwrap_or_else(|_| {
                panic!(
                    "ERROR: failed to parse blast length from parts: {:?}",
                    parts
                );
            }) as f32
            * 100.0;

        Self {
            blast_id: String::new(), // INFO: set on the fly
            blast_idx_id,
            blast_pid,
            blast_e_value,
            blast_offset,
            blast_alignment_len,
            percent_aligned: percent_aligned, // INFO: set on the fly
        }
    }

    /// Sets the `blast_id` for the `BlastRecord`.
    ///
    /// This method is used to assign a specific identifier to the record after
    /// its initial creation, typically combining information from the original
    /// sequence and its genomic location.
    ///
    /// # Arguments
    ///
    /// * `id` - A `String` representing the ID to be set for the record.
    pub fn set_id(&mut self, id: String) {
        self.blast_id = id;
    }
}

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
) -> (PathBuf, PathBuf) {
    let seqs =
        parse_fa(fasta).unwrap_or_else(|e| panic!("ERROR: failed to parse FASTA file -> {e}"));

    let mut mapper = HashMap::new();
    let mut index = create_index(fasta);
    let mut dedup = create_fasta(fasta, "dedup.fa")
        .unwrap_or_else(|| panic!("ERROR: could not create file {:?}", fasta));

    let regex = regex::Regex::new(HEADER_REGEX).unwrap_or_else(|e| {
        panic!("ERROR: failed to compile regex for header parsing -> {e}");
    });

    // INFO: loops through sequences and populates mapper
    for (header, seq) in seqs.iter() {
        let len = seq.len();
        let mut key = Vec::new();
        let header = header.to_string();

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
            )
        }

        let record = Arc::<[u8]>::from(header.into_bytes());
        mapper.entry(key).or_insert(Vec::new()).push(record);
    }

    let _ = make_index(mapper, &mut index, &mut dedup);

    // INFO: cleaning footprint on the fly
    // if !seqs.is_empty() {
    //     std::fs::remove_file(fasta).expect("ERROR: failed to remove original FASTA file");
    // }

    return (
        fasta.with_extension("dedup.fa"),
        fasta.with_extension("dedup.index"),
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
fn split_record(
    header: &String,
    seq: &Vec<u8>,
    seq_length: usize,
    min_len: usize,
    min_percent: f32,
    mapper: &mut HashMap<Vec<u8>, Vec<Arc<[u8]>>>,
    needle: &[u8],        // b"ATG" or b"M"
    seq_type: SeqType,    // determines scanning logic
    regex: &regex::Regex, // regex for parsing header
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
                        let (orig_start, orig_end, strand) = split_header(&header, regex)
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
                        let cannonical_id = header.split(' ').next().unwrap_or_else(|| {
                            panic!(
                                "ERROR: failed to parse cannonical ID from header: {}",
                                header
                            );
                        });
                        let sub_id = format!(
                            "{} [{}-{}]({})@{}",
                            cannonical_id, nested_start, nested_end, strand, orf_count
                        );

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
fn split_header<'a>(header: &'a str, capture: &regex::Regex) -> Option<(usize, usize, char)> {
    let caps = capture.captures(header)?;

    let start = caps[1].parse().ok()?;
    let end = caps[2].parse().ok()?;
    let strand = caps[3].chars().next()?;
    Some((start, end, strand))
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
fn make_index(
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
fn get_read_encoding(header: &[u8]) -> Option<(u16, Vec<u8>, u16, u16, u32, u32)> {
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
