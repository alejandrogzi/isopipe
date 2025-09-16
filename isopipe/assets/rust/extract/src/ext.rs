//! Core module for extracting sequences using .2bit from a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for efficiently getting whole
//! sequences from .2bit genomic files using plain coordinates. Additionally,
//! it provides the option of deduplicate repeated entries in order to save
//! space.
//!
//! In short, every sequence is extracted from a .2bit and holded in memory
//! as plain bytes for every read in the query set [rev comp for neg strand].
//! The command line interface provides the option of specifying indexing, leading
//! to a one-pass deduplication step and the creation of an index where
//! simple integers map to read identifiers [all of them as plain bytes].
//! The process is heavily parallelized to offer fast performance on large datasets.

use config::{OverlapType, SCALE, Sequence, Strand};
use dashmap::DashMap;
use iso_polya::utils::get_sequences;
use log::debug;
use packbed::{GenePred, unpack};
use rayon::prelude::*;

use std::collections::{HashMap, hash_map::Entry};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use crate::cli::{ExtractArgs, SeqMode};

/// Extracts sequences from a 2bit genome file based on BED records,
/// chunks them, and writes them to temporary FASTA and BED files.
///
/// It performs the following steps:
/// 1. Creates a temporary directory for chunked output files.
/// 2. Unpacks the input BED file into `GenePred` structures, potentially
///    handling overlapping exons.
/// 3. Loads the entire genome sequences from the 2bit file into memory.
/// 4. Chunks the `GenePred` records by chromosome and then into smaller
///    sub-chunks based on `args.chunk_size`.
/// 5. For each chunk, it extracts the corresponding sequences from the genome
///    and writes them to a new FASTA file. The original BED entries for that
///    chunk are written to a new BED file.
///
/// # Arguments
///
/// * `args` - An `Args` struct containing:
///   - `bed`: Path to the input BED file.
///   - `twobit`: Path to the 2bit genome file.
///   - `output_dir`: Base directory for temporary chunked output.
///   - `dir_prefix`: Prefix for the temporary directory name.
///   - `suffix`: Suffix for the temporary directory name.
///   - `chunk_size`: Maximum number of records per chunk.
///   - `mode`: A boolean indicating the extraction mode (true for `Indexed`, false for `Raw`).
///
/// # Returns
///
/// A `Vec<(PathBuf, PathBuf)>` where each tuple contains the path to a
/// chunked FASTA file and its corresponding chunked BED file.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to create the temporary output directory.
/// - It fails to unpack the BED file.
/// - It fails to load sequences from the 2bit file.
/// - It fails to create chunk-specific directories or files.
/// - It fails to write sequences or BED lines to the chunked files.
/// - A chromosome is missing in the 2bit genome during sequence extraction.
pub fn extract(args: ExtractArgs) -> Vec<(PathBuf, PathBuf)> {
    let mode = ExtractMode::from(args.mode);

    log::info!(
        "INFO: Extracting mapped read sequences [{}] from .2bit file...",
        args.bed.display()
    );

    let tmp_dir = args
        .output_dir
        .join(format!("{}_{}", args.dir_prefix, args.suffix));
    std::fs::create_dir_all(&tmp_dir).unwrap_or_else(|e| {
        panic!(
            "ERROR: could not creat temporary directory in {} -> {e}",
            &tmp_dir.display()
        )
    });

    let bed = unpack::<GenePred, _>(vec![args.bed.clone()], OverlapType::Exon, false)
        .unwrap_or_else(|e| {
            panic!(
                "ERROR: could not unpack reads -> {}. {e}",
                args.bed.display()
            )
        });
    log::debug!("DEBUG: packed bed -> {bed:#?}");

    let (genome, _) = get_sequences(args.twobit.clone()).unwrap_or_else(|| {
        panic!(
            "ERROR: could not get sequences from .2bit -> {}",
            args.twobit.display()
        )
    });
    debug!("DEBUG: loaded genome with {} chromosomes", genome.len());

    // INFO: define the chunk size for parallel processing
    // INFO: if chunk size > bed records -> symlink
    let paths: Vec<(PathBuf, PathBuf)> = bed
        .into_par_iter()
        .flat_map(|(chrom, records)| {
            records
                .chunks(args.chunk_size)
                .enumerate()
                .map(move |(chunk_id, chunk)| (format!("{}:{}", chrom, chunk_id), chunk.to_vec()))
                .collect::<Vec<_>>()
        })
        .map(|(chunk_id, transcripts)| {
            let chunk_path = tmp_dir.join(&chunk_id);
            std::fs::create_dir_all(&chunk_path).unwrap();

            let chunk_fa = chunk_path.join(format!("tmp_chunk_{}.fa", &chunk_id));
            let chunk_bed = chunk_path.join(format!("tmp_chunk_{}.bed", &chunk_id));

            let writer_fa = BufWriter::new(File::create(&chunk_fa).unwrap());
            let writer_bed = BufWriter::new(File::create(&chunk_bed).unwrap());

            let chr = chunk_id.split(':').next().unwrap_or(&chunk_id);

            match mode {
                ExtractMode::Raw => {
                    raw(
                        transcripts,
                        &genome,
                        chr,
                        writer_fa,
                        writer_bed,
                        &args.seq_mode,
                    );
                }
                ExtractMode::Indexed => {
                    index(
                        transcripts,
                        &genome,
                        chr,
                        &chunk_path,
                        &chunk_id,
                        writer_fa,
                        writer_bed,
                        &args.seq_mode,
                    );
                }
            }

            (chunk_fa, chunk_bed)
        })
        .collect();

    return paths;
}

/// An enum representing the mode for sequence extraction.
///
/// This enum determines whether sequences are extracted directly ("Raw")
/// or if an indexing approach is used (though `Indexed`)
enum ExtractMode {
    Raw,
    Indexed,
}

impl ExtractMode {
    /// Creates an `ExtractMode` from a boolean value.
    ///
    /// # Arguments
    ///
    /// * `mode` - A boolean value. `true` maps to `ExtractMode::Indexed`,
    ///            `false` maps to `ExtractMode::Raw`.
    ///
    /// # Returns
    ///
    /// An `ExtractMode` variant.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// let raw_mode = ExtractMode::from(false);
    /// assert!(matches!(raw_mode, ExtractMode::Raw));
    ///
    /// let indexed_mode = ExtractMode::from(true);
    /// assert!(matches!(indexed_mode, ExtractMode::Indexed));
    /// ```
    fn from(mode: bool) -> Self {
        match mode {
            true => Self::Indexed,
            false => Self::Raw,
        }
    }
}

/// Retrieves a specific sequence from the genome based on chromosome and `GenePred` transcript information.
///
/// This function extracts a subsequence from the provided `genome` (a `DashMap` of chromosome sequences)
/// corresponding to the coordinates defined in the `transcript`. It handles both forward and
/// reverse strands, performing reverse complementation for reverse-strand transcripts.
///
/// # Arguments
///
/// * `genome` - A reference to a `DashMap<String, Vec<u8>>` containing chromosome sequences.
/// * `chr` - A string slice representing the chromosome name.
/// * `transcript` - A reference to a `GenePred` struct containing the transcript's
///                  genomic coordinates (start, end) and strand information.
///
/// # Returns
///
/// A `config::Sequence` object representing the extracted sequence.
///
/// # Panics
///
/// This function will panic if:
/// - The specified `chr` is not found as a key in the `genome` `DashMap`.
/// - The `SCALE` constant (used for reverse strand coordinate transformation)
///   is not defined or accessible.
fn get_sequence(
    genome: &DashMap<String, Vec<u8>>,
    chr: &str,
    transcript: &GenePred,
    seq_mode: &SeqMode,
) -> config::Sequence {
    let chr_seq = genome
        .get_mut(chr)
        .unwrap_or_else(|| panic!("ERROR: missing chromosome in .2bit -> {chr}"));

    log::debug!(
        "DEBUG: extracting sequence for {}:{}-{} ({:?})",
        chr,
        transcript.start,
        transcript.end,
        transcript.strand
    );

    let seq = match seq_mode {
        SeqMode::Genome => match transcript.strand {
            Strand::Forward => {
                Sequence::new(&chr_seq[transcript.start as usize..transcript.end as usize])
            }
            Strand::Reverse => Sequence::new(
                &chr_seq[(SCALE - transcript.end) as usize..(SCALE - transcript.start) as usize],
            )
            .reverse_complement(),
        },

        SeqMode::Exon => {
            // INFO: extract and concatenate exon sequences
            let mut exonic_seq: Vec<u8> = Vec::with_capacity(transcript.exon_len as usize);
            for (exon_start, exon_end) in &transcript.exons {
                match transcript.strand {
                    Strand::Forward => {
                        let start = *exon_start as usize;
                        let end = *exon_end as usize;
                        exonic_seq.extend_from_slice(&chr_seq[start..end]);
                    }
                    Strand::Reverse => {
                        let start = (SCALE - *exon_end) as usize;
                        let end = (SCALE - *exon_start) as usize;
                        let mut buf = chr_seq[start..end].to_vec();
                        __rev_complement_u8(&mut buf);

                        log::debug!(
                            "DEBUG: reversing exon seq -> {}:{}-{} ({})",
                            chr,
                            start,
                            end,
                            Sequence::new(&buf)
                        );

                        exonic_seq.extend_from_slice(&buf);
                    }
                }
            }

            let exonic = Sequence::new(&exonic_seq);
            log::debug!("DEBUG: extracted exonic seq -> {}", exonic);

            exonic
        }

        SeqMode::Intron => {
            // optional: for completeness, handle introns similarly
            let mut intronic_seq: Vec<u8> = Vec::new();
            for (intron_start, intron_end) in &transcript.introns {
                match transcript.strand {
                    Strand::Forward => {
                        let start = *intron_start as usize;
                        let end = *intron_end as usize;
                        intronic_seq.extend_from_slice(&chr_seq[start..end]);
                    }
                    Strand::Reverse => {
                        let start = (SCALE - *intron_end) as usize;
                        let end = (SCALE - *intron_start) as usize;
                        let mut buf = chr_seq[start..end].to_vec();
                        __rev_complement_u8(&mut buf);

                        intronic_seq.extend_from_slice(&buf);
                    }
                }
            }

            let intronic = Sequence::new(&intronic_seq);
            log::debug!("DEBUG: extracted intronic seq -> {}", intronic);

            intronic
        }
    };

    seq
}

/// Writes extracted raw sequences and their corresponding BED lines to respective files.
///
/// This function iterates through a vector of `GenePred` transcripts. For each transcript,
/// it extracts its sequence from the provided `genome` using `get_sequence` and writes
/// the sequence to a FASTA file. It also writes the original BED line of the transcript
/// to a BED file.
///
/// # Arguments
///
/// * `transcripts` - A `Vec<GenePred>` containing the transcript records to process.
/// * `genome` - A reference to a `DashMap<String, Vec<u8>>` containing chromosome sequences.
/// * `chr` - A string slice representing the chromosome name for the current set of transcripts.
/// * `writer_fa` - A mutable `BufWriter<File>` for writing FASTA entries.
/// * `writer_bed` - A mutable `BufWriter<File>` for writing BED entries.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to write to `writer_fa` or `writer_bed`.
/// - The `get_sequence` function panics (e.g., due to a missing chromosome).
fn raw(
    transcripts: Vec<GenePred>,
    genome: &DashMap<String, Vec<u8>>,
    chr: &str,
    mut writer_fa: BufWriter<File>,
    mut writer_bed: BufWriter<File>,
    seq_mode: &SeqMode,
) {
    for tx in transcripts {
        let seq = get_sequence(genome, chr, &tx, seq_mode);

        writeln!(writer_fa, ">{}\n{}", tx.name, seq).unwrap();
        writeln!(writer_bed, "{}", tx.line).unwrap();
    }
}

/// Indexes and deduplicates `GenePred` transcripts, writing unique sequences
/// to a FASTA file, a reduced BED file, and an index file.
///
/// This function processes a chunk of `GenePred` transcripts. It extracts the
/// sequence for each transcript, uses the sequence as a key to deduplicate
/// records, and maintains a mapping of original transcript IDs to a new,
/// sequential ID for unique sequences.
///
/// For each unique sequence:
/// - It writes the sequence to a FASTA file, with the new sequential ID as the header.
/// - It writes a modified BED line to a "reduced" BED file, where the name
///   field is replaced by the new sequential ID.
/// - It records the mapping of original IDs to the new sequential ID in a binary index file.
///
/// For all transcripts (including duplicates):
/// - Their original BED lines are written to a separate BED file.
///
/// # Arguments
///
/// * `transcripts` - A `Vec<GenePred>` containing the transcript records for the current chunk.
/// * `genome` - A reference to a `DashMap<String, Vec<u8>>` containing chromosome sequences.
/// * `chr` - A string slice representing the chromosome name for the current chunk.
/// * `path` - A `PathBuf` representing the base directory for the chunk's output files.
/// * `chunk_id` - A string slice representing the unique identifier for the current chunk.
/// * `writer_fa` - A mutable `BufWriter<File>` for writing FASTA entries for unique sequences.
/// * `writer_bed` - A mutable `BufWriter<File>` for writing all original BED entries.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to create the index file or the reduced BED file.
/// - It fails to write to any of the `BufWriter`s (`writer_fa`, `writer_bed`, `writer_reduced_bed`, `index`).
/// - `get_sequence` or `encode_id` functions panic due to invalid data or missing resources.
/// - It fails to split or join BED line fields.
fn index(
    mut transcripts: Vec<GenePred>,
    genome: &DashMap<String, Vec<u8>>,
    chr: &str,
    path: &PathBuf,
    chunk_id: &String,
    mut writer_fa: BufWriter<File>,
    mut writer_bed: BufWriter<File>,
    seq_mode: &SeqMode,
) {
    let mut mapper = HashMap::new();

    let idx = path.join(format!("tmp_chunk_{}.index", chunk_id));
    let mut index = BufWriter::new(File::create(&idx).unwrap());

    let chunk_reduced_bed = path.join(format!("tmp_chunk_{}_reduced.bed", &chunk_id));
    let mut writer_reduced_bed = BufWriter::new(File::create(&chunk_reduced_bed).unwrap());

    let mut count = 0usize;
    for tx in transcripts.iter_mut() {
        let seq = get_sequence(genome, chr, &tx, seq_mode);
        let key = seq.seq.as_bytes().to_vec();
        let encoded = encode_id(&tx.name);

        match mapper.entry(key) {
            Entry::Vacant(v) => {
                // INFO: first one for this seq
                v.insert(vec![count as u32, encoded]);

                // INFO: only for unseen seqs
                let mut fields: Vec<String> =
                    tx.line_mut().split('\t').map(|s| s.to_string()).collect();

                fields[3] = count.to_string();
                let line = fields.join("\t");

                writeln!(writer_reduced_bed, "{}", line).unwrap_or_else(|e| {
                    panic!("ERROR: could not write line from -> {:?}. {e}", tx)
                });
                writeln!(writer_fa, ">{}\n{}", count, seq).unwrap_or_else(|e| {
                    panic!(
                        "ERROR: could not write sequence to .fa from -> {:?}. {e}",
                        tx
                    )
                });

                log::debug!(
                    "NEW -> ENCODE: {encoded}, COUNT: {count}, NAME: {}",
                    &tx.name
                );

                count += 1;
            }
            Entry::Occupied(mut o) => {
                o.get_mut().push(encoded); // INFO: append to existing

                log::debug!(
                    "SEEN -> ENCODE: {encoded}, COUNT: {count}, NAME: {}",
                    &tx.name
                );
            }
        }

        writeln!(writer_bed, "{}", tx.line)
            .unwrap_or_else(|e| panic!("ERROR: could not write line from -> {:?}. {e}", tx));
    }

    // INFO: for every element in mapper -> write encoded id and encoded group
    for (_, group) in mapper {
        let header = &group[0];

        log::debug!("INSERTING: {header} as index for group: {group:?}");

        let n_ids = group.len() as u16 - 1;
        index.write_all(&n_ids.to_be_bytes()).unwrap();

        index.write_all(&header.to_be_bytes()).unwrap();

        for read in group.iter().skip(1) {
            index.write_all(&read.to_be_bytes()).unwrap();
        }
    }
}

/// Encodes a read ID string into a `u32` integer.
///
/// This function expects read IDs to be in a specific format, typically starting
/// with 'R' followed by a number, and then potentially more information
/// separated by underscores. It extracts the numeric part immediately
/// following 'R' and parses it as a `u32`.
///
/// # Arguments
///
/// * `id` - A reference to a `String` containing the read ID.
///          Expected format example: "R9834_chr16__FC37#TC0#PA0#PR0#IY887"
///
/// # Returns
///
/// A `u32` integer representing the numeric part of the read ID.
///
/// # Panics
///
/// This function will panic if:
/// - The `id` string does not contain an underscore (i.e., `split('_').next()` fails).
/// - The part before the first underscore does not start with 'R'
///   (i.e., `strip_prefix('R')` fails).
/// - The remaining string after stripping 'R' cannot be successfully parsed as a `u32` integer.
///
/// # Example
///
/// ```rust, ignore
/// let id1 = "R9834_chr16__FC37#TC0#PA0#PR0#IY887".to_string();
/// let encoded1 = encode_id(&id1);
/// assert_eq!(encoded1, 9834);
///
/// let id2 = "R123_abc".to_string();
/// let encoded2 = encode_id(&id2);
/// assert_eq!(encoded2, 123);
///
/// // Example of a string that would cause a panic:
/// // let bad_id1 = "9834_chr16".to_string();
/// // let _ = encode_id(&bad_id1); // Panics: "ERROR: could not preserve numeric part from 9834_chr16"
///
/// // let bad_id2 = "Rabc_chr16".to_string();
/// // let _ = encode_id(&bad_id2); // Panics: "ERROR: could not parse as number: Rabc_chr16 -> invalid digit found in string"
/// ```
fn encode_id(id: &String) -> u32 {
    // INFO: R9834_chr16__FC37#TC0#PA0#PR0#IY887
    id.split('_')
        .next()
        .unwrap_or_else(|| panic!("ERROR: could not get read number from {id}"))
        .strip_prefix('R')
        .unwrap_or_else(|| panic!("ERROR: could not preserve numeric part from {id}"))
        .parse::<u32>()
        .unwrap_or_else(|e| panic!("ERROR: could not parse as number: {id} -> {e}"))
}

/// Locates or extracts read IDs from a binary index file.
///
/// This function serves two primary purposes: finding a specific group of read IDs
/// associated with a given header, or iterating through the entire index and
/// writing all header-ID group pairs to an output file. The function reads a
/// custom binary format, which consists of a 2-byte count for the number of IDs
/// in a group, followed by a 4-byte header, and then a sequence of 4-byte IDs.
///
/// # Arguments
///
/// * `args` - A `crate::cli::IndexArgs` struct containing command-line arguments,
///            including the path to the index file (`args.index`), an optional
///            header ID to search for (`args.id`), and flags to control the
///            operation, such as whether to write the output (`args.write`)
///            and the output directory (`args.output_dir`).
///
/// # Panics
///
/// This function will panic if:
/// - It fails to open the specified index file.
/// - The `--write` flag is used, but the specified output directory cannot be created.
/// - The `--write` flag is used, but the output file cannot be created within that directory.
/// - A read operation from the index file fails unexpectedly before the end of the file is reached.
/// - The `--write` flag is not used and no header ID is provided via `--id`.
///
/// # Example
///
/// ```rust, ignore
/// // Example of finding a specific header ID:
/// let args = IndexArgs {
///     index: PathBuf::from("path/to/my/index.bin"),
///     id: Some(1234),
///     write: false,
///     output_dir: PathBuf::new(),
/// };
/// find(args); // This would print the group of IDs associated with header 1234
///
/// // Example of writing all header-ID groups to a file:
/// let args = IndexArgs {
///     index: PathBuf::from("path/to/my/index.bin"),
///     id: None,
///     write: true,
///     output_dir: PathBuf::from("output_dir"),
/// };
/// find(args); // This would create a file named 'index' in 'output_dir'
/// ```
pub fn find(args: crate::cli::IndexArgs) {
    use std::io::Read;

    let mut reader = std::io::BufReader::new(
        File::open(args.index).unwrap_or_else(|e| panic!("ERROR: failed to open index -> {e}")),
    );

    let mut writer = None;
    if args.write {
        std::fs::create_dir_all(&args.output_dir).unwrap();

        writer = Some(std::io::BufWriter::new(
            File::create(args.output_dir.join("index")).unwrap(),
        ));
    }

    loop {
        let mut id_count_buf = [0u8; 2];
        if reader.read_exact(&mut id_count_buf).is_err() {
            break; // EOF reached cleanly
        }
        let n_ids = u16::from_be_bytes(id_count_buf);

        let mut header_buf = [0u8; 4];
        reader.read_exact(&mut header_buf).unwrap();
        let header = u32::from_be_bytes(header_buf);

        if args.write {
            let mut id_buf = [0u8; 4];
            let mut group = Vec::with_capacity(n_ids as usize);
            for _ in 0..n_ids {
                reader.read_exact(&mut id_buf).unwrap();
                let id = u32::from_be_bytes(id_buf);

                // WARN: id fmt -> R{int}_chr{chr}, skipping tags!
                let name = format!("R{}", id);
                group.push(name);
            }

            let _ = writeln!(writer.as_mut().unwrap(), "{}\t{:?}", header, group);
        } else {
            let id = args.id.unwrap_or_else(|| {
                panic!("ERROR: you forgot to pass --id <ID>, otherwise use --write")
            });

            if header == id {
                let mut id_buf = [0u8; 4];
                let mut group = Vec::with_capacity(n_ids as usize);
                for _ in 0..n_ids {
                    reader.read_exact(&mut id_buf).unwrap();
                    let id = u32::from_be_bytes(id_buf);

                    // WARN: id fmt -> R{int}_chr{chr}, skipping tags!
                    let name = format!("R{}", id);
                    group.push(name);
                }

                let rs = format!("ID: {} -> {:?}", header, group);
                print!("{}", rs);
            }
        }
    }
}

/// Computes the reverse complement of a DNA sequence in-place.
///
/// This function takes a mutable reference to a vector of bytes representing a DNA sequence
/// and transforms it into its reverse complement. The sequence is reversed and each base is
/// complemented according to Watson-Crick base pairing rules: A ↔ T and C ↔ G. The operation
/// is performed in-place, modifying the original vector. Ambiguous bases (like 'N') and
/// unrecognized characters remain unchanged.
///
/// # Arguments
///
/// * `seq` - A mutable reference to a `Vec<u8>` containing the DNA sequence as ASCII bytes.
///           The sequence can contain uppercase or lowercase nucleotide characters (A, T, C, G).
///           Other characters (such as 'N' for ambiguous bases) are preserved unchanged.
///
/// # Panics
///
/// This function does not panic under normal circumstances.
///
/// # Example
///
/// ```rust, ignore
/// // Example with a simple DNA sequence:
/// let mut sequence = b"ATCG".to_vec();
/// __rev_complement_u8(&mut sequence);
/// assert_eq!(sequence, b"CGAT".to_vec());
///
/// // Example with mixed case:
/// let mut sequence = b"AtcG".to_vec();
/// __rev_complement_u8(&mut sequence);
/// assert_eq!(sequence, b"CGAT".to_vec());
///
/// // Example with ambiguous bases:
/// let mut sequence = b"ATCGN".to_vec();
/// __rev_complement_u8(&mut sequence);
/// assert_eq!(sequence, b"NCGAT".to_vec());
/// ```
///
/// # Notes
///
/// - The function handles both uppercase and lowercase nucleotide characters
/// - Ambiguous bases (N) and unrecognized characters remain in their original positions but reversed
/// - The algorithm is optimized to work in-place with O(1) additional memory usage
/// - Time complexity is O(n) where n is the length of the sequence
pub fn __rev_complement_u8(seq: &mut Vec<u8>) {
    // INFO: A <-> T, C <-> G
    let complement = |base: u8| match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        _ => base, // INFO: N or other ambiguous bases remain unchanged
    };

    let len = seq.len();
    for i in 0..(len / 2) {
        let j = len - 1 - i;
        let temp = complement(seq[i]);
        seq[i] = complement(seq[j]);
        seq[j] = temp;
    }

    // INFO: if the sequence length is odd, complement the middle base
    if len % 2 == 1 {
        seq[len / 2] = complement(seq[len / 2]);
    };
}
