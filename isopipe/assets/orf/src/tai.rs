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

use config::{BedColumn, BedColumnValue, SCALE, bed_to_struct_collection};
use dashmap::{DashMap, DashSet};
use hashbrown::HashMap;
use isopipe::config::depure;
use log::warn;
use packbed::{
    reader as bed_reader,
    record::{Bed6, GenePred},
};
use rayon::prelude::*;
use smol_str::ToSmolStr;

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::PathBuf;

use crate::{cli::TaiArgs, utils::*};

const PREDICTIONS: &str = "predictions";
const RESULT: &str = "tai.result";
pub const TAI_VENV: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tai/.venv/bin/activate");

const COLUMNS: &[BedColumn] = &[
    BedColumn::Chrom,
    BedColumn::Start,
    BedColumn::End,
    BedColumn::Name,
    BedColumn::Strand,
];

/// Runs Translation AI on a set of query reads
///
/// This function orchestrates the entire Tai analysis, which involves:
/// 1. Creating an output directory.
/// 2. Reformatting the input FASTA and BED files into a format suitable for Tai.
/// 3. Executing the `translationai` command-line tool.
/// 4. Reading the Tai index.
/// 5. Unrolling and processing the TAI predictions, merging them with original
///    genomic alignment information, and writing the final results to an output file.
///
/// # Arguments
///
/// * `args` - A `TaiArgs` struct containing paths to input files, output directory,
///            and TAI-specific parameters like the prediction threshold.
///
/// # Panics
///
/// This function will panic if:
/// - The output directory cannot be created.
/// - The `refmt` function fails to reformat the input files.
/// - The `translationai` command fails to execute.
/// - The TAI index cannot be read.
/// - The `unroll_tai` function encounters any errors during prediction processing.
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
///
/// struct CommonArgs {
///     fasta: PathBuf,
///     alignments: PathBuf,
///     outdir: PathBuf,
/// }
///
/// struct TaiArgs {
///     common: CommonArgs,
///     threshold: f32,
/// }
///
/// let args = TaiArgs {
///     common: CommonArgs {
///         fasta: PathBuf::from("/path/to/input.fasta"),
///         alignments: PathBuf::from("/path/to/alignments.bed"),
///         outdir: PathBuf::from("/path/to/output_dir"),
///     },
///     threshold: 0.5,
/// };
///
/// run_tai(args);
/// ```
pub fn run_tai(args: TaiArgs) {
    let mode = Mode::from(&args.common.index);

    let dir = args.common.outdir.join("tai");
    std::fs::create_dir_all(&dir)
        .unwrap_or_else(|e| panic!("ERROR: could not create directory -> {e}!"));

    let (fasta, index) = refmt(
        &args.common.fasta,
        &args.common.alignments,
        &dir,
        &mode,
        &args.common.index,
    );

    let cmd = format!(
        "source {} && translationai -I {} -t {},{} -O {}",
        TAI_VENV,
        fasta.display(),
        args.threshold,
        args.threshold,
        fasta.with_extension(PREDICTIONS).display()
    );

    std::process::Command::new("bash")
        .arg("-c")
        .arg(cmd)
        .status()
        .unwrap_or_else(|e| panic!("ERROR: failed to execute orfipy command -> {e}"));

    unroll_tai(index, fasta, &args.common.alignments, &mode);

    isopipe::depure!(dir, "result");
}

/// Unrolls and processes Translation AI (TAI) predictions, merging them with
/// original genomic alignment information.
///
/// This function reads the TAI prediction file (which contains predicted ORFs
/// and their scores), and an index that maps original sequence IDs to their
/// reformatted counterparts. It then iterates through the predictions,
/// retrieves corresponding genomic alignment records, maps the predicted
/// ORF coordinates to absolute genomic coordinates, and writes the inflated
/// results to a new output file.
///
/// # Arguments
///
/// * `index` - A `HashMap<String, Vec<String>>` mapping reformatted sequence IDs
///             (references) to a list of original query IDs that map to them.
/// * `fasta` - A `PathBuf` pointing to the reformatted FASTA file, which is
///             used to derive the path to the TAI predictions file.
/// * `alignments` - A `PathBuf` pointing to the original BED file containing
///                  genomic alignments, used to retrieve `GenePred` records.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to read the TAI predictions file or the original BED alignments file.
/// - It fails to create the output result file.
/// - It cannot parse parts of the prediction line (e.g., start, stop, scores).
/// - A predicted ORF's start position is greater than its stop position.
/// - A chromosome or canonical ID from the prediction is not found in the
///   `records` (gene prediction collection).
/// - It fails to write to the output file.
///
/// # Example
///
/// ```rust, no_run
/// use std::collections::HashMap;
/// use std::path::PathBuf;
/// use std::fs::File;
/// use std::io::Write;
/// use tempfile::NamedTempFile;
///
/// // Create dummy files and data for the example
/// let mut predictions_file = NamedTempFile::new().unwrap();
/// predictions_file.write_all(b">chr16:91343975-91360783 (+)(R9834_chr16__FC37#TC0#PA0#PR0#IY887)	13190,13667,0.641,0.921\n").unwrap();
/// let predictions_path = predictions_file.path().to_path_buf();
///
/// let mut alignments_file = NamedTempFile::new().unwrap();
/// alignments_file.write_all(b"chr16\t91343975\t91360783\tR9834_chr16__FC37#TC0#PA0#PR0#IY887\t0\t+\t91343975\t91360783\t0,0,0\t1\t16808,\t0,\n").unwrap();
/// let alignments_path = alignments_file.path().to_path_buf();
///
/// let mut index_map: HashMap<String, Vec<String>> = HashMap::new();
/// index_map.insert("R9834_chr16".to_string(), vec!["R9834_chr16__FC37#TC0#PA0#PR0#IY887".to_string()]);
///
/// // Assuming `PREDICTIONS` and `RESULT` constants are defined
/// const PREDICTIONS: &str = "predictions";
/// const RESULT: &str = "result";
///
/// // Rename the dummy predictions file to match the expected extension
/// let final_predictions_path = predictions_path.with_extension(PREDICTIONS);
/// std::fs::rename(&predictions_path, &final_predictions_path).unwrap();
///
/// unroll_tai(index_map, final_predictions_path.clone(), &alignments_path);
///
/// // Clean up dummy files
/// std::fs::remove_file(&final_predictions_path).unwrap();
/// std::fs::remove_file(&final_predictions_path.with_extension(RESULT)).unwrap();
/// std::fs::remove_file(&alignments_path).unwrap();
/// ```
fn unroll_tai(index: PathBuf, fasta: PathBuf, alignments: &PathBuf, mode: &Mode) {
    // INFO: predictions in fmt -> id st,sp,st_score,sp_score ...
    let predictions = bed_reader(fasta.with_extension(PREDICTIONS))
        .unwrap_or_else(|e| panic!("ERROR: failed to read predictions file -> {e}"));

    let output = fasta.with_extension(RESULT);
    let writer = BufWriter::new(
        File::create(&output)
            .unwrap_or_else(|e| panic!("ERROR: cannot create index from sequences -> {e}")),
    );

    let accumulator = DashSet::new();

    let records = bed_to_struct_collection::<GenePred>(
        bed_reader(alignments)
            .unwrap_or_else(|e| panic!("ERROR: failed to read BED file -> {e}"))
            .into(),
        config::BedColumn::Name,
    )
    .unwrap_or_else(|e| panic!("ERROR: failed construct BED to GenePred collection -> {e}"));

    match mode {
        Mode::Raw => cannonical(predictions, records, index, accumulator, writer),
        Mode::Indexed => indexed(predictions, records, index, accumulator, writer),
    }
}

fn cannonical(
    predictions: String,
    records: DashMap<String, HashMap<String, GenePred>>,
    index: PathBuf,
    accumulator: DashSet<String>,
    mut writer: BufWriter<File>,
) {
    let index = read_tai_index(&index);

    // INFO: inflate results!
    predictions
        .par_lines()
        .filter(|line| !line.starts_with('#'))
        .for_each(|line| {
            let parts: Vec<&str> = line.split('\t').collect();

            // INFO: >chr16:91343975-91360783 +) R9834_chr16__FC37#TC0#PA0#PR0#IY887) 0, 0,)
            // INFO: >chr16:91343975-91360783 +) 0) 0, 0,) -> if extract mode on!
            let name = parts[0].split("(").collect::<Vec<&str>>();

            let strand = name[1].trim_end_matches(')').to_string();

            // INFO: R9834_chr16__FC37#TC0#PA0#PR0#IY887
            // INFO: 0 -> if extract mode on!
            let cannonical_id = name[2].trim_end_matches(')').to_string();

            let id = cannonical_id.split(config::BIG_SEP).collect::<Vec<&str>>()[0]; // INFO: R9834_chr16
            let chr = id.split("_").collect::<Vec<&str>>()[1];

            // INFO: for each query all orfs in the current record!
            // WARN: leaving option because unique ids are not mapped to an index!
            let queries = index.get(id);

            for (orf_idx, orf) in parts.iter().skip(1).enumerate() {
                // INFO: 13190,13667,0.6413461491465569,0.921319767832756
                let orf_parts: Vec<&str> = orf.split(',').collect();
                if orf_parts.len() < 2 {
                    panic!("ERROR: ORF does not have enough parts to parse: {}", orf);
                }

                let start = orf_parts[0].parse::<u64>().unwrap_or_else(|_| {
                    panic!("ERROR: failed to parse start position from ORF: {}", orf);
                });
                let stop = orf_parts[1].parse::<u64>().unwrap_or_else(|_| {
                    panic!("ERROR: failed to parse stop position from ORF: {}", orf);
                }) + 3; // INFO: stop is inclusive, so we add 3 to include the stop codon
                let start_score = orf_parts[2].parse::<f32>().unwrap_or_else(|_| {
                    panic!("ERROR: failed to parse start position from ORF: {}", orf);
                });
                let stop_score = orf_parts[3].parse::<f32>().unwrap_or_else(|_| {
                    panic!("ERROR: failed to parse stop position from ORF: {}", orf);
                });

                if start > stop {
                    panic!(
                        "ERROR: start position is greater than stop position in ORF: {}",
                        orf
                    );
                }

                // INFO: retrieving the reference gene prediction record
                let (orf_start, orf_end) = records
                    .get_mut(chr)
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
                    .map_absolute_cds(start as u64, stop as u64)
                    .unwrap_or_default();

                // WARN: skipping unreliable ORFs for the current alignment
                if orf_start == 0 && orf_end == 0 {
                    warn!(
                        "WARN: ORF start and end are zero for ID: {}.p{}, skipping!",
                        id, orf_idx
                    );
                    continue;
                }

                // INFO: retrieving the reference gene prediction record
                // INFO: since indexing groups exact similar records
                // INFO: we safely assume ref gp record could be applied to all queries
                let ref_id = format!("{}.p{}", id, orf_idx + 1);
                let ref_line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    ref_id, start, stop, start_score, stop_score, strand, orf_start, orf_end,
                );

                accumulator.insert(ref_line);

                if let Some(queries) = queries {
                    // INFO: queries are the orfs in the current record
                    for query in queries {
                        let query_id = format!("{}.p{}", query, orf_idx + 1);
                        let query_line = format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            query_id,
                            start,
                            stop,
                            start_score,
                            stop_score,
                            strand,
                            orf_start,
                            orf_end
                        );

                        accumulator.insert(query_line);
                    }
                } else {
                    // here we should do: if log level is debug -> warn
                    warn!("WARN: no queries found for ID: {}", id);
                    continue;
                }
            }
        });

    accumulator.into_iter().for_each(|line| {
        writeln!(writer, "{}", line).unwrap();
    });
}

fn indexed(
    predictions: String,
    records: DashMap<String, HashMap<String, GenePred>>,
    index: PathBuf,
    accumulator: DashSet<String>,
    mut writer: BufWriter<File>,
) {
    let chr = get_chr_from_path(&index);
    let index = extract::read::read_index(&index, &chr);

    // INFO: inflate results!
    predictions
        .par_lines()
        .filter(|line| !line.starts_with('#'))
        .for_each(|line| {
            let parts: Vec<&str> = line.split('\t').collect();

            // INFO: >chr16:91343975-91360783 +) R9834_chr16__FC37#TC0#PA0#PR0#IY887) 0, 0,)
            // INFO: >chr16:91343975-91360783 +) 0) 0, 0,) -> if extract mode on!
            let name = parts[0].split("(").collect::<Vec<&str>>();

            let strand = name[1].trim_end_matches(')').to_string();

            // INFO: R9834_chr16__FC37#TC0#PA0#PR0#IY887
            // INFO: 0 -> if extract mode on!
            let cannonical_id = name[2].trim_end_matches(')').to_string();
            let _id = cannonical_id
                .parse::<u32>()
                .unwrap_or_else(|e| panic!("ERROR: could parse ID to u32 -> {cannonical_id}. {e}"));

            // INFO: unpacking index reference -> queries
            let queries = index.get(&_id).unwrap_or_else(|| {
                panic!("ERROR: could not find index match for id {cannonical_id}")
            });

            for (orf_idx, orf) in parts.iter().skip(1).enumerate() {
                // INFO: 13190,13667,0.6413461491465569,0.921319767832756
                let orf_parts: Vec<&str> = orf.split(',').collect();
                if orf_parts.len() < 2 {
                    panic!("ERROR: ORF does not have enough parts to parse: {}", orf);
                }

                let start = orf_parts[0].parse::<u64>().unwrap_or_else(|_| {
                    panic!("ERROR: failed to parse start position from ORF: {}", orf);
                });
                let stop = orf_parts[1].parse::<u64>().unwrap_or_else(|_| {
                    panic!("ERROR: failed to parse stop position from ORF: {}", orf);
                }) + 3; // INFO: stop is inclusive, so we add 3 to include the stop codon
                let start_score = orf_parts[2].parse::<f32>().unwrap_or_else(|_| {
                    panic!("ERROR: failed to parse start position from ORF: {}", orf);
                });
                let stop_score = orf_parts[3].parse::<f32>().unwrap_or_else(|_| {
                    panic!("ERROR: failed to parse stop position from ORF: {}", orf);
                });

                if start > stop {
                    panic!(
                        "ERROR: start position is greater than stop position in ORF: {}",
                        orf
                    );
                }

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
                    .map_absolute_cds(start as u64, stop as u64)
                    .unwrap_or_default();

                // WARN: skipping unreliable ORFs for the current alignment
                if orf_start == 0 && orf_end == 0 {
                    // warn!(
                    //     "WARN: ORF start and end are zero for ID: {}.p{}, skipping!",
                    //     _id, orf_idx
                    // );
                    continue;
                }

                // INFO: queries are the orfs in the current record
                for query in queries {
                    let query_id = format!("{}.p{}", query, orf_idx + 1);
                    let query_line = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        query_id, start, stop, start_score, stop_score, strand, orf_start, orf_end
                    );

                    accumulator.insert(query_line);
                }
            }
        });

    accumulator.into_iter().for_each(|line| {
        writeln!(writer, "{}", line).unwrap();
    });
}

/// Reformats FASTA headers based on corresponding BED alignment information
/// and creates a new FASTA file and an associated index.
///
/// This function reads both a FASTA file and a BED file. It then iterates
/// through the BED records, uses the information to reformat the corresponding
/// FASTA headers, and writes the reformatted sequences to a new FASTA file.
/// Additionally, it generates a binary index file that maps the original
/// sequence IDs (from the reformatted headers) to their deduplicated counterparts
/// and their chromosome.
///
/// The reformatted header format is: `>chr:start-end(strand)(original_name)(0, 0,)`
///
/// # Arguments
///
/// * `fasta` - A `PathBuf` representing the path to the input FASTA file.
/// * `bed` - A `PathBuf` representing the path to the input BED file containing
///           alignment information.
/// * `outdir` - A `PathBuf` representing the output directory where the
///              reformatted FASTA and index files will be created.
///
/// # Returns
///
/// A tuple containing two `PathBuf` instances:
/// - The path to the reformatted FASTA file.
/// - The path to the generated index file.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to parse the BED file or the FASTA file.
/// - It cannot create the reformatted FASTA file or the index file.
/// - A sequence from the BED file's header is not found in the parsed FASTA sequences.
/// - It fails to encode a header for the TAI index.
/// - It fails to write to the reformatted FASTA file or the index file.
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
/// use std::fs::File;
/// use std::io::Write;
/// use tempfile::NamedTempFile;
/// use crate::config::{self, BedColumn, BedColumnValue, BedParser};
///
/// // Assume Bed6, parse_bed, parse_fa, create_fasta, create_index,
/// // get_header_from_values, and encode_for_tai are available.
/// // Also assume COLUMNS constant is defined.
/// const COLUMNS: &[BedColumn] = &[
///     BedColumn::Chrom, BedColumn::Start, BedColumn::End, BedColumn::Name, BedColumn::Strand
/// ];
///
/// // Create dummy FASTA and BED files
/// let mut fasta_file = NamedTempFile::new().unwrap();
/// fasta_file.write_all(b">R1_chr1__FC#TC#PA#PR#IY_ORF.1 [10-20](+)_type:complete\nATGC\n").unwrap();
/// let fasta_path = fasta_file.path().to_path_buf();
///
/// let mut bed_file = NamedTempFile::new().unwrap();
/// bed_file.write_all(b"chr1\t10\t20\tR1_chr1__FC#TC#PA#PR#IY\t0\t+\n").unwrap();
/// let bed_path = bed_file.path().to_path_buf();
///
/// let outdir = PathBuf::from("temp_tai_out");
/// std::fs::create_dir_all(&outdir).unwrap();
///
/// let (reformatted_fasta, tai_index) = refmt(&fasta_path, &bed_path, &outdir);
///
/// assert!(reformatted_fasta.exists());
/// assert!(tai_index.exists());
///
/// // Clean up dummy files and directory
/// std::fs::remove_file(&reformatted_fasta).unwrap();
/// std::fs::remove_file(&tai_index).unwrap();
/// std::fs::remove_dir_all(&outdir).unwrap();
/// ```
fn refmt(
    fasta: &PathBuf,
    bed: &PathBuf,
    outdir: &PathBuf,
    mode: &Mode,
    index_path: &Option<PathBuf>,
) -> (PathBuf, PathBuf) {
    let records = parse_bed::<Bed6>(bed, COLUMNS.to_vec());
    let seqs =
        parse_fa(fasta).unwrap_or_else(|e| panic!("ERROR: failed to parse FASTA file -> {e}"));

    let fmt = outdir.join(format!(
        "{}.fmt.fa",
        fasta.file_stem().unwrap().to_str().unwrap()
    ));
    let mut writer = create_fasta(&fmt, "fa")
        .unwrap_or_else(|| panic!("ERROR: could not create file {:?}", fmt));
    let mut index = create_index(&fmt);

    let accumulator = DashMap::new();

    records.into_par_iter().for_each(|(chr, rows)| {
        for (header, values) in rows.iter() {
            let fmt_header = get_header_from_values(values, header);
            let seq = seqs.get(&header.clone().to_smolstr()).unwrap_or_else(|| {
                panic!("ERROR: chromosome from {header} not found in sequences -> {chr}!");
            });

            accumulator
                .entry(seq)
                .or_insert_with(Vec::new)
                .push(fmt_header)
        }
    });

    // INFO: will write the first header on collection and will point other headers to it
    accumulator.into_iter().for_each(|(seq, headers)| {
        writeln!(writer, "{}", &headers[0]).unwrap();
        writer.write_all(seq).unwrap();
        writeln!(writer).unwrap();

        match mode {
            // INFO: nothing bc we already have an index!
            Mode::Indexed => {}
            Mode::Raw => {
                if headers.len() > 1 {
                    let mut count = 0;
                    for header in &headers {
                        let (id, chr) = if let Some((read, chr)) = encode_for_tai(header) {
                            (read, chr)
                        } else {
                            panic!("ERROR: failed to encode header: {}", header);
                        };

                        // INFO: first record -> chr byte len + chr bytes
                        // INFO: chr-len chr n_ids id[ref] id1 id2
                        if count < 1 {
                            index.write_all(&[chr.len() as u8]).unwrap();
                            index.write_all(chr.as_bytes()).unwrap();

                            let n_ids = headers.len() as u16;
                            index.write_all(&n_ids.to_be_bytes()).unwrap();
                        }

                        index.write_all(&id.to_be_bytes()).unwrap();

                        count += 1;
                    }
                }
            }
        }
    });

    match mode {
        Mode::Raw => (fmt.clone(), fmt.with_extension("dedup.index")),
        Mode::Indexed => (fmt.clone(), index_path.as_ref().unwrap().clone()),
    }
}

/// Reads a TAI index file into a `HashMap` that maps reference IDs to a list of query IDs.
///
/// The TAI index file is a custom binary format used to store mappings between
/// a "reference" sequence ID (typically the first one in a deduplicated group)
/// and all other "query" sequence IDs that are identical to it.
///
/// The binary format for each entry is:
/// - `chr_len`: `u8` (length of the chromosome name in bytes)
/// - `chr_bytes`: `[u8; chr_len]` (chromosome name as bytes)
/// - `n_ids`: `u16` (total number of IDs in this group, including the reference)
/// - Followed by `n_ids` number of `u16` IDs:
///     - `id_reference`: `u16` (the ID of the reference sequence)
///     - `id_1`: `u16` (ID of the first query sequence)
///     - ...
///     - `id_n`: `u16` (ID of the last query sequence)
///
/// The function reconstructs `R{id}_chr{chr}` strings for both reference and query IDs.
///
/// # Arguments
///
/// * `index` - A `PathBuf` representing the path to the TAI index file.
///
/// # Returns
///
/// A `HashMap<String, Vec<String>>` where:
/// - Keys are the reference IDs (e.g., "R123_chrX").
/// - Values are `Vec<String>` containing the query IDs (e.g., ["R456_chrX", "R789_chrX"])
///   that map to the reference ID. Note that the reference itself is not included in the `Vec`.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to open the index file.
/// - It encounters an `io::Error` during reading, unless it's an EOF.
/// - It fails to convert chromosome bytes to a UTF-8 string.
/// - It fails to convert byte arrays to `u16` integers.
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
/// use std::fs::File;
/// use std::io::Write;
/// use tempfile::NamedTempFile;
/// use std::collections::HashMap;
///
/// // Create a dummy TAI index file
/// let mut temp_file = NamedTempFile::new().unwrap();
/// // Entry 1: Reference R1_chr1, Query R2_chr1
/// temp_file.write_all(&1u8.to_be_bytes()).unwrap(); // chr_len = 1 ('1')
/// temp_file.write_all(b"1").unwrap();                // chr_bytes
/// temp_file.write_all(&2u16.to_be_bytes()).unwrap(); // n_ids = 2 (ref + 1 query)
/// temp_file.write_all(&1u16.to_be_bytes()).unwrap(); // id_reference = 1
/// temp_file.write_all(&2u16.to_be_bytes()).unwrap(); // id_1 = 2
///
/// // Entry 2: Reference R3_chr2
/// temp_file.write_all(&1u8.to_be_bytes()).unwrap(); // chr_len = 1 ('2')
/// temp_file.write_all(b"2").unwrap();                // chr_bytes
/// temp_file.write_all(&1u16.to_be_bytes()).unwrap(); // n_ids = 1 (ref only)
/// temp_file.write_all(&3u16.to_be_bytes()).unwrap(); // id_reference = 3
///
/// let index_path = temp_file.path().to_path_buf();
///
/// let index_map = read_tai_index(&index_path);
///
/// assert!(index_map.contains_key("R1_chr1"));
/// assert_eq!(index_map["R1_chr1"], vec!["R2_chr1".to_string()]);
///
/// assert!(index_map.contains_key("R3_chr2"));
/// assert!(index_map["R3_chr2"].is_empty()); // No queries for R3_chr2
///
/// // Clean up dummy file
/// std::fs::remove_file(&index_path).unwrap();
/// ```
fn read_tai_index(index: &PathBuf) -> HashMap<String, Vec<String>> {
    let mut mapper = HashMap::new();

    let mut reader = BufReader::new(
        File::open(index).unwrap_or_else(|e| panic!("ERROR: failed to open index -> {e}")),
    );

    loop {
        let mut chr_len_buf = [0u8; 1];
        if reader.read_exact(&mut chr_len_buf).is_err() {
            break; // EOF reached cleanly
        }
        let chr_len = chr_len_buf[0] as usize;

        let mut chr_buf = vec![0u8; chr_len];
        reader.read_exact(&mut chr_buf).unwrap();
        let chr = String::from_utf8(chr_buf)
            .unwrap_or_else(|e| panic!("ERROR: failed to convert chr bytes to string -> {e}"));

        let mut id_count_buf = [0u8; 2];
        reader.read_exact(&mut id_count_buf).unwrap();
        let n_ids = u16::from_be_bytes(id_count_buf);

        let mut id_buf = [0u8; 2];
        let mut group = Vec::with_capacity(n_ids as usize);
        for _ in 0..n_ids {
            reader.read_exact(&mut id_buf).unwrap();
            let id = u16::from_be_bytes(id_buf);

            // WARN: id fmt -> R{int}_chr{chr}, skipping tags!
            let name = format!("R{}_chr{}", id, chr);
            group.push(name);
        }

        let reference = group[0].clone();
        let queries = group[1..].to_vec();

        mapper
            .entry(reference)
            .or_insert_with(Vec::new)
            .extend(queries);
    }

    mapper
}

/// Encodes a reformatted FASTA header into a tuple containing the read ID and chromosome string.
///
/// This function parses a header string that is expected to follow a specific format
/// generated by `get_header_from_values` and extracts the numeric read ID and the
/// chromosome string.
///
/// Expected header format: `>chr:start-end(strand)(R<read_id>_chr<chr_num>__...)`
///
/// # Arguments
///
/// * `header` - A `String` reference to the reformatted FASTA header.
///
/// # Returns
///
/// An `Option<(u16, &str)>`:
/// - `Some((read_id, chr_str))`: If the header is successfully parsed,
///   containing the read ID as `u16` and the chromosome as a string slice.
/// - `None`: If the header does not match the expected format or if parsing
///   of the read ID fails.
///
/// # Panics
///
/// This function will panic if:
/// - The header string does not have enough parts after splitting by `'('`.
/// - The "R" prefix is missing from the read ID part.
/// - The "chr" prefix is missing from the chromosome part.
/// - Parsing of the read ID into a `u16` fails.
///
/// # Example
///
/// ```rust, no_run
/// let header = ">chr16:45612921-45619040(+)(R6713_chr16__FC48#TC40#PA0#PR0#IY876)(0, 0,)".to_string();
/// let (read_id, chr) = encode_for_tai(&header).unwrap();
/// assert_eq!(read_id, 6713);
/// assert_eq!(chr, "16");
///
/// let malformed_header = ">malformed_header".to_string();
/// assert!(encode_for_tai(&malformed_header).is_none());
/// ```
fn encode_for_tai(header: &String) -> Option<(u16, &str)> {
    // INFO: ">chr16:45612921-45619040(+)(R6713_chr16__FC48#TC40#PA0#PR0#IY876)(0, 0,)",
    // INFO: 6713 16
    let parts: Vec<&str> = header.split('(').collect();
    if parts.len() < 2 {
        panic!(
            "ERROR: header does not have enough parts to parse: {}",
            header
        );
    }

    let id = parts[2].trim_end_matches(')');
    let id_parts: Vec<&str> = id.split('_').collect();

    let read = id_parts[0]
        .strip_prefix('R')?
        .parse::<u16>()
        .ok()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse read ID from header: {}", header);
        });
    let chr = id_parts[1].strip_prefix("chr").unwrap_or_else(|| {
        panic!("ERROR: failed to parse chromosome from header: {}", header);
    });

    Some((read, chr))
}

/// Reformats a FASTA header string based on extracted `BedColumnValue`s.
///
/// This function takes a vector of parsed BED column values and the original
/// header string, then constructs a new header string following a specific format:
/// `>chr:start-end(strand)(name)(0, 0,)`.
/// It also handles strand-specific coordinate transformations if the strand is '-'.
///
/// # Arguments
///
/// * `values` - A reference to a `Vec<BedColumnValue>` containing the parsed
///              values for chromosome, start, end, name, and strand.
/// * `header` - A reference to the original header `String`. Used primarily for error messages.
///
/// # Returns
///
/// A `String` representing the newly formatted header.
///
/// # Panics
///
/// This function will panic if:
/// - Any of the required `BedColumnValue`s (chromosome, start, end, name, strand)
///   cannot be extracted or are of an unexpected type.
/// - The `strand` value is not `"+"` or `"-"`.
/// - `SCALE` constant is not defined or accessible for strand transformation.
///
/// # Example
///
/// ```rust, no_run
/// use crate::config::{BedColumnValue, BedColumn};
///
/// // Assume SCALE constant is defined, e.g., const SCALE: u64 = 1_000_000_000;
/// const SCALE: u64 = 1000;
///
/// let values_plus_strand = vec![
///     BedColumnValue::Str("chr1".to_string()),
///     BedColumnValue::Number(100),
///     BedColumnValue::Number(200),
///     BedColumnValue::Str("geneA".to_string()),
///     BedColumnValue::Str("+".to_string()),
/// ];
/// let header_plus = "original_header_plus".to_string();
/// let formatted_plus = get_header_from_values(&values_plus_strand, &header_plus);
/// assert_eq!(formatted_plus, ">chr1:100-200(+)(geneA)(0, 0,)");
///
/// let values_minus_strand = vec![
///     BedColumnValue::Str("chr2".to_string()),
///     BedColumnValue::Number(100),
///     BedColumnValue::Number(200),
///     BedColumnValue::Str("geneB".to_string()),
///     BedColumnValue::Str("-".to_string()),
/// ];
/// let header_minus = "original_header_minus".to_string();
/// let formatted_minus = get_header_from_values(&values_minus_strand, &header_minus);
/// // For '-' strand, start and end are transformed: start = SCALE - original_end, end = SCALE - original_start
/// assert_eq!(formatted_minus, ">chr2:800-900(-)(geneB)(0, 0,)");
/// ```
fn get_header_from_values(values: &Vec<BedColumnValue>, header: &String) -> String {
    let chr = values[0].as_str().unwrap_or_else(|| {
        panic!(
            "ERROR: start position not found for header: {} - {:?}",
            header, values
        )
    });
    let mut start = values[1].as_number().unwrap_or_else(|| {
        panic!(
            "ERROR: start position not found for header: {} - {:?}",
            header, values
        )
    });
    let mut end = values[2].as_number().unwrap_or_else(|| {
        panic!(
            "ERROR: end position not found for header: {} - {:?}",
            header, values
        )
    });
    let name = values[3].as_str().unwrap_or_else(|| {
        panic!(
            "ERROR: name not found for header: {} - {:?}",
            header, values
        )
    });
    let strand = values[4].as_str().unwrap_or_else(|| {
        panic!(
            "ERROR: strand not found for header: {} - {:?}",
            header, values
        )
    });

    match strand {
        "+" => {}
        "-" => {
            let tmp = start;
            start = SCALE - end;
            end = SCALE - tmp;
        }
        _ => {
            panic!(
                "ERROR: strand not recognized for header: {} - {:?}",
                header, values
            );
        }
    }

    let hdr = format!(">{}:{}-{}({})({})(0, 0,)", chr, start, end, strand, name);
    hdr
}
