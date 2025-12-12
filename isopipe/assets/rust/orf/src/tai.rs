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

use config::{bed_to_custom_struct_collection, BedColumnValue, Strand, SCALE};
use dashmap::{DashMap, DashSet};
use hashbrown::HashMap;
use isopipe::config::depure;
use log::warn;
use memchr::memchr;
use memmap2::Mmap;
use packbed::{reader as bed_reader, record::GenePred};
use rayon::prelude::*;
use smol_str::{SmolStr, ToSmolStr};

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::str::from_utf8;

use crate::{cli::TaiArgs, utils::*};

const PREDICTIONS: &str = "predictions";
pub const TAI_VENV: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tai/.venv/bin/activate");

pub const START_CODON: &str = "ATG";
pub const STOP_CODONS: [&str; 3] = ["TAA", "TAG", "TGA"];
pub const STOP_CODONS_BYTES: [&[u8]; 3] = [b"TAA", b"TAG", b"TGA"];

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
#[allow(unused_assignments)]
pub fn run_tai(args: TaiArgs) {
    let dir = args.common.outdir.join("tai");
    std::fs::create_dir_all(&dir)
        .unwrap_or_else(|e| panic!("ERROR: could not create directory -> {e}!"));

    let chr = get_chr_from_filename(&args.common.index);
    let index = extract::read::read_index(&args.common.index, &chr);

    let records = flatten_dashmap(
        bed_to_custom_struct_collection::<GenePred>(
            bed_reader(&args.common.alignments)
                .unwrap_or_else(|e| panic!("ERROR: failed to read BED file -> {e}"))
                .into(),
            config::BedColumn::Name,
            config::BedOperation::SplitName("__", 0),
        )
        .unwrap_or_else(|e| panic!("ERROR: failed construct BED to GenePred collection -> {e}")),
    );

    let (fasta, sequences) = refmt(&args.common.fasta, &records, &dir, &index);

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

    let predictions = bed_reader(fasta.with_extension(PREDICTIONS))
        .unwrap_or_else(|e| panic!("ERROR: failed to read predictions file -> {e}"));

    let writer = BufWriter::new(
        File::create(fasta.with_extension("result"))
            .unwrap_or_else(|e| panic!("ERROR: cannot create index from sequences -> {e}")),
    );

    unroll_tai(index, predictions, records, writer, sequences);
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
#[allow(unused_assignments)]
fn unroll_tai(
    index: std::collections::HashMap<u32, Vec<String>>,
    predictions: String,
    mut records: HashMap<String, GenePred>,
    mut writer: BufWriter<File>,
    sequences: HashMap<smol_str::SmolStr, Vec<u8>>,
) {
    let accumulator = DashSet::new();

    for line in predictions.lines() {
        let parts: Vec<&str> = line.split('\t').collect();

        // INFO: >chr6:128352418-128362107(+)(ENSMUST00000100926.4)(0, 0,)
        // INFO: >chr16:91343975-91360783 +) R9834_chr16__FC37#TC0#PA0#PR0#IY887) 0, 0,)
        // INFO: >chr16:91343975-91360783 +) 0) 0, 0,) -> if extract mode on!
        let mut id = parts
            .first()
            .unwrap_or_else(|| {
                panic!("ERROR: ID not found in line: {}", line);
            })
            .to_string();

        // INFO: ENSMUST00000100926.4
        // INFO: R9834_chr16__FC37#TC0#PA0#PR0#IY887
        // INFO: 0 -> if extract mode on!
        id = id.split("(").collect::<Vec<&str>>()[2]
            .split(")")
            .collect::<Vec<&str>>()[0]
            .to_string();

        let u_id = id
            .parse::<u32>()
            .unwrap_or_else(|e| panic!("ERROR: could not parse ID to u32 -> {id}. {e}"));

        // INFO: unpacking index reference -> queries
        let queries = index
            .get(&u_id)
            .unwrap_or_else(|| panic!("ERROR: could not find index match for id {id}"));

        for (orf_idx, orf) in parts.iter().skip(1).enumerate() {
            // INFO: 13190,13667,0.6413461491465569,0.921319767832756
            let orf_parts: Vec<&str> = orf.split(',').collect();
            if orf_parts.len() < 2 {
                panic!("ERROR: ORF does not have enough parts to parse: {}", orf);
            }

            let start = orf_parts[0].parse::<u64>().unwrap_or_else(|_| {
                panic!("ERROR: failed to parse start position from ORF: {}", orf);
            });
            let mut stop = orf_parts[1].parse::<u64>().unwrap_or_else(|_| {
                panic!("ERROR: failed to parse stop position from ORF: {}", orf);
            }); // INFO: stop is inclusive, so we add 3 to include the stop codon
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

            // queries -> { 1 : [ R1_chr1, R2_chr1 ]}
            for query in queries {
                // INFO: retrieving the reference gene prediction record from original bed
                // INFO: this means bed has the original name as a name (likely splitted by BedOperation)
                let record = records.get_mut(query).unwrap_or_else(|| {
                    panic!("ERROR: id not found in BED, this is a bug -> {}!", id);
                });

                let strand = record.strand.clone();

                let (mut orf_start, mut orf_end) = record.map_absolute_cds(start, stop);

                // WARN: skipping unreliable ORFs for the current alignment
                if orf_start == 0 && orf_end == 0 {
                    continue;
                }

                let mut stop_codon = "";

                // INFO: before adding up we need to check the stop codon to see if its a real one or not
                let sequence = sequences.get(&id.to_smolstr()).unwrap_or_else(|| {
                    panic!("ERROR: sequence not found for {id:?}");
                });
                let mut orf_sequence = sequence.clone();

                let start_codon =
                    from_utf8(&sequence[start as usize..(start + 3) as usize])
                        .unwrap_or_else(|e| {
                            panic!("ERROR: failed to parse start codon -> {e} -> {sequence:?} from {id:?} using {orf_start:?}");
                        });

                __check_start_codon(start_codon, &id, start);

                // INFO: stop is inclusive, so we add 3 to include the stop codon
                match strand.as_str() {
                    // WARN: some weird cases where the tool predicts a non-stopped ORF:
                    // WARN: R146001_manual_scaffold_1.p1    102     2199
                    // WARN: sizes -> 144,117,332,112,120,117,138,78,66,270,246,137,79,156,87 = 2199
                    // WARN: record -> manual_scaffold_1 189532046 189543938 R146001 60 + 189532148 189543941
                    // WARN: would be out-of-bounds if we add +3 in this case
                    "+" => {
                        if orf_end + 3 > record.end {
                            dbg!(
                                "WARN: translationAi predicted a non-stop ORF: {orf:?} for {gp:?}"
                            );
                            // orf_end = record.end
                            stop_codon =
                                from_utf8(&sequence[(stop - 3) as usize..(stop) as usize])
                                    .unwrap_or_else(|e| {
                                        panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {id:?} using {orf_end:?}");
                                    });
                            dbg!("WARN: non-stop ORF stop_codon picked -> {:?}", stop_codon);

                            orf_sequence = sequence[start as usize..stop as usize].to_vec();
                        } else {
                            stop_codon =
                                from_utf8(&sequence[(stop) as usize..(stop + 3) as usize])
                                    .unwrap_or_else(|e| {
                                        panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {id:?} using {orf_end:?}");
                                    });

                            if !STOP_CODONS.contains(&stop_codon) {
                                // WARN: if stop_codon is not cannonical this is probably a case where the tool is wrong
                                dbg!(
                                    "WARN: stop codon is not TAA, TAG, or TGA -> {:?} from {:?} using {:?}",
                                    stop_codon,
                                    &id,
                                    stop
                                );

                                // INFO: taking stop as last base and going back 2 nt to capture last codon
                                stop_codon =
                                    from_utf8(&sequence[(stop - 3) as usize..(stop) as usize])
                                        .unwrap_or_else(|e| {
                                            panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {id:?} using {orf_end:?}");
                                        });

                                dbg!("WARN: non-stop ORF stop_codon picked -> {:?}", stop_codon);
                                orf_sequence = sequence[start as usize..stop as usize].to_vec();
                            } else {
                                // INFO: stop_codon is cannonical so we can safely add 3 to the end
                                orf_end += 3;
                                orf_sequence =
                                    sequence[start as usize..(stop + 3) as usize].to_vec();
                                stop += 3;
                            }
                        }
                    }
                    "-" => {
                        if orf_start - 3 < config::SCALE - record.end {
                            dbg!(
                                "WARN: translationAi predicted a non-stop ORF: {orf:?} for {gp:?}"
                            );
                            // orf_start = config::SCALE - record.end
                            stop_codon =
                                from_utf8(&sequence[(stop - 3) as usize..(stop) as usize])
                                    .unwrap_or_else(|e| {
                                        panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {id:?} using {orf_end:?}");
                                    });
                            dbg!("WARN: non-stop ORF stop_codon picked -> {:?}", stop_codon);
                            orf_sequence = sequence[start as usize..stop as usize].to_vec();
                        } else {
                            stop_codon =
                                from_utf8(&sequence[(stop) as usize..(stop + 3) as usize])
                                    .unwrap_or_else(|e| {
                                        panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {id:?} using {orf_end:?}");
                                    });

                            if !STOP_CODONS.contains(&stop_codon) {
                                // WARN: if stop_codon is not cannonical this is probably a case where the tool is wrong
                                dbg!(
                                    "WARN: stop codon is not TAA, TAG, or TGA -> {:?} from {:?} using {:?}",
                                    stop_codon,
                                    &id,
                                    stop
                                );

                                // INFO: taking stop as last base and going back 2 nt to capture last codon
                                stop_codon =
                                    from_utf8(&sequence[(stop - 3) as usize..(stop) as usize])
                                        .unwrap_or_else(|e| {
                                            panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {id:?} using {orf_end:?}");
                                        });

                                dbg!("WARN: non-stop ORF stop_codon picked -> {:?}", stop_codon);
                                orf_sequence = sequence[start as usize..stop as usize].to_vec();
                            } else {
                                // INFO: stop_codon is cannonical so we can safely add 3 to the end
                                orf_start -= 3;
                                orf_sequence =
                                    sequence[start as usize..(stop + 3) as usize].to_vec();
                                stop += 3;
                            }
                        }
                    }
                    _ => panic!("ERROR: unexpected strand value: {}", strand),
                }

                let pep = translate(&orf_sequence);
                let inner_stops = scan_stops(orf_sequence);

                // INFO: queries are the orfs in the current record
                let query_id = format!("{}.p{}", query, orf_idx + 1);
                let query_line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    record.chrom,
                    orf_start,
                    orf_end,
                    query_id,
                    strand,
                    start_score,
                    stop_score,
                    start,
                    stop,
                    start_codon,
                    stop_codon,
                    inner_stops,
                    pep
                );

                accumulator.insert(query_line);
            }
        }
    }

    accumulator.into_iter().for_each(|line| {
        writeln!(writer, "{}", line).unwrap();
    });
}

/// Flattens a `DashMap` into a `HashMap`.
///
/// # Arguments
///
/// * `original` - The `DashMap` to flatten.
///
/// # Returns
///
/// A `HashMap` containing the flattened `DashMap`.
///
/// # Example
///
/// ```rust
/// use std::path::PathBuf;
/// use std::fs::File;
/// use std::io::BufWriter;
///
/// let original = DashMap::new();
/// let mut inner_map = HashMap::new();
/// inner_map.insert("key".to_string(), GenePred::default());
/// original.insert("id".to_string(), inner_map);
///
/// let result = flatten_dashmap(original);
/// ```
fn flatten_dashmap(
    original: DashMap<String, HashMap<String, GenePred>>,
) -> HashMap<String, GenePred> {
    let mut result = HashMap::new();
    let read_only = original;

    for (_, inner_map) in read_only {
        for (key, value) in inner_map {
            result.insert(key, value);
        }
    }

    result
}

/// Processes TAI (Translation AI) predictions in a "canonical" mode, which is used
/// when input sequences were not indexed during a preceding extraction step.
///
/// This function reads a TAI-specific index to map canonical IDs to their
/// original queries. It then parses the prediction file in parallel, translates
/// ORF coordinates to absolute CDS coordinates using a `GenePred` reference, and
/// "inflates" each prediction by generating output records for both the canonical ID
/// and all associated original queries. The results are collected in a thread-safe
/// accumulator and written to a file.
///
/// # Arguments
///
/// * `predictions` - A `String` containing the full content of the TAI prediction output.
/// * `records` - A `DashMap<String, HashMap<String, GenePred>>` holding the reference
///   gene prediction records, keyed by chromosome and then by ID.
/// * `index` - A `PathBuf` representing the path to the TAI index file.
/// * `accumulator` - A `DashSet<String>` used as a thread-safe set to collect the
///   output lines before writing to the final file.
/// * `writer` - A `BufWriter<File>` used for writing the final output to a file.
///
/// # Returns
///
/// This function does not return a value. It writes its output directly to the provided `writer`.
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
/// use std::collections::HashMap;
/// use std::io::BufWriter;
/// use std::fs::File;
/// use dashmap::DashMap;
/// use dashset::DashSet;
///
/// // Dummy types for the example
/// struct GenePred;
///
/// // Assume a function `cannonical` exists with the correct signature.
/// // cannonical(predictions: String, records: DashMap<String, HashMap<String, GenePred>>, ...);
///
/// let predictions = "example predictions...".to_string();
/// let records = DashMap::new();
/// let index_path = PathBuf::from("path/to/index.txt");
/// let accumulator = DashSet::new();
/// let writer = BufWriter::new(File::create("output.txt").unwrap());
///
/// // Calling the function would look like this:
/// // cannonical(predictions, records, index_path, accumulator, writer);
/// ```
#[deprecated(note = "Functionality replaced by extract indexing")]
fn _cannonical(
    predictions: String,
    records: DashMap<String, HashMap<String, GenePred>>,
    index: PathBuf,
    accumulator: DashSet<String>,
    mut writer: BufWriter<File>,
) {
    let index = _read_tai_index(&index);

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
                }); // INFO: stop is inclusive, so we add 3 to include the stop codon
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
                let mut chr_records = records.get_mut(chr).unwrap_or_else(|| {
                    panic!(
                        "ERROR: chromosome from {} not found in sequences -> {}!",
                        cannonical_id, chr
                    );
                });

                let gp = chr_records.get_mut(&cannonical_id).unwrap_or_else(|| {
                    panic!(
                        "ERROR: id not found in BED, this is a bug -> {}!",
                        cannonical_id
                    );
                });

                let (mut orf_start, mut orf_end) = gp.map_absolute_cds(start, stop);

                // WARN: skipping unreliable ORFs for the current alignment
                if orf_start == 0 && orf_end == 0 {
                    log::debug!(
                        "WARN: ORF start and end are zero for ID: {}.p{}, skipping!",
                        id,
                        orf_idx
                    );
                    continue;
                }

                // INFO: stop is inclusive, so we add 3 to include the stop codon
                match strand.as_str() {
                    // WARN: some weird cases where the tool predicts a non-stopped ORF:
                    // WARN: R146001_manual_scaffold_1.p1    102     2199
                    // WARN: sizes -> 144,117,332,112,120,117,138,78,66,270,246,137,79,156,87 = 2199
                    // WARN: record -> manual_scaffold_1 189532046 189543938 R146001 60 + 189532148 189543941
                    // WARN: would be out-of-bounds if we add +3 in this case
                    "+" => {
                        if orf_end + 3 > gp.end {
                            log::warn!(
                                "WARN: translationAi predicted a non-stop ORF: {orf:?} for {gp:?}"
                            );
                            orf_end = gp.end
                        } else {
                            orf_end += 3
                        }
                    }
                    "-" => {
                        //
                        if orf_start - 3 < SCALE - gp.end {
                            log::warn!(
                                "WARN: translationAi predicted a non-stop ORF: {orf:?} for {gp:?}"
                            );
                            orf_start = gp.end
                        } else {
                            orf_start -= 3;
                        }
                    }
                    _ => panic!("ERROR: unexpected strand value: {}", strand),
                }

                // INFO: retrieving the reference gene prediction record
                // INFO: since indexing groups exact similar records
                // INFO: we safely assume ref gp record could be applied to all queries
                let ref_id = format!("{}.p{}", id, orf_idx + 1);
                let ref_line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    ref_id,
                    start,
                    stop,
                    start_score,
                    stop_score,
                    strand,
                    orf_start,
                    orf_end,
                    gp.chrom
                );

                accumulator.insert(ref_line);

                if let Some(queries) = queries {
                    // INFO: queries are the orfs in the current record
                    for query in queries {
                        let query_id = format!("{}.p{}", query, orf_idx + 1);
                        let query_line = format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            query_id,
                            start,
                            stop,
                            start_score,
                            stop_score,
                            strand,
                            orf_start,
                            orf_end,
                            gp.chrom
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

/// Processes TAI (Translation AI) predictions when input sequences **were indexed**
/// during a prior extraction step. This function leverages a pre-computed binary
/// index to efficiently map numeric IDs back to their original genomic queries.
///
/// It operates in parallel, parsing the prediction file and using the numeric ID
/// from each line to retrieve the original `GenePred` record and its associated queries.
/// It then performs coordinate translation and inflates the results, generating
/// a record for each original query. These records are then collected in a
/// thread-safe accumulator and written to the output file.
///
/// # Arguments
///
/// * `predictions` - A `String` containing the full content of the TAI prediction output.
/// * `records` - A `DashMap<String, HashMap<String, GenePred>>` holding the reference
///   gene prediction records, keyed by chromosome and then by ID.
/// * `index` - A `PathBuf` representing the path to the binary index file from the extraction step.
/// * `accumulator` - A `DashSet<String>` used as a thread-safe set to collect the
///   output lines before writing to the final file.
/// * `writer` - A `BufWriter<File>` used for writing the final output to a file.
///
/// # Returns
///
/// This function does not return a value. It writes its output directly to the provided `writer`.
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
/// use std::collections::HashMap;
/// use std::io::BufWriter;
/// use std::fs::File;
/// use dashmap::DashMap;
/// use dashset::DashSet;
///
/// // Dummy types for the example
/// struct GenePred;
///
/// // Assume a function `indexed` exists with the correct signature.
/// // indexed(predictions: String, records: DashMap<String, HashMap<String, GenePred>>, ...);
///
/// let predictions = "example predictions...".to_string();
/// let records = DashMap::new();
/// let index_path = PathBuf::from("path/to/extract.index");
/// let accumulator = DashSet::new();
/// let writer = BufWriter::new(File::create("output.txt").unwrap());
///
/// // Calling the function would look like this:
/// // indexed(predictions, records, index_path, accumulator, writer);
/// ```
#[allow(unused_assignments)]
#[deprecated(note = "Functionality replaced by extract indexing")]
fn _indexed(
    predictions: String,
    records: DashMap<String, HashMap<String, GenePred>>,
    index: PathBuf,
    accumulator: DashSet<String>,
    mut writer: BufWriter<File>,
    sequences: HashMap<smol_str::SmolStr, Vec<u8>>,
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
                let mut stop = orf_parts[1].parse::<u64>().unwrap_or_else(|_| {
                    panic!("ERROR: failed to parse stop position from ORF: {}", orf);
                });
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
                let mut chr_records = records.get_mut(&chr).unwrap_or_else(|| {
                    panic!(
                        "ERROR: chromosome from {} not found in sequences -> {}!",
                        cannonical_id, chr
                    );
                });

                let gp = chr_records.get_mut(&cannonical_id).unwrap_or_else(|| {
                    panic!(
                        "ERROR: id not found in BED, this is a bug -> {}!",
                        cannonical_id
                    );
                });

                let (mut orf_start, mut orf_end) = gp.map_absolute_cds(start, stop);

                // WARN: skipping unreliable ORFs for the current alignment
                if orf_start == 0 && orf_end == 0 {
                    log::debug!(
                        "WARN: ORF start and end are zero for ID: {}.p{}, skipping!",
                        _id,
                        orf_idx
                    );
                    continue;
                }

                let mut stop_codon = "";

                // INFO: before adding up we need to check the stop codon to see if its a real one or not
                let sequence = sequences
                    .get(&cannonical_id.to_smolstr())
                    .unwrap_or_else(|| {
                        panic!("ERROR: sequence not found for {cannonical_id:?}");
                    });
                let mut orf_sequence = sequence.clone();

                let start_codon =
                    std::str::from_utf8(&sequence[start as usize..(start + 3) as usize])
                        .unwrap_or_else(|e| {
                            panic!("ERROR: failed to parse start codon -> {e} -> {sequence:?} from {cannonical_id:?} using {orf_start:?}");
                        });

                __check_start_codon(start_codon, &cannonical_id, start);

                // INFO: stop is inclusive, so we add 3 to include the stop codon
                match strand.as_str() {
                    // WARN: some weird cases where the tool predicts a non-stopped ORF:
                    // WARN: R146001_manual_scaffold_1.p1    102     2199
                    // WARN: sizes -> 144,117,332,112,120,117,138,78,66,270,246,137,79,156,87 = 2199
                    // WARN: record -> manual_scaffold_1 189532046 189543938 R146001 60 + 189532148 189543941
                    // WARN: would be out-of-bounds if we add +3 in this case
                    "+" => {
                        if orf_end + 3 > gp.end {
                            log::warn!(
                                "WARN: translationAi predicted a non-stop ORF: {orf:?}
                                with mapped coords {orf_start}:{orf_end} for {gp:?}
                                because {orf_end} + 3 > {}",
                                gp.end,
                            );
                            // orf_end = gp.end

                                stop_codon =
                                    std::str::from_utf8(&sequence[(stop - 3) as usize..(stop) as usize])
                                        .unwrap_or_else(|e| {
                                            panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {cannonical_id:?} using {orf_end:?}");
                                        });

                                log::debug!("WARN: non-stop ORF stop_codon picked -> {:?}", stop_codon);
                                orf_sequence = sequence[start as usize..stop as usize].to_vec();
                        } else {
                            stop_codon =
                                std::str::from_utf8(&sequence[(stop) as usize..(stop + 3) as usize])
                                    .unwrap_or_else(|e| {
                                        panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {cannonical_id:?} using {orf_end:?}");
                                    });

                            if !STOP_CODONS.contains(&stop_codon) {
                                // WARN: if stop_codon is not cannonical this is probably a case where the tool is wrong
                                log::debug!(
                                    "WARN: stop codon is not TAA, TAG, or TGA -> {:?} from {:?} using {:?}",
                                    stop_codon,
                                    &cannonical_id,
                                    stop
                                );

                                // INFO: taking stop as last base and going back 2 nt to capture last codon
                                stop_codon =
                                    std::str::from_utf8(&sequence[(stop - 3) as usize..(stop) as usize])
                                        .unwrap_or_else(|e| {
                                            panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {cannonical_id:?} using {orf_end:?}");
                                        });

                                dbg!("WARN: non-stop ORF stop_codon picked -> {:?}", stop_codon);
                                orf_sequence = sequence[start as usize..stop as usize].to_vec();
                            } else {
                                // INFO: stop_codon is cannonical so we can safely add 3 to the end
                                orf_end += 3;
                                orf_sequence = sequence[start as usize..(stop + 3) as usize].to_vec();
                                stop += 3;
                            }
                        }
                    }
                    "-" => {
                        // WARN: scalling down gp.end because orf_start represents scaled orf_end
                        if orf_start - 3 < SCALE - gp.end {
                            log::warn!(
                                "WARN: translationAi predicted a non-stop ORF: {orf:?}
                                with mapped coords {orf_start}:{orf_end} for {gp:?}
                                because {orf_start} - 3 < {}",
                                SCALE - gp.end
                            );
                            // orf_start = SCALE - gp.end

                                stop_codon =
                                    std::str::from_utf8(&sequence[(stop - 3) as usize..(stop) as usize])
                                        .unwrap_or_else(|e| {
                                            panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {cannonical_id:?} using {orf_end:?}");
                                        });
                                log::debug!("WARN: non-stop ORF stop_codon picked -> {:?}", stop_codon);
                                orf_sequence = sequence[start as usize..stop as usize].to_vec();
                        } else {
                            // orf_start -= 3;
                            stop_codon =
                                std::str::from_utf8(&sequence[(stop) as usize..(stop + 3) as usize])
                                    .unwrap_or_else(|e| {
                                        panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {cannonical_id:?} using {orf_end:?}");
                                    });

                            if !STOP_CODONS.contains(&stop_codon) {
                                // WARN: if stop_codon is not cannonical this is probably a case where the tool is wrong
                                log::debug!(
                                    "WARN: stop codon is not TAA, TAG, or TGA -> {:?} from {:?} using {:?}",
                                    stop_codon,
                                    &cannonical_id,
                                    stop
                                );

                                // INFO: taking stop as last base and going back 2 nt to capture last codon
                                stop_codon =
                                    std::str::from_utf8(&sequence[(stop - 3) as usize..(stop) as usize])
                                        .unwrap_or_else(|e| {
                                            panic!("ERROR: failed to parse stop codon -> {e} -> {sequence:?} from {cannonical_id:?} using {orf_end:?}");
                                        });

                                log::debug!("WARN: non-stop ORF stop_codon picked -> {:?}", stop_codon);
                                orf_sequence = sequence[start as usize..stop as usize].to_vec();
                            } else {
                                // INFO: stop_codon is cannonical so we can safely add 3 to the end
                                orf_start -= 3;
                                orf_sequence = sequence[start as usize..(stop + 3) as usize].to_vec();
                                stop += 3;
                            }
                        }
                    }
                    _ => panic!("ERROR: unexpected strand value: {}", strand),
                }

                let pep = translate(&orf_sequence);
                let inner_stops = scan_stops(orf_sequence);

                // INFO: queries are the orfs in the current record
                for query in queries {
                    let query_id = format!("{}.p{}", query, orf_idx + 1);
                    let query_line = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        gp.chrom,
                        start,
                        stop,
                        query_id,
                        strand,
                        start_score,
                        stop_score,
                        orf_start,
                        orf_end,
                        start_codon,
                        stop_codon,
                        inner_stops,
                        pep
                    );

                    accumulator.insert(query_line);
                }
            }
        });

    accumulator.into_iter().for_each(|line| {
        writeln!(writer, "{}", line).unwrap();
    });
}

/// Translates a sequence into amino acids.
///
/// # Arguments
///
/// * `sequence` - The sequence to translate.
///
/// # Returns
///
/// The translated sequence as a string.
///
/// # Panics
///
/// This function will panic if:
/// - A codon is not a valid codon.
///
/// # Example
///
/// ```rust
/// let sequence = vec![b'A', b'T', b'G', b'C'];
///
/// let aa = translate(&sequence);
/// ```
fn translate(sequence: &[u8]) -> String {
    let mut aa = String::new();

    for codon in sequence.chunks(3) {
        if codon.len() != 3 {
            break;
        }

        if codon.iter().any(|&b| !is_unambiguous_dna_base(b)) {
            aa.push('X');
            continue;
        }

        let amino_acid = translate_codon(codon);
        if amino_acid == "X" {
            panic!(
                "ERROR: codon -> {:?} is not a valid codon from sequence -> {:?}!",
                std::str::from_utf8(codon).unwrap(),
                std::str::from_utf8(sequence).unwrap()
            );
        }
        aa.push_str(amino_acid);
    }

    aa
}

/// Checks if a base is unambiguous DNA.
///
/// # Arguments
///
/// * `b` - The base to check.
///
/// # Returns
///
/// A boolean indicating whether the base is unambiguous DNA.
fn is_unambiguous_dna_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Translates a codon into an amino acid.
///
/// # Arguments
///
/// * `codon` - The codon to translate.
///
/// # Returns
///
/// The amino acid corresponding to the given codon.
///
/// # Panics
///
/// This function will panic if:
/// - The codon is not a valid codon.
///
/// # Example
///
/// ```rust
/// let codon = vec![b'A', b'T', b'G'];
///
/// let amino_acid = translate_codon(codon);
/// ```
fn translate_codon(codon: &[u8]) -> &'static str {
    for (table_codon, amino_acid) in &CODON_TABLE {
        if codon == *table_codon {
            return amino_acid;
        }
    }

    "X" // INFO: unknown codon
}

/// Scans a sequence for stop codons
///
/// This function scans a sequence for stop codons by checking if the windows
/// of size 3 contain any of the expected stop codons (TAA, TAG, or TGA).
/// The function returns the number of stop codons found.
///
/// # Arguments
///
/// * `sequence` - A `Vec<u8>` containing the sequence to scan.
///
/// # Returns
///
/// A `usize` representing the number of stop codons found.
fn scan_stops(sequence: Vec<u8>) -> usize {
    sequence
        .windows(3)
        .enumerate()
        .filter(|(idx, window)| idx % 3 == 0 && STOP_CODONS_BYTES.contains(window))
        .count()
}

/// Checks if the start codon is correct
///
/// This function checks if the start codon is correct by comparing it to the
/// expected value of ATG. If the start codon is not ATG, a warning is logged.
///
/// # Arguments
///
/// * `start_codon` - A `String` slice containing the start codon.
/// * `id` - A `String` slice containing the ID of the sequence.
/// * `orf_start` - A `u64` representing the start position of the ORF.
///
/// # Panics
///
/// This function will panic if:
/// - The start codon is not long enough.
/// - The start codon is not ATG.
fn __check_start_codon(start_codon: &str, id: &str, orf_start: u64) {
    if start_codon.len() < 3 {
        panic!(
            "ERROR: start codon is not long enough -> {start_codon:?} from {id:?} using {orf_start:?}"
        );
    }

    if start_codon != START_CODON {
        log::warn!(
            "WARN: start codon is not ATG -> {:?} from {:?} using {:?}",
            start_codon,
            &id,
            orf_start
        );
    }
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
    records: &HashMap<String, GenePred>,
    outdir: &Path,
    index: &std::collections::HashMap<u32, Vec<String>>,
) -> (PathBuf, HashMap<smol_str::SmolStr, Vec<u8>>) {
    let seqs = fasta_to_hashmap(fasta)
        .unwrap_or_else(|e| panic!("ERROR: failed to parse FASTA file -> {e}"));

    let fmt = outdir.join(format!(
        "{}.fmt.fa",
        fasta.file_stem().unwrap().to_str().unwrap()
    ));
    let mut writer = create_fasta(&fmt, "fa")
        .unwrap_or_else(|| panic!("ERROR: could not create file {:?}", fmt));

    // INFO: not using further deduplication (DashMap in HL version) -> seqs are already deduplicated by extract
    let mut accumulator = Vec::with_capacity(seqs.len());

    seqs.iter().for_each(|(header, seq)| {
        // INFO: >0 -> 0; opens index { 0 : [ R1_chr1, R2_chr1 ] }
        let u_header = header.parse::<u32>().unwrap_or_else(|e| {
            panic!("ERROR: could not parse header -> {e}");
        });

        // INFO: [ R1_chr1, R2_chr1 ] -> only use first index to build fmt_header
        let queries = index
            .get(&u_header)
            .unwrap_or_else(|| panic!("ERROR: could not find index match for id {header}"));

        let gp = records.get(&queries[0]).unwrap_or_else(|| {
            panic!("ERROR: could not find bed records for {} using index {u_header} that unroll to {queries:?}", queries[0]);
        });

        let mut start = gp.start;
        let mut end = gp.end;

        match gp.strand {
            Strand::Forward => {}
            Strand::Reverse => {
                let tmp = start;
                start = config::SCALE - end;
                end = config::SCALE - tmp;
            }
        }
        accumulator.push((
            seq,
            format!(
                ">{}:{}-{}({})({})(0, 0,)",
                gp.chrom, start, end, gp.strand, header
            ),
        ));
    });

    accumulator.into_iter().for_each(|(seq, header)| {
        writeln!(writer, "{}", &header).unwrap_or_else(|e| {
            panic!("ERROR: failed to write header -> {e}");
        });
        writer.write_all(seq).unwrap_or_else(|e| {
            panic!("ERROR: failed to write sequence -> {e}");
        });
        writeln!(writer).unwrap_or_else(|e| {
            panic!("ERROR: failed to write newline -> {e}");
        });
    });

    (fmt.clone(), seqs)
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
fn _read_tai_index(index: &PathBuf) -> HashMap<String, Vec<String>> {
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
fn _encode_for_tai(header: &String) -> Option<(u16, &str)> {
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
fn _get_header_from_values(values: &Vec<BedColumnValue>, header: &String) -> String {
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

/// Parses a FASTA file into a hash map.
///
/// # Arguments
///
/// * `f` - The path to the FASTA file.
///
/// # Returns
///
/// A `HashMap<SmolStr, Vec<u8>>` instance containing the parsed FASTA data.
///
/// # Panics
///
/// This function will panic if:
/// - The FASTA file cannot be read.
///
/// # Example
///
/// ```rust
/// use std::path::PathBuf;
/// use std::fs::File;
/// use std::io::BufWriter;
///
/// let fasta = PathBuf::from("path/to/fasta.fa");
///
/// let seqs = parse_fa(&fasta);
/// ```
pub fn fasta_to_hashmap<F: AsRef<Path>>(
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
        let header = from_utf8(&entry[..header_end])
            .unwrap()
            .trim()
            .split(' ')
            .next()
            .unwrap();
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
