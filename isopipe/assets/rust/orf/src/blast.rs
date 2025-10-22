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

use config::bed_to_custom_struct_collection;
use hashbrown::HashMap;
use isopipe::config::depure;
use packbed::{reader as map_to_string, record::GenePred};

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::str::from_utf8;

use crate::{blast::core::deduplicate, cli::BlastArgs, utils::*};

pub mod core;

const ORF_PEP: &str = "orfs.pep.fa";
const RESULT: &str = "result";
pub const VENV: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../py/orfipy/.venv/bin/activate"
);

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
    let dir = &args.common.outdir.join("orf");
    std::fs::create_dir_all(dir)
        .unwrap_or_else(|e| panic!("ERROR: could not create directory -> {e}!"));

    let pep = run_orfipy(&args.common.fasta, dir);
    let records = get_bed(&args.common.alignments);

    let orfs = deduplicate(
        &pep,
        true,
        args.orf_min_len,
        args.orf_min_percent,
        "M".as_bytes(),
    );

    let index = args.common.index;

    let table = get_table(orfs, records, &index, dir, args.tai);
    __run_diamond(table, &args.database, dir);

    if !args.common.keep_temp {
        isopipe::depure!(dir, "result");
    }
}

/// Gets the table of ORF records and their corresponding lines
///
/// # Arguments
///
/// * `mapper` - The mapper of ORF records and their sequences
/// * `bed` - The bed file
/// * `outdir` - The output directory
/// * `tai` - The path to the translationAi output
///
/// # Returns
///
/// A `HashMap` of `usize` and their corresponding `Vec<String>`
///
/// # Example
///
/// ```rust
/// let table = get_table(mapper, bed, outdir, tai).unwrap();
/// ```
pub fn get_table(
    mapper: HashMap<HashHead, Vec<OrfRecord>>,
    mut bed: HashMap<String, GenePred>,
    index: &PathBuf,
    outdir: &Path,
    tai: Option<PathBuf>,
) -> HashMap<usize, Vec<String>> {
    let chr = get_chr_from_path(index);
    let index = extract::read::read_index(index, &chr);

    // INFO: sequence index -> vec of constructed lines
    let mut table = HashMap::new();

    // INFO: read tai values if provided
    // WARM: tai table key fmt -> R1231_chr1:26054102-26060105
    let mut tai = if let Some(tai) = tai {
        read_tai_table(tai)
    } else {
        HashMap::new()
    };

    let filename = &outdir.join("orf");
    let mut pep = create_file(filename, "dedup.pep")
        .unwrap_or_else(|| panic!("ERROR: could not create file {:?}", filename));

    let mut unique_tai_idx = 0;

    // INFO: iterate over mapper -> write index + mapper to records and then tai if provided
    for (head, records) in mapper.iter() {
        let idx = records[0].seq_idx; // INFO: sequence index

        // INFO: table holds seq idx -> vec of ORF records
        table.entry(idx).or_insert(Vec::new());

        // WARN: record id is extract index -> 0_ORF.1 [1-10](+)
        for record in records {
            unique_tai_idx += 1;

            let mut id_parts = record.id.split("_ORF.");
            let record_extract_id = id_parts.next().unwrap_or_else(|| {
                panic!(
                    "ERROR: failed to parse cannonical ID from header: {}",
                    record.id
                );
            });

            // INFO: now id_parts holds ORF number + nested; eg. 0_.ORF.1@2 -> 1@2 -> OR1#NE2
            let tags = format!(
                "__OR{}",
                id_parts
                    .next()
                    .unwrap_or_else(|| {
                        panic!("ERROR: failed to parse tags from header: {}", record.id);
                    })
                    .replace("@", "#NE")
            );

            let gp = bed.get_mut(record_extract_id).unwrap_or_else(|| {
                panic!(
                    "ERROR: could not find gene prediction for ID: {}",
                    record_extract_id
                );
            });

            // INFO: extract index to ids -> 0 : [R123_chr12, R124_chr12]
            let record_cannonical_ids = index
                .get(&record_extract_id.parse::<u32>().unwrap_or_else(|e| {
                    panic!(
                        "ERROR: could not parse {:?} into u32 -> {e}",
                        record_extract_id
                    );
                }))
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: could not find sequence index for ID: {}",
                        record_extract_id
                    )
                });

            let (orf_start, orf_end) = gp.map_absolute_cds(record.start as u64, record.end as u64);

            for id in record_cannonical_ids {
                let record_key = format!("{}:{}-{}", id, orf_start, orf_end);

                // INFO: only important thing to inherit from orfipy is orf_type; strand is inherited from bed!
                if tai.contains_key(&record_key) {
                    // INFO: ORF predicted by translationAi -> merge with orfipy and remove from tai
                    let tai_prediction = tai.get(&record_key).unwrap_or_else(|| {
                        panic!(
                            "ERROR: could not find tai prediction for key: {}. This is a bug!",
                            record_key
                        );
                    });

                    if tai_prediction.relative_orf_start != record.start
                        || tai_prediction.relative_orf_end != record.end
                        || tai_prediction.start as u64 != orf_start
                        || tai_prediction.end as u64 != orf_end
                        || (record.start_codon != "NA"
                            && tai_prediction.start_codon != record.start_codon)
                        || (record.stop_codon != "NA"
                            && tai_prediction.stop_codon != record.stop_codon)
                    {
                        panic!(
                            "ERROR: tai prediction does not match orfipy prediction for key: {}. This is a bug -> {:?} vs {:?}!",
                            record_key, tai_prediction, record
                        );
                    }

                    let line = format!(
                        "{}\t{}\t{}\t{}{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        gp.chrom,
                        orf_start,
                        orf_end,
                        id,   // INFO: R123_chr12
                        tags, // #OR1#NE2
                        gp.strand,
                        tai_prediction.start_score,
                        tai_prediction.stop_score,
                        record.start,
                        record.end,
                        record.start_codon,
                        tai_prediction.stop_codon, // INFO: on purpose to avoid orfipy NA stop
                        tai_prediction.inner_stops,
                        record.orf_type
                    );
                    table.entry(idx).or_insert(Vec::new()).push(line);

                    // INFO: remove from tai -> tai retains tai-only
                    tai.remove(&record_key);
                } else {
                    // INFO: ORF not predicted by translationAi -> only skipped if it is not complete
                    if record.orf_type != OrfType::Complete
                        && record.orf_type != OrfType::CompleteNested
                    {
                        continue;
                    };

                    // INFO: orfipy complete or complete-nested -> follow tai fmt + orf_type + strand
                    let line = format!(
                        "{}\t{}\t{}\t{}{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        gp.chrom,
                        orf_start,
                        orf_end,
                        id,   // INFO: R123_chr12
                        tags, // #OR1#NE2
                        gp.strand,
                        0.0, // INFO: no tai start score
                        0.0, // INFO: no tai stop score
                        record.start,
                        record.end,
                        record.start_codon,
                        record.stop_codon,
                        1, // INFO: complete or complete-nested only have 1 stop codon, otherwise would not be complete
                        record.orf_type
                    );
                    table.entry(idx).or_insert(Vec::new()).push(line);
                }
            }
        }

        let handle = table.get(&idx).unwrap_or_else(|| {
            panic!("ERROR: could not find table for head -> {:?}", head);
        });

        if handle.is_empty() {
            println!(
                "WARN: empty table for head -> {:?}. Orfipy ORF not supported by tai + TP or FP",
                idx
            );
            table.remove(&idx);
        } else {
            let hh = format!(
                ">{}\n{}",
                idx,
                from_utf8(&head.seq).unwrap_or_else(|_| {
                    panic!("ERROR: failed to convert key to UTF-8 -> {:?}", head.seq);
                })
            );
            writeln!(pep, "{}", hh).unwrap_or_else(|e| {
                panic!("ERROR: failed to write record to file -> {e} -> {:?}", hh);
            });
        }
    }

    if !tai.is_empty() {
        for (_, prediction) in tai.iter() {
            unique_tai_idx += 1;

            let rc = format!(">{}\n{}", unique_tai_idx, prediction.aa);
            writeln!(pep, "{}", rc).unwrap_or_else(|e| {
                panic!("ERROR: failed to write record to file -> {e} -> {:?}", rc);
            });

            // INFO: R221_chrX.p2 -> R221_chrX__UT2
            let id = prediction.id.replace(".p", "__UT");

            let line = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                prediction.chr,
                prediction.start,
                prediction.end,
                id,
                prediction.strand,
                prediction.start_score,
                prediction.stop_score,
                prediction.relative_orf_start,
                prediction.relative_orf_end,
                prediction.start_codon,
                prediction.stop_codon,
                prediction.inner_stops,
                "UN", // INFO: unknown orf_type -> can be guessed with codons but its helpful for tai-only
            );
            table.entry(unique_tai_idx).or_insert(Vec::new()).push(line);
        }
    }

    table
}

/// Read a bed file and convert it to a hashmap of genepred records
///
/// # Arguments
///
/// * `bed` - The path to the bed file
///
/// # Returns
///
/// * `HashMap<String, GenePred>` - A hashmap of genepred records
///
/// # Example
///
/// ```rust
/// let bed = get_bed(&args.bed);
/// ```
fn get_bed(bed: &PathBuf) -> HashMap<String, GenePred> {
    let records = bed_to_custom_struct_collection::<GenePred>(
        map_to_string(bed)
            .unwrap_or_else(|e| panic!("ERROR: failed to read BED file -> {e}"))
            .into(),
        config::BedColumn::Name,
        config::BedOperation::SplitName(config::BIG_SEP, 0), // INFO: R9834_chr16__FC37#TC0#PA0#PR0#IY887 or 0
    )
    .unwrap_or_else(|e| panic!("ERROR: failed construct BED to GenePred collection -> {e}"));

    let combined: HashMap<String, GenePred> = records
        .into_iter() // consumes DashMap
        .flat_map(|(_, inner)| inner) // drop outer key, keep inner HashMap
        .collect();

    combined
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
fn run_orfipy(fasta: &Path, dir: &Path) -> PathBuf {
    let cmd = format!(
        "source {} && orfipy {} --pep {} --partial-5 --partial-3 --strand f --include-stop --min 100 --ignore-case --outdir {}",
        VENV,
        fasta.display(),
        ORF_PEP,
        &dir.display()
    );

    std::process::Command::new("bash")
        .arg("-c")
        .arg(&cmd)
        .status()
        .unwrap_or_else(|e| panic!("ERROR: failed to execute orfipy command -> {e}"));

    // INFO: checking if orfipy output is empty and make it run again!
    let outfile = dir.join(ORF_PEP);
    if !outfile.exists() || outfile.metadata().unwrap().len() == 0 {
        // INFO: forcing orfipy to run again
        std::process::Command::new("bash")
            .arg("-c")
            .arg(cmd)
            .status()
            .unwrap_or_else(|e| panic!("ERROR: failed to execute orfipy command -> {e}"));
    }

    outfile
}

/// Run diamond on the orfipy output
///
/// # Arguments
///
/// * `table` - The orfipy output
/// * `database` - The path to the translationAi fasta
/// * `outdir` - The output directory
///
/// # Returns
///
/// * `None`
///
/// # Example
///
/// ```rust
/// let table = get_table(orfs, bed, &dir, args.tai);
/// __run_diamond(table, DATABASE, &dir);
/// ```
fn __run_diamond(mut table: HashMap<usize, Vec<String>>, database: &Path, outdir: &Path) {
    let diamond = outdir.join("orf.diamond");
    let orfs = outdir.join("orf.dedup.pep");

    let cmd = format!(
        "~/Documents/binaries/diamond blastp --query {} --db {} --out {} --outfmt 6 qseqid pident qlen slen length qstart qend sstart send evalue --threads 8 --sensitive -e 1e-10",
        orfs.display(),
        database.display(),
        diamond.display()
    );

    log::info!("INFO: Executing -> {}", cmd);

    std::process::Command::new("bash")
        .arg("-c")
        .arg(cmd)
        .status()
        .unwrap_or_else(|e| panic!("ERROR: failed to execute diamond command -> {e}"));

    let mut seen = HashSet::new();
    let predictions = map_to_string(&diamond)
        .unwrap_or_else(|e| panic!("ERROR: failed to read blast predictions file -> {e}"));

    let mut writer = BufWriter::new(File::create(diamond.with_extension(RESULT)).unwrap_or_else(
        |e| {
            panic!("ERROR: failed to create output file for blast results -> {e}");
        },
    ));

    for line in predictions.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        let header = parts[0];
        let u_header = header.parse::<usize>().unwrap_or_else(|_| {
            panic!("ERROR: failed to parse header as usize -> {header} in parts: {parts:?}")
        });

        // INFO: only storing first row -> sorted by default
        if seen.contains(header) {
            continue;
        }

        seen.insert(header);

        let lines = table.get(&u_header).unwrap_or_else(|| {
            panic!("ERROR: could not find header {header} in table -> {parts:?}!")
        });

        for line in lines {
            let blast = BlastRecord::from_parts(&parts);
            let output = format!(
                "{}\t{:.2}\t{:e}\t{}\t{}\t{}\n",
                line,
                blast.blast_pid,
                blast.blast_e_value,
                blast.blast_offset,
                blast.blast_alignment_len,
                blast.percent_aligned
            );

            writer.write_all(output.as_bytes()).unwrap_or_else(|e| {
                panic!("ERROR: failed to write to file -> {e}");
            });
        }

        // INFO: table records with no blast hits remain
        table.remove(&u_header);
    }

    for (_, lines) in table.iter() {
        for line in lines {
            let output = format!("{}\t0\t1\t0\t0\t0\n", line);

            writer.write_all(output.as_bytes()).unwrap_or_else(|e| {
                panic!("ERROR: failed to write to file -> {e}");
            });
        }
    }
}
