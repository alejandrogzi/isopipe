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

use std::{
    fs::File,
    io::{BufWriter, Write},
};

use packbed::reader;

use crate::{cli::SambaArgs, utils::*};

const WEIGHTS: &str =
    "https://github.com/apcamargo/RNAsamba/raw/refs/heads/master/data/full_length_weights.hdf5";

pub fn run_samba(args: SambaArgs) {
    let dir = &args.outdir.join("samba");
    std::fs::create_dir_all(dir)
        .unwrap_or_else(|e| panic!("ERROR: could not create directory -> {e}!"));

    let weights = args.weights.unwrap_or_else(|| {
        let weights = dir.join("full_length_weights.hdf5");

        log::info!("INFO: downloading weights from {WEIGHTS}");
        let cmd = format!(
            "wget {WEIGHTS} -O {weights}",
            WEIGHTS = WEIGHTS,
            weights = weights.display()
        );

        std::process::Command::new("bash")
            .arg("-c")
            .arg(&cmd)
            .status()
            .unwrap_or_else(|e| panic!("ERROR: could not run command -> {cmd} -> {e}!"));

        weights
    });

    let output = dir.join(
        args.fasta
            .with_extension("tmp.samba.tsv")
            .file_name()
            .unwrap(),
    );

    let cmd = format!(
        "rnasamba classify {output} {input} {weights}",
        output = output.display(),
        input = args.fasta.display(),
        weights = weights.display()
    );

    log::info!("INFO: running command -> {cmd}");
    let status = std::process::Command::new("bash")
        .arg("-c")
        .arg(&cmd)
        .status()
        .unwrap_or_else(|e| panic!("ERROR: could not run command -> {cmd} -> {e}!"));

    if !status.success() {
        panic!("ERROR: command failed -> {cmd}");
    }

    let chr = get_chr_from_path(&args.index);
    let index = extract::read::read_index(&args.index, &chr);

    let scores =
        reader(&output).unwrap_or_else(|e| panic!("ERROR: could not read {output:?} -> {e}!"));

    let expanded_output = dir.join(
        args.fasta
            .with_extension("rnasamba.tsv")
            .file_name()
            .unwrap(),
    );
    let mut writer = BufWriter::new(
        File::open(&expanded_output)
            .unwrap_or_else(|e| panic!("ERROR: could not open {expanded_output:?} -> {e}!")),
    );

    for line in scores.lines() {
        if line.starts_with("sequence_name") {
            continue;
        }

        let parts = line.split('\t').collect::<Vec<_>>();
        let u_id = parts[0].parse::<u32>().unwrap_or_else(|e| {
            panic!(
                "ERROR: could not parse u_id from {line} -> {e}! \
                This is probably a bug in the program."
            )
        });
        let score = parts[1].parse::<f32>().unwrap_or_else(|e| {
            panic!(
                "ERROR: could not parse score from {line} -> {e}! \
                This is probably a bug in the program."
            )
        });

        let queries = index
            .get(&u_id)
            .unwrap_or_else(|| panic!("ERROR: could not find u_id {u_id} in index {index:?}!"));

        for query in queries {
            let new_line = format!("{}\t{}", query, score);
            writeln!(writer, "{}", new_line).unwrap_or_else(|e| {
                panic!(
                    "ERROR: could not write to {expanded_output:?} -> {e}! \
                    This is probably a bug in the program."
                )
            });
        }
    }

    if !args.keep_temp {
        isopipe::config::remove_any(&output);
    }
}
