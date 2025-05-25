use crate::{
    config::*,
    consts::*,
    executor::{job::Job, manager::__get_assets_dir},
};

use config::{OverlapType, Sequence, Strand, SCALE};
use iso_polya::utils::get_sequences;
use packbed::{record::Bed6, unpack};
use rayon::prelude::*;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

/// Run ORF prediction on fusion files
///
/// # Arguments
///
/// * `step` - Pipeline step
/// * `config` - Pipeline configuration
/// * `input_dir` - Path to the input directory
/// * `step_output_dir` - Path to the output directory
///
/// # Returns
///
/// * Vec<Job> - List of jobs to be executed
///
/// # Example
///
/// ```rust, no_run
/// let jobs = orf(
///     &step,
///     &config,
///     &input_dir,
///     &step_output_dir,
/// );
/// ```
pub fn orf(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
) -> Vec<Job> {
    let mut jobs = Vec::new();
    let executable = __get_assets_dir().join(ORF).join(ORF_EXE);

    let args = config.get_step_args(
        step,
        vec![INPUT_DIR, OUTPUT_DIR, MEMORY, TIME, GENOME, NUM_THREADS],
    );

    let twobit = PathBuf::from(config.get_step_custom_fields(step, vec![GENOME])[0].clone());

    // INFO: looping through all fusion outputs
    for file in FUSION_FILES {
        let bed = input_dir.join(file);

        if !bed.exists() || std::fs::metadata(&bed).unwrap().len() == 0 {
            log::warn!("WARNING: {} does not exist or its empty!", bed.display());
            continue;
        }

        let filename = file.replace(".bed", "");
        let fasta = extract(&bed, &twobit, step_output_dir, filename);

        let cmd = format!(
            "source {} && {} --fasta {} --alignments {} --output_dir {} --suffix {} {}",
            ORFPY_ENV,
            executable.display(),
            fasta.display(),
            bed.display(),
            step_output_dir.display(),
            file.replace(".bed", ".orf"),
            args
        );

        jobs.push(Job::from(cmd));
    }

    return jobs;
}

/// Extract sequences for fusion-free predicted
/// transcripts from a .2bit file
///
/// # Arguments
///
/// * `reads` - Path to the reads file
/// * `twobit` - Path to the .2bit file
/// * `step_output_dir` - Path to the output directory
///
/// # Returns
///
/// * None
///
/// # Example
///
/// ```rust
/// /// let jobs = orf(
///     &step,
///     &config,
///     &input_dir,
///     &step_output_dir,
/// );
/// ```
pub fn extract(
    reads: &PathBuf,
    twobit: &PathBuf,
    step_output_dir: &PathBuf,
    filename: String,
) -> PathBuf {
    log::info!(
        "INFO: Extracting mapped read sequences [{}] from .2bit file...",
        reads.display()
    );

    let tmp_dir = step_output_dir.join(SEQS);
    std::fs::create_dir_all(&tmp_dir).unwrap_or_else(|e| {
        panic!(
            "ERROR: could not creat temporary directory in {} -> {e}",
            &tmp_dir.display()
        )
    });

    let fasta = step_output_dir.join(format!("{}.{}", filename, TRANSCRIPTS_FA));

    let bed = unpack::<Bed6, _>(vec![reads.clone()], OverlapType::Exon, false)
        .unwrap_or_else(|e| panic!("ERROR: could not unpack reads -> {}. {e}", reads.display()));
    let (genome, _) = get_sequences(twobit.clone()).unwrap_or_else(|| {
        panic!(
            "ERROR: could not get sequences from .2bit -> {}",
            twobit.display()
        )
    });

    bed.par_iter().for_each(|(chr, transcripts)| {
        let chr_fa = tmp_dir.join(format!("tmp_chunk_{}.fa", chr));
        let mut writer =
            BufWriter::new(File::create(&chr_fa).expect("Could not create temp FASTA file"));

        for tx in transcripts {
            let seq = match tx.strand {
                Strand::Forward => Sequence::new(
                    genome
                        .get(chr)
                        .unwrap_or_else(|| panic!("ERROR: Could not chromosome from genome!"))
                        [tx.coord.0 as usize..tx.coord.1 as usize]
                        .as_ref(),
                ),
                Strand::Reverse => Sequence::new(
                    genome
                        .get(chr)
                        .unwrap_or_else(|| panic!("ERROR: Could not chromosome from genome!"))
                        [(SCALE - tx.coord.1) as usize..(SCALE - tx.coord.0) as usize]
                        .as_ref(),
                )
                .reverse_complement(),
            };

            writeln!(writer, ">{}\n{}", tx.id, seq)
                .unwrap_or_else(|e| panic!("ERROR: could not write sequence {} -> {e}", seq));
        }
    });

    let seqs = std::fs::read_dir(&tmp_dir)
        .unwrap_or_else(|e| panic!("ERROR: failed to read {:?} directory -> {e}", &tmp_dir))
        .flatten()
        .map(|entry| entry.path())
        .collect::<Vec<_>>();

    crate::cat!(&seqs, &fasta);
    crate::rm!(tmp_dir);

    return fasta;
}
