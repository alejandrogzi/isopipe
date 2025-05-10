use crate::{
    config::*,
    consts::*,
    executor::{job::Job, manager::__get_assets_dir},
};

use config::{write_objs, OverlapType, Sequence, Strand, SCALE};
use dashmap::DashSet;
use iso_polya::utils::get_sequences;
use packbed::{record::Bed6, unpack};
use rayon::prelude::*;
use std::path::PathBuf;

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
            "{} --fasta {} --alignments {} --output_dir {} {}",
            executable.display(),
            fasta.display(),
            bed.display(),
            step_output_dir.display(),
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

    let fasta = step_output_dir.join(format!("{}.{}", filename, TRANSCRIPTS_FA));
    let accumulator: DashSet<String> = DashSet::new();

    let bed = unpack::<Bed6, _>(vec![reads.clone()], OverlapType::Exon, false).expect(&format!(
        "ERROR: could not unpack reads -> {}",
        reads.display(),
    ));
    let (genome, _) = get_sequences(twobit.clone()).expect(&format!(
        "ERROR: could not get sequences from .2bit -> {}",
        twobit.display(),
    ));

    bed.par_iter().for_each(|(chr, transcripts)| {
        for tx in transcripts {
            let seq = match tx.strand {
                Strand::Forward => Sequence::new(
                    genome
                        .get(chr)
                        .expect("ERROR: Could not chromosome from genome!")
                        [tx.coord.0 as usize..tx.coord.1 as usize]
                        .as_ref(),
                ),
                Strand::Reverse => Sequence::new(
                    genome
                        .get(chr)
                        .expect("ERROR: Could not read donor context!")
                        [(SCALE - tx.coord.1) as usize..(SCALE - tx.coord.0) as usize]
                        .as_ref(),
                )
                .reverse_complement(),
            };

            accumulator.insert(format!(">{}\n{}", tx.id, seq.to_string()));
        }
    });

    write_objs(
        &accumulator,
        fasta
            .to_str()
            .expect("ERROR: could not convert path to str!"),
    );

    return fasta;
}
