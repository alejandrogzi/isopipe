use crate::{
    config::*,
    consts::*,
    executor::{job::Job, manager::__get_assets_dir},
};

use config::{OverlapType, Sequence, Strand, SCALE};
use iso_polya::utils::get_sequences;
use packbed::{unpack, GenePred};
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
        vec![
            INPUT_DIR,
            OUTPUT_DIR,
            MEMORY,
            TIME,
            GENOME,
            NUM_THREADS,
            CHUNK,
        ],
    );

    let twobit = PathBuf::from(config.get_step_custom_fields(step, vec![GENOME])[0].clone());
    let chunk_size = crate::numerical!(config.get_step_custom_field(step, CHUNK) => usize)
        .unwrap_or_else(|e| panic!("ERROR: could not convert chunk to num -> {e}!"));

    // INFO: looping through all fusion outputs? -> free + fakes + review [other color + tag]; fusions
    for file in FUSION_FILES {
        let bed = input_dir.join(file);

        if !bed.exists() || std::fs::metadata(&bed).unwrap().len() == 0 {
            log::warn!("WARNING: {} does not exist or its empty!", bed.display());
            continue;
        }

        let basename = file.replace(".bed", "");
        let suffix = basename
            .split(".")
            .last()
            .unwrap_or_else(|| panic!("ERROR: could not get suffix from file -> {}", file));
        let paths = extract(&bed, &twobit, step_output_dir, chunk_size, suffix);

        for (chunked_fa, chunked_bed) in paths {
            let cmd = format!(
                "source {} && {} --fasta {} --alignments {} --output_dir {} --suffix {} {}",
                ORFPY_ENV,
                executable.display(),
                chunked_fa.display(),
                chunked_bed.display(),
                step_output_dir.display(),
                chunked_bed.with_extension("orf").display(),
                args
            );

            jobs.push(Job::from(cmd));
        }
    }

    return jobs;
}

/// Extract sequences for fusion-free predicted
/// transcripts from a .2bit file into a .fa file
/// chunked
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
    chunk_size: usize,
    suffix: &str,
) -> Vec<(PathBuf, PathBuf)> {
    log::info!(
        "INFO: Extracting mapped read sequences [{}] from .2bit file...",
        reads.display()
    );

    let tmp_dir = step_output_dir.join(format!("{}_{}", SEQS, suffix));
    std::fs::create_dir_all(&tmp_dir).unwrap_or_else(|e| {
        panic!(
            "ERROR: could not creat temporary directory in {} -> {e}",
            &tmp_dir.display()
        )
    });

    let bed = unpack::<GenePred, _>(vec![reads.clone()], OverlapType::Exon, false)
        .unwrap_or_else(|e| panic!("ERROR: could not unpack reads -> {}. {e}", reads.display()));
    let (genome, _) = get_sequences(twobit.clone()).unwrap_or_else(|| {
        panic!(
            "ERROR: could not get sequences from .2bit -> {}",
            twobit.display()
        )
    });

    // INFO: define the chunk size for parallel processing
    // INFO: if chunk size > bed records -> symlink
    let paths: Vec<(PathBuf, PathBuf)> = bed
        .into_par_iter()
        .flat_map(|(chrom, records)| {
            records
                .chunks(chunk_size)
                .enumerate()
                .map(move |(chunk_id, chunk)| (format!("{}:{}", chrom, chunk_id), chunk.to_vec()))
                .collect::<Vec<_>>()
        })
        .map(|(chunk_id, transcripts)| {
            let chunk_path = tmp_dir.join(&chunk_id);
            std::fs::create_dir_all(&chunk_path).unwrap();

            let chunk_fa = chunk_path.join(format!("tmp_chunk_{}.fa", &chunk_id));
            let chunk_bed = chunk_path.join(format!("tmp_chunk_{}.bed", &chunk_id));

            let mut writer_fa = BufWriter::new(File::create(&chunk_fa).unwrap());
            let mut writer_bed = BufWriter::new(File::create(&chunk_bed).unwrap());

            let chr = chunk_id.split(':').next().unwrap_or(&chunk_id);
            for tx in transcripts {
                let seq = match tx.strand {
                    Strand::Forward => Sequence::new(
                        genome.get(chr).expect("ERROR: missing chromosome")
                            [tx.start as usize..tx.end as usize]
                            .as_ref(),
                    ),
                    Strand::Reverse => Sequence::new(
                        genome.get(chr).expect("ERROR: missing chromosome")
                            [(SCALE - tx.end) as usize..(SCALE - tx.start) as usize]
                            .as_ref(),
                    )
                    .reverse_complement(),
                };

                writeln!(writer_fa, ">{}\n{}", tx.name, seq).unwrap();
                writeln!(writer_bed, "{}", tx.line).unwrap();
            }

            (chunk_fa, chunk_bed)
        })
        .collect();

    return paths;
}
