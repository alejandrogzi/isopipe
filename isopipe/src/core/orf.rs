use crate::{config::*, consts::*, executor::job::Job};

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
    let args = config.get_step_custom_fields(step, vec![GENOME, ORFIPY, ORF_MIN_LEN, DATABASE]);

    let twobit = PathBuf::from(args[0].clone());
    let chunk_size = crate::numerical!(config.get_step_custom_field(step, CHUNK) => usize)
        .unwrap_or_else(|e| panic!("ERROR: could not convert chunk to num -> {e}!"));

    log::info!(
        "INFO: Merging TOGA predictions in a single file here: {}...",
        step_output_dir.display()
    );
    __merge_toga(step_output_dir, config, step);

    // INFO: looping through all fusion outputs? -> free + fakes + review [other color + tag]; fusions
    for file in FUSION_FILES {
        let bed = input_dir.join(file);

        if !bed.exists() || std::fs::metadata(&bed).unwrap().len() == 0 {
            log::warn!("WARN: {} does not exist or its empty!", bed.display());
            continue;
        }

        let basename = file.replace(".bed", "");
        let suffix = basename
            .split(".")
            .last()
            .unwrap_or_else(|| panic!("ERROR: could not get suffix from file -> {}", file));

        // INFO: inflection point -> chunking fusion files [2nd chunking step in the pipeline]
        let paths = extract(&bed, &twobit, step_output_dir, chunk_size, suffix);

        for (chunked_fa, chunked_bed) in paths {
            let chunked_dir = &chunked_bed.parent().unwrap_or_else(|| {
                panic!(
                    "ERROR: could not get parent directory for chunked_bed -> {}",
                    chunked_bed.display()
                )
            });

            let blast = format!(
                "{} blast -e {} --fasta {} --alignments {} --outdir {} --orf-min-len {} --db {}",
                ORF_RELEASE,
                args[1], // INFO: orfipy
                &chunked_fa.display(),
                &chunked_bed.display(),
                chunked_dir.display(),
                args[2], // INFO: orf_min_len,
                args[3], // INFO: database
            );

            let tai = format!(
                "{} tai --fasta {} --alignments {} --outdir {}",
                ORF_RELEASE,
                chunked_fa.display(),
                chunked_bed.display(),
                chunked_dir.display(),
            );

            jobs.push(Job::from(blast));
            jobs.push(Job::from(tai));
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

/// Merges TOGA predictions into a single file by invoking the `orf toga` command.
///
/// This function constructs and executes a shell command to run the `orf toga`
/// subcommand. It retrieves the path to the TOGA results from the provided
/// `Config` and specifies the output directory for the merged file.
///
/// # Arguments
///
/// * `step_output_dir` - A `PathBuf` representing the output directory for the current pipeline step.
/// * `config` - A reference to a `Config` struct, used to retrieve the path to TOGA results.
/// * `step` - A reference to a `PipelineStep` enum, used to identify the current step
///            and retrieve its custom fields from the `config`.
///
/// # Panics
///
/// This function will panic if:
/// - It fails to retrieve the TOGA path from the `config`.
/// - The `shell` function (which executes the command) encounters an error.
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
///
/// struct Config;
/// impl Config {
///     fn get_step_custom_field(&self, step: &PipelineStep, field_name: &str) -> PathBuf {
///         // In a real scenario, this would return a path based on config and step
///         PathBuf::from("/path/to/toga_raw_results")
///     }
/// }
///
/// enum PipelineStep {
///     MergeToga,
/// }
///
/// fn shell(cmd: String, msg: &str, tool: &str) {
///     println!("Executing: {}", cmd);
///     println!("Message: {}", msg);
///     println!("Tool: {}", tool);
/// }
///
/// const TOGA: &str = "toga_field_name"; // Placeholder for the actual field name in Config
/// const ORF_RELEASE: &str = "orf"; // Placeholder for the actual orf executable path
///
/// let output_dir = PathBuf::from("/tmp/pipeline_output/merge_toga_step");
/// let app_config = Config;
/// let current_step = PipelineStep::MergeToga;
///
/// // This would execute a command similar to:
/// // orf --path /path/to/toga_raw_results --outdir /tmp/pipeline_output/merge_toga_step
/// __merge_toga(&output_dir, &app_config, &current_step);
/// ```
fn __merge_toga(step_output_dir: &PathBuf, config: &Config, step: &PipelineStep) {
    let toga = config.get_step_custom_field(step, TOGA);
    let msg = "INFO: Merging TOGA predictions in a single file...";
    let tool = "iso-orf";

    let cmd = format!(
        "{} --path {} --outdir {}",
        ORF_RELEASE,
        toga,
        step_output_dir.display()
    );

    shell(cmd, msg, tool);
}
