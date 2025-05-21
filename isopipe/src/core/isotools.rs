use config::CHUNK_SIZE;
use iso_fusion::lib_iso_fusion;
use iso_polya::{
    cli::AparentArgs,
    core::apa::{calculate_polya, create_joblist, write_bed, RAM_PER_SITE},
};
use iso_split::lib_iso_split;
use isotools::lib;
use packbed::par_reader;
use rayon::prelude::*;

use std::path::PathBuf;
use std::sync::Arc;

use crate::executor::manager::ParallelExecutor;
use crate::{config::*, consts::*, executor::job::Job};

/// Run isotools iso-fusion
///
/// # Arguments
/// * `step` - The pipeline step
/// * `config` - The configuration
/// * `input_dir` - The input directory
/// * `step_output_dir` - The output directory
///
/// # Returns
/// A vector of jobs
///
/// # Example
/// ```
/// let jobs = iso_fusion(&step, &config, &input_dir, &step_output_dir);
/// ```
pub fn iso_fusion(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
) -> Vec<Job> {
    let mut non_cannonical = true;
    let parts = config
        .get_step_args(
            step,
            vec![INPUT_DIR, PREFIX, OUTPUT_DIR, MEMORY, TIME, PREFIX],
        )
        .split(" ")
        .map(String::from)
        .filter(|s| !s.is_empty())
        .collect::<Vec<String>>();

    for category in CLUSTERING_CATEGORIES {
        if *category == "lq" {
            continue;
        }

        let mut args = Vec::new();

        let query = format!(
            "{}/{}.{}.{}",
            input_dir.display(),
            CU_ALN,
            category,
            CORR_MINIMAP_GOOD_BED
        );

        if !std::path::Path::new(&query).exists() {
            log::warn!("WARN: {} does not exist, skipping...", query);
            continue;
        }

        args.extend(vec![String::from("--query"), query]);

        let prefix = step_output_dir.join(format!("{}", category));
        args.extend(vec![String::from("--prefix"), prefix.display().to_string()]);

        std::fs::create_dir_all(&prefix).expect(&format!(
            "ERROR: Failed to create directory {}",
            prefix.display()
        ));

        args.extend(parts.clone());

        if *category == "singletons" {
            args.extend(vec![
                String::from("--suffix"),
                SINGLETONS.to_string(),
                String::from("--colorize"),
                SGN_COLOR.to_string(),
            ]);
        }

        let _ = lib_iso_fusion(Arc::new(args));
        non_cannonical = false;
    }

    log::info!("INFO [STEP 7]: Pre-processing completed -> Running...");

    if non_cannonical {
        __build_non_cannonical_fusions(input_dir, step_output_dir, parts);
    }

    let jobs = __aggregate_fusions(step_output_dir);
    return jobs;
}

/// Aggregate fusions from all categories into a single file.
///
/// # Arguments
/// * `step_output_dir` - The output directory
///
/// # Returns
/// A vector of jobs
///
/// # Example
/// ```
/// let jobs = aggregate_fusions(&step_output_dir);
/// ```
fn __aggregate_fusions(step_output_dir: &PathBuf) -> Vec<Job> {
    FUSION_TYPES
        .iter()
        .map(|ty| {
            let pattern = if *ty == "fusions" {
                format!("{0}/*/{1}.bed", step_output_dir.display(), ty)
            } else {
                format!("{0}/*/*.{1}.bed", step_output_dir.display(), ty)
            };
            let output = format!("{0}/fusions.{1}.bed", step_output_dir.display(), ty);
            Job::from(format!("cat {} > {}", pattern, output))
        })
        .collect()
}

/// Build non-cannonical fusions
///
/// # Arguments
///
/// * `input_dir` - The input directory
/// * `step_output_dir` - The output directory
/// * `parts` - The parts to be used
///
/// # Example
///
/// ```rust, no_run
/// let input_dir = PathBuf::from("/path/to/input_dir");
/// let step_output_dir = PathBuf::from("/path/to/step_output_dir");
/// let parts = vec![String::from("--part1"), String::from("--part2")];
///
/// __build_non_cannonical_fusions(&input_dir, &step_output_dir, parts);
/// ```
fn __build_non_cannonical_fusions(
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
    parts: Vec<String>,
) {
    log::warn!(
        "WARN: No cannonical jobs found to run for isotools fusion in {} -> trying to grab any .corrected.good.bed!",
        input_dir.display()
    );

    for (idx, entry) in std::fs::read_dir(input_dir)
        .expect("Failed to read assets directory")
        .flatten()
        .filter(|entry| {
            entry
                .path()
                .file_name()
                .and_then(|ext| ext.to_str())
                .map(|name| name.ends_with(CORR_MINIMAP_GOOD_BED))
                .unwrap_or(false)
        })
        .enumerate()
    {
        let mut args = Vec::new();
        let query = entry.path();

        args.extend(vec![
            String::from("--query"),
            query.to_string_lossy().to_string(),
        ]);

        let prefix = step_output_dir.join(format!("reads_{}", idx));
        args.extend(vec![String::from("--prefix"), prefix.display().to_string()]);

        std::fs::create_dir_all(&prefix).expect(&format!(
            "ERROR: Failed to create directory {}",
            prefix.display()
        ));

        args.extend(parts.clone());
        let _ = lib_iso_fusion(Arc::new(args));
    }
}

/// Run iso-polya aparent
///
/// # Arguments
///
/// * `executor` - The parallel executor
/// * `config` - The configuration
/// * `input_dir` - The input directory
/// * `step_output_dir` - The output directory
/// * `step` - The pipeline step
///
/// # Example
///
/// ```rust, no_run
/// let executor = ParallelExecutor::new();
/// let config = Config::new();
/// let input_dir = PathBuf::from("/path/to/input_dir");
/// let step_output_dir = PathBuf::from("/path/to/step_output_dir");
/// let step = PipelineStep::new();
///
/// iso_polya_aparent(&executor, &config, &input_dir, &step_output_dir, &step);
/// ```
fn __iso_polya_aparent(
    executor: &mut ParallelExecutor,
    config: &Config,
    step_output_dir: &PathBuf,
    step: &PipelineStep,
    bed: &String,
    twobit: &String,
) -> PathBuf {
    let mem = CHUNK_SIZE as f32 * RAM_PER_SITE * 1024.0;
    let package = config.get_package_from_step(step);

    let args = vec![
        String::from("--bed"),
        bed.clone(),
        String::from("--twobit"),
        twobit.clone(),
        String::from("--outdir"),
        step_output_dir.display().to_string(),
    ];

    let accumulator = calculate_polya(AparentArgs::from(args))
        .expect("ERROR: Failed to calculate polyA apparent");
    let jobs = create_joblist(&accumulator);

    log::info!("INFO [STEP 9a]: Pre-processing completed -> Running APARENT...");
    executor.__para(config, &step.to_unique_str(), &jobs, 1, mem as u32, package);

    merge_aparent(step_output_dir.clone())
}

/// Run isotools on final ORF-called data
/// [iso-classify, iso-intron, iso-polya, iso-utr]
///
/// # Arguments
///
/// * `step` - The pipeline step
/// * `config` - The configuration
/// * `input_dir` - The input directory
/// * `step_output_dir` - The output directory
/// * `executor` - The parallel executor
///
/// # Returns
///
/// A vector of jobs
///
/// # Example
///
/// ```rust
/// let jobs = polish(&step, &config, &input_dir, &step_output_dir, &executor);
/// ```
pub fn polish(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
    executor: &mut ParallelExecutor,
) -> Vec<Job> {
    __aggregate_orfs(input_dir);
    let jobs = Vec::new();

    let mut args = config
        .get_step_args(step, vec![INPUT_DIR, OUTPUT_DIR, GENOME])
        .split(" ")
        .map(String::from)
        .filter(|s| !s.is_empty())
        .collect::<Vec<String>>();

    let twobit = config.get_step_custom_fields(step, vec![GENOME])[0].clone();
    let bed = format!("{}/{}", input_dir.display(), ORF_OUTPUT);
    let aparent = __iso_polya_aparent(executor, config, step_output_dir, step, &bed, &twobit);

    log::info!("INFO [STEP 9b]: Pre-processing completed -> Polishing...");

    args.extend(vec![
        String::from("--query"),
        bed,
        String::from("--aparent"),
        aparent.display().to_string(),
        String::from("--twobit"),
        twobit,
        step_output_dir.display().to_string(),
    ]);

    // TODO: move "coding.fusions.fusions.orf.bed" from orf step to final dir
    lib(args);

    return jobs;
}

/// Aggregate ORF data [only free + fake]
///
/// # Arguments
///
/// * `input_dir` - The input directory
///
/// # Example
///
/// ```rust, no_run
/// let input_dir = PathBuf::from("/path/to/input_dir");
/// __aggregate_orfs(&input_dir);
/// ```
fn __aggregate_orfs(input_dir: &PathBuf) {
    // INFO: review file has _RVW tag!
    let cmd = format!(
        "cat {}/*fakes*bed {}/*review*bed >> {}/*free*bed && rm {}/*fakes*bed {}/*review*bed",
        input_dir.display(),
        input_dir.display(),
        input_dir.display(),
        input_dir.display(),
        input_dir.display()
    );

    let _ = std::process::Command::new("bash")
        .arg("-c")
        .arg(cmd.clone())
        .output()
        .expect(&format!("ERROR: Failed to aggregate ORF data -> {}", cmd));
}

/// Merge chunks of APARENT data
///
/// # Arguments
///
/// * `outdir` - The output directory
///
/// # Example
///
/// ```rust, no_run
/// let outdir = PathBuf::from("/path/to/output_dir");
/// merge_aparent(&outdir);
/// ```
fn merge_aparent(outdir: PathBuf) -> PathBuf {
    let assets = outdir.join(APARENT_CHUNKS);
    let mut beds = Vec::new();

    for entry in std::fs::read_dir(assets.clone())
        .expect("ERROR: Failed to read assets directory")
        .flatten()
    {
        let path = entry.path();
        if let Some(ext) = path.extension() {
            match ext.to_str() {
                Some("bed") => beds.push(path),
                _ => {}
            }
        }
    }

    let bed = par_reader(beds).expect("ERROR: Failed to merge bed files");
    let bed_dest = assets.join(APARENT_OUTPUT);
    write_bed(bed_dest.clone(), bed);

    log::info!("INFO: Merged chunks and cleaning...");
    for entry in std::fs::read_dir(assets).expect("ERROR: Failed to read assets directory") {
        if let Ok(entry) = entry {
            let path = entry.path();
            if path
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .starts_with("polya")
            {
                let _ = std::fs::remove_file(path);
            }
        }
    }

    log::info!("SUCCESS: APPARENT finished successfully!");
    return bed_dest;
}

/// Run isotools iso-split
///
/// # Arguments
/// * `step` - The pipeline step
/// * `config` - The configuration
/// * `input_dir` - The input directory
/// * `step_output_dir` - The output directory
///
/// # Returns
/// A vector of jobs
///
/// # Example
/// ```
/// let jobs = iso_fusion(&step, &config, &input_dir, &step_output_dir);
/// ```
pub fn iso_split(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
) -> PathBuf {
    log::info!(
        "INFO [ISO-SPLIT]: Preparing to split files in {} and writing to {}/chunks",
        input_dir.display(),
        step_output_dir.display()
    );

    let chunk_val = crate::numerical!(config.get_step_custom_field(step, CHUNK) => usize)
        .expect("Missing or invalid chunk size");
    let chunks = format!("--chunks {}", chunk_val);
    let outdir = step_output_dir.join(CHUNKS);

    let entries = std::fs::read_dir(input_dir)
        .expect("ERROR: Failed to read input directory")
        .flatten()
        .map(|entry| entry.path())
        .filter(|path| {
            path.file_name()
                .and_then(|name| name.to_str())
                .map(|filename| {
                    filename.ends_with(FASTQ_GZ)
                        || filename.ends_with(FQ_GZ)
                        || filename.ends_with(HQ_FASTA_GZ)
                        || filename.ends_with(HQ_FA_GZ)
                        || filename.ends_with(SGN_FASTA_GZ)
                        || filename.ends_with(SGN_FA_GZ)
                })
                .unwrap_or(false)
        })
        .collect::<Vec<_>>();

    if entries.is_empty() {
        log::error!(
            "ERROR: could not find any .fa/.fq in {}!",
            input_dir.display()
        );
        std::process::exit(1);
    } else {
        log::info!(
            "INFO [ISO-SPLIT]: Found {} files in {}",
            entries.len(),
            input_dir.display(),
        );
    }

    let (fq, fa): (Vec<_>, Vec<_>) = entries
        .into_iter()
        .partition(|path| path.ends_with(FASTQ_GZ) || path.ends_with(FQ_GZ));

    process_fa(&fa, &chunks, &outdir).expect("ERROR: something failed while processing .fa files!");
    process_fq(&fq, &chunks, &outdir).expect("ERROR: something failed while processing .fq files!");

    return outdir;
}

/// Process a collection of fasta files
///
/// # Arguments
/// * `files` - Collection of fasta files
/// * `chunks` - Number of chunks
/// * `outdir` - Output directory
///
/// # Example
///
/// ```rust, no_run
/// process_fa(&files, &chunks, &outdir);
/// ```
fn process_fa(
    files: &[PathBuf],
    chunks: &str,
    outdir: &PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
    if files.is_empty() {
        return Ok(());
    }

    for file in files {
        let suffix = file
            .file_stem()
            .expect("Could not get file stem")
            .to_str()
            .expect("ERROR: Invalid UTF-8 in file name");

        let args = vec![
            format!("--file {}", file.display()),
            chunks.to_string(),
            format!("--outdir {}", outdir.display()),
            format!("--suffix {}", suffix),
            "--threads 16".to_string(),
        ];

        let _ = lib_iso_split(args);
    }

    Ok(())
}

/// Process a collection of fastq files
///
/// # Arguments
/// * `files` - Collection of fasta files
/// * `chunks` - Number of chunks
/// * `outdir` - Output directory
///
/// # Example
///
/// ```rust, no_run
/// process_fq(&files, &chunks, &outdir);
/// ```
fn process_fq(
    files: &[PathBuf],
    chunks: &str,
    outdir: &PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
    if files.is_empty() {
        return Ok(());
    }

    log::info!(
        "INFO [ISO-SPLIT]: Splitting {} in chunks of size {}",
        files.len(),
        chunks
    );

    files.par_iter().for_each(|file| {
        let suffix = file
            .file_stem()
            .expect("ERROR: Could not get file stem")
            .to_str()
            .expect("ERROR: Invalid UTF-8 in file name");

        let args = vec![
            format!("--file {}", file.display()),
            chunks.to_string(),
            format!("--outdir {}", outdir.display()),
            format!("--suffix {}", suffix),
            "--threads 1".to_string(),
        ];

        let _ = lib_iso_split(args);
    });

    Ok(())
}
