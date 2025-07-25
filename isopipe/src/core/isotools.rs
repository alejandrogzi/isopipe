use config::CHUNK_SIZE;
use iso_polya::{
    cli::AparentArgs,
    core::apa::{calculate_polya, create_joblist, write_bed, RAM_PER_SITE},
};
use isotools::lib;
use packbed::par_reader;

use std::path::PathBuf;

use crate::{config::*, consts::*, executor::job::Job};
use crate::{executor::manager::ParallelExecutor, isotools};

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
    let mut jobs = Vec::new();

    let refs = config.get_step_custom_field(step, TOGA);
    let keep_rejected = config
        .get_step_custom_field(step, KEEP_REJECTED)
        .parse::<bool>()
        .unwrap_or(false);

    for file in std::fs::read_dir(input_dir)
        .expect("Failed to read assets directory")
        .flatten()
        .filter(|entry| {
            entry
                .path()
                .file_name()
                .and_then(|ext| ext.to_str())
                .map(|name| name.ends_with(BED))
                .unwrap_or(false)
        })
    {
        if file.path().ends_with(ALN_POLYA_REJECT) && !keep_rejected {
            continue;
        }

        // WARN: will need to change for the corrected minimap step suffix!
        let query = file.path();

        if !std::path::Path::new(&query).exists() {
            log::warn!("WARN: {} does not exist, skipping...", query.display());
            continue;
        }

        let prefix = step_output_dir.join(
            file.path()
                .file_stem()
                .unwrap_or_else(|| panic!("ERROR: could not build prefix for {:?}", file.path())),
        );

        std::fs::create_dir_all(&prefix).expect(&format!(
            "ERROR: Failed to create directory {}",
            prefix.display()
        ));

        if file.path().ends_with(ALN_POLYA_SGN) {
            let cmd = format!(
                "{} --ref {} --query {} --prefix {} --tag {} --colorize {}",
                isotools!(ISO_FUSION).display(),
                refs,
                query.display(),
                prefix.display(),
                SINGLETONS,
                SGN_COLOR
            );

            jobs.push(Job::from(cmd));
        } else {
            let cmd = format!(
                "{} --ref {} --query {} --prefix {}",
                isotools!(ISO_FUSION).display(),
                refs,
                query.display(),
                prefix.display(),
            );

            jobs.push(Job::from(cmd));
        }

        non_cannonical = false;
    }

    log::info!("INFO [STEP 7]: Pre-processing completed -> Running...");

    if non_cannonical {
        let nc_jobs = __build_non_cannonical_fusions(input_dir, step_output_dir, refs);
        return nc_jobs;
    }

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
pub fn agg_fusions(dir: &PathBuf, executor: &mut ParallelExecutor, config: &Config) {
    let jobs = FUSION_TYPES
        .iter()
        .map(|ty| {
            let (output, pattern) = if *ty == "fusions" {
                (
                    format!("{0}/{1}.bed", dir.display(), ty),
                    format!("{0}/*/{1}.bed", dir.display(), ty),
                )
            } else {
                // WARN: will aggregate free/fakes/review into fusions.free.bed!
                // INFO: fakes have :FK tag, review has :RW tag
                (
                    format!("{0}/fusions.free.bed", dir.display()),
                    format!("{0}/*/*.{1}.bed", dir.display(), ty),
                )
            };

            Job::from(format!("cat {} >> {}", pattern, output))
        })
        .collect();

    executor
        .add_jobs(jobs)
        .and_send(config, "agg-fusions", dir.clone(), 4, 8, None);
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
    refs: String,
) -> Vec<Job> {
    log::warn!(
        "WARN: No cannonical jobs found to run for isotools fusion in {} -> trying to grab any .bed!",
        input_dir.display()
    );

    let mut jobs = Vec::new();

    for (idx, entry) in std::fs::read_dir(input_dir)
        .expect("Failed to read assets directory")
        .flatten()
        .filter(|entry| {
            entry
                .path()
                .file_name()
                .and_then(|ext| ext.to_str())
                .map(|name| name.ends_with(BED))
                .unwrap_or(false)
        })
        .enumerate()
    {
        let query = entry.path();
        let prefix = step_output_dir.join(format!("reads_{}", idx));

        std::fs::create_dir_all(&prefix).expect(&format!(
            "ERROR: Failed to create directory {}",
            prefix.display()
        ));

        let cmd = format!(
            "{} --ref {} --query {} --prefix {}",
            isotools!(ISO_FUSION).display(),
            refs,
            query.display(),
            prefix.display(),
        );

        jobs.push(Job::from(cmd));
    }

    if jobs.is_empty() {
        log::error!(
            "ERROR: No .bed files found in {} to run non-cannonical fusion!",
            input_dir.display()
        );
        std::process::exit(1);
    } else {
        log::info!(
            "INFO [STEP 7]: Found {} non-cannonical fusion files in {}",
            jobs.len(),
            input_dir.display(),
        );
    }

    return jobs;
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
    executor.__para(
        config,
        &step.to_unique_str(),
        &jobs,
        1,
        mem as u32,
        Some(package),
    );

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
    executor: &mut ParallelExecutor,
) -> PathBuf {
    log::info!(
        "INFO [ISO-SPLIT]: Preparing to split files in {} and writing to {}/chunks",
        input_dir.display(),
        step_output_dir.display()
    );

    let chunks = crate::numerical!(config.get_step_custom_field(step, CHUNK) => usize)
        .expect("Missing or invalid chunk size");
    let outdir = step_output_dir.join(CHUNKS);

    // WARN: if any of the entries is FASTQ_GZ we assume the whole run is non-cannonical
    // WARN: would not make any sense combining both -> user should make two different runs and merge
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

    split_reads(entries, &chunks, &outdir, executor, config, step, input_dir);

    return outdir;
}

/// Generic processor for a collection of input
/// FASTA/FASTQ files
///
/// # Arguments
///
/// * `files` - A slice of `PathBuf`s representing the input files.
/// * `chunks` - The number of chunks to split each file into.
/// * `outdir` - The output directory where the chunks should be written.
fn split_reads(
    files: Vec<PathBuf>,
    chunks: &usize,
    outdir: &PathBuf,
    executor: &mut ParallelExecutor,
    config: &Config,
    step: &PipelineStep,
    input_dir: &PathBuf,
) {
    log::info!(
        "INFO [ISO-SPLIT]: Running iso-split on {} files!",
        files.len()
    );

    let mut jobs = Vec::new();
    let mut global_threads = 0;

    files.iter().for_each(|file| {
        let threads = if file
            .file_name()
            .and_then(|name| name.to_str())
            .map(|name| name.ends_with(FASTQ_GZ))
            .unwrap_or(false)
        {
            1
        } else {
            16
        };

        // WARN: if any of the entries is FASTQ_GZ we assume the whole run is non-cannonical
        // WARN: would not make any sense combining both -> user should make two different runs and merge
        // WARN: here enforcing the amount of threads dependent on extension!
        if global_threads < 1 {
            global_threads = threads
        }

        let bind = file.with_extension("");
        let suffix = bind
            .file_stem()
            .unwrap_or_else(|| panic!("ERROR: could not build prefix for {:?}", file));

        let cmd = format!(
            "{} {} {} {} {} {} {} {} {} {} {}",
            isotools!(ISO_SPLIT).display(),
            "--file".to_string(),
            file.display().to_string(),
            "--chunks".to_string(),
            chunks.to_string(),
            "--outdir".to_string(),
            outdir.clone().display().to_string(),
            "--suffix".to_string(),
            suffix.to_string_lossy(),
            "--threads".to_string(),
            format!("{}", threads),
        );

        jobs.push(Job::from(cmd));
    });

    executor.add_jobs(jobs).and_send(
        config,
        "iso-split",
        input_dir.clone(),
        global_threads,
        8,
        Some(config.get_package_from_step(step)),
    );
}
