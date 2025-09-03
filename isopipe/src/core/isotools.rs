use config::{CHUNK_SIZE, FUSION_FAKES, FUSION_FREE, FUSION_REVIEW};
use iso_polya::{
    cli::AparentArgs,
    core::apa::{calculate_polya, write_bed, RAM_PER_SITE},
};
use packbed::par_reader;

use std::{fs::remove_dir, path::PathBuf};

use crate::{config::*, consts::*, executor::job::Job};
use crate::{executor::manager::ParallelExecutor, isotools};

/// This function processes input directories to identify specific prediction files
/// (`tmp_predictions.bed`) and generates a series of commands (jobs) for the
/// `isotools iso-nmd` utility. It iterates through a hierarchical directory
/// structure, extracts chromosome and chunk information from directory names,
/// and constructs a command for each relevant file.
///
/// # Arguments
/// * `step` - A reference to a `PipelineStep` object, representing the current step in the pipeline.
///   It's used to retrieve custom configuration fields, such as `PARALLEL_MODE`.
/// * `config` - A reference to a `Config` object, which holds the overall configuration for the
///   pipeline. This is used to access step-specific custom fields.
/// * `input_dir` - A reference to a `PathBuf` indicating the root directory where the function
///   should search for input data. This directory is expected to contain subdirectories structured
///   as `{step_orf}/seqs_{suffix}` or `toga`.
/// * `step_output_dir` - A reference to a `PathBuf` specifying the directory where the output of
///   the `iso-nmd` commands will be stored.
///
/// # Returns
/// A `Vec<Job>`: A vector containing `Job` objects, where each `Job` represents a command to be
/// executed by the pipeline. These commands are generated based on the `tmp_predictions.bed`
/// files found in the input directory structure.
pub fn iso_nmd(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
) -> Vec<Job> {
    let _mode = ParallelMode::from_str(&config.get_step_custom_field(step, PARALLEL_MODE));
    let mut jobs = Vec::new();

    let keep_temp = config
        .get_step_custom_field(step, KEEP_TEMP)
        .parse::<bool>()
        .unwrap_or(false);

    if keep_temp {
        log::info!("INFO [STEP 9]: Keeping temporary files -> keep_temp set to true!");
    }

    // INFO: path would look like: {step_orf}/seqs_{suffix} or toga
    for entry in std::fs::read_dir(input_dir)
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", input_dir))
        .flatten()
        .filter(|e| e.path().is_dir())
    {
        let subdir = entry.path(); // INFO: seqs_{suffix}
        let suffix = subdir
            .file_name()
            .unwrap_or_else(|| panic!("ERROR: could not get file name from {:?}", subdir))
            .to_str()
            .unwrap_or_else(|| panic!("ERROR: could not convert file name to str"));

        // INFO: structure {step_orf}/seqs_{suffix}/{chr}:{chunk}
        for chunk in std::fs::read_dir(&subdir)
            .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", subdir))
            .flatten()
            .filter(|e| e.path().is_dir())
        {
            let chunked_dir = chunk.path(); // INFO: {chr}:{chunk}
            let chr = chunked_dir
                .file_name()
                .unwrap_or_else(|| panic!("ERROR: could not get file name from {:?}", chunked_dir))
                .to_str()
                .unwrap_or_else(|| panic!("ERROR: could not convert file name to str"))
                .split(':')
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: could not get chromosome from {:?}", chunked_dir)
                });

            for file in std::fs::read_dir(&chunked_dir)
                .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", chunked_dir))
                .flatten()
            {
                let file = file.path();

                if file.is_file() {
                    // INFO: matching specific name
                    if file.ends_with("tmp_predictions.bed") {
                        log::debug!("DEBUG: found predictions.bed for {chr:?} -> {file:?}");

                        let chr_outdir = step_output_dir.join(chr).join(suffix);
                        let mut cmd = format!(
                            "{} --ref {} --outdir {} --prefix {} && cat {} >> {}",
                            isotools!(ISO_NMD).display(),
                            file.display(),
                            chr_outdir.display(),
                            format!("tmp_{}", chunked_dir.file_name().unwrap().to_str().unwrap()),
                            file.with_extension("tsv").display(),
                            chr_outdir
                                .join(format!("{}.predictions.{}.tsv", chr, suffix))
                                .display(),
                        );

                        // INFO: takes care of step8 tmp files
                        if !keep_temp {
                            let rest = format!(
                                " && rm {}",
                                file.parent().unwrap().join("tmp*").display(),
                            );

                            cmd += &rest;
                        }

                        log::debug!("DEBUG: executing cmd: {cmd}");
                        jobs.push(Job::from(cmd));
                    }
                }
            }
        }
    }

    log::info!("INFO [STEP 9]: Pre-processing completed -> Running...");
    jobs
}

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
/// ```rust, ignore
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

    // INFO: input_dir would have per-chr .bed files with polya suffixes
    for file in std::fs::read_dir(input_dir)
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", input_dir))
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
        if file
            .path()
            .file_name()
            .and_then(|ext| ext.to_str())
            .unwrap()
            .ends_with(ALN_POLYA_REJECT)
            && !keep_rejected
        {
            continue;
        }

        // WARN: will need to change for the corrected minimap step suffix!
        let query = file.path();

        if !std::path::Path::new(&query).exists() {
            log::warn!("WARN: {} does not exist, skipping...", query.display());
            continue;
        }

        // INFO: file stem should be 'chrX_all.aligned.accept' -> dir name
        let prefix = step_output_dir.join(
            file.path()
                .file_stem()
                .unwrap_or_else(|| panic!("ERROR: could not build prefix for {:?}", file.path())),
        );

        std::fs::create_dir_all(&prefix).expect(&format!(
            "ERROR: Failed to create directory {}",
            prefix.display()
        ));

        if file
            .path()
            .file_name()
            .and_then(|ext| ext.to_str())
            .unwrap()
            .ends_with(ALN_POLYA_SGN)
        {
            let cmd = format!(
                "{} --ref {} --query {} --prefix {} --tag {} --colorize {} --recover",
                isotools!(ISO_FUSION).display(),
                refs,
                query.display(),
                prefix.display(),
                SINGLETONS,
                SGN_COLOR
            );

            log::debug!("debug: processing {file:?} with cmd: {cmd:?}...");
            jobs.push(Job::from(cmd));
        } else {
            let cmd = format!(
                "{} --ref {} --query {} --prefix {} --recover",
                isotools!(ISO_FUSION).display(),
                refs,
                query.display(),
                prefix.display(),
            );

            log::debug!("debug: processing {file:?} with cmd: {cmd:?}...");
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
pub fn agg_fusions(
    dir: &PathBuf,
    executor: &mut ParallelExecutor,
    config: &Config,
    step: &PipelineStep,
) {
    let mode = ParallelMode::from_str(&config.get_step_custom_field(step, PARALLEL_MODE));

    let jobs: Vec<Job> = match mode {
        ParallelMode::Chromosome => {
            log::info!(
                "INFO [STEP 8]: Aggregating fusions per-chr -> {}...",
                dir.display()
            );

            let mut jobs = Vec::new();

            // INFO: loop through all subdirectories in {dir}
            for entry in std::fs::read_dir(dir)
                .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", dir))
                .flatten()
                .filter(|e| e.path().is_dir())
            {
                let subdir = entry.path();

                // INFO: structure should be {step_fusion}/chr{chr}_all.aligned.accept/fusions*
                // INFO: expected: fusions.bed / fusions.free.bed; optional: fusions.review + fusions.fakes
                for file in std::fs::read_dir(&subdir)
                    .unwrap()
                    .flatten()
                    .filter(|e| e.path().is_file())
                {
                    let file = file.path();

                    if file.ends_with(FUSION_FAKES) || file.ends_with(FUSION_REVIEW) {
                        let target = file
                            .parent()
                            .unwrap_or_else(|| {
                                panic!("ERROR: could not get parent from {:?}", file)
                            })
                            .join(FUSION_FREE);

                        let cmd = format!(
                            "cat {} >> {} && rm {}",
                            file.display(),
                            target.display(),
                            file.display()
                        );

                        jobs.push(Job::from(cmd))
                    }
                }
            }

            jobs
        }
        ParallelMode::Genome => {
            let jobs = FUSION_TYPES
                .iter()
                .map(|ty| {
                    log::info!(
                        "INFO [STEP 8]: Aggregating {} in single files -> {}...",
                        ty,
                        dir.display()
                    );

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

            jobs
        }
    };

    executor
        .add_jobs(jobs)
        .and_send(config, "agg-fusions", dir.clone(), 1, 8, None, None);
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
fn iso_polya_aparent(
    step_output_dir: &PathBuf,
    bed: &String,
    twobit: &String,
    prefix: &str,
) -> Vec<Job> {
    log::info!("INFO [STEP 10a]: Pre-processing for APARENT...");

    let mut jobs = Vec::new();

    let outdir = step_output_dir.join(prefix);

    std::fs::create_dir_all(&outdir).expect(&format!(
        "ERROR: Failed to create directory {}",
        outdir.display()
    ));

    let args = vec![
        String::from("aparent"),
        String::from("--bed"),
        bed.clone(),
        String::from("--twobit"),
        twobit.clone(),
        String::from("--outdir"),
        outdir.display().to_string(),
    ];

    let accumulator = calculate_polya(AparentArgs::from(args))
        .expect("ERROR: Failed to calculate polyA apparent");

    for path in accumulator.paths.iter() {
        let cmd = format!(
            "source {} && python3 {} -p {}",
            TAI_VENV,
            APARENT_PY,
            path.as_str()
        );

        jobs.push(Job::from(cmd));
    }

    jobs
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
    let mut jobs = Vec::new();
    let mut inner_jobs = Vec::new();

    let args = config.get_step_args(
        step,
        vec![INPUT_DIR, OUTPUT_DIR, MEMORY, TIME, GENOME, NUM_THREADS],
    );

    let twobit = config.get_step_custom_fields(step, vec![GENOME])[0].clone();
    let mem = CHUNK_SIZE as f32 * RAM_PER_SITE; // INFO: will be converted to MB by executor

    // INFO: path would look like: {step_nmd}/{chr} -> looping for each chr
    for entry in std::fs::read_dir(input_dir)
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", input_dir))
        .flatten()
        .filter(|e| e.path().is_dir())
    {
        let subdir = entry.path(); // INFO: {chr}

        let chr =
            subdir.file_name().unwrap().to_str().unwrap_or_else(|| {
                panic!("ERROR: could not convert chr name to str -> {subdir:?}")
            });

        // INFO: {chr}/seqs_{suffix}
        for seqs in std::fs::read_dir(&subdir)
            .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", subdir))
            .flatten()
            .filter(|e| e.path().is_dir())
        {
            let seqs_dir = seqs.path();
            let suffix = seqs_dir
                .file_name()
                .unwrap_or_else(|| panic!("ERROR: could not get file name from {:?}", seqs_dir))
                .to_str()
                .unwrap_or_else(|| panic!("ERROR: could not convert file name to str"));

            let chr_beds = seqs_dir.join("*reads.bed");
            let chr_nmds = seqs_dir.join("*nmd.bed");

            let cat = format!(
                "cat {} > {chr}.reads.bed && cat {} > {chr}.nmd.bed",
                chr_beds.display(),
                chr_nmds.display()
            );

            std::process::Command::new("bash")
                .arg("-c")
                .arg(cat)
                .current_dir(&seqs_dir)
                .output()
                .expect("ERROR: Failed to concatenate bed files");

            if suffix == "seqs_free" {
                let bed = seqs_dir.join(format!("{chr}.reads.bed"));

                let apa_jobs = iso_polya_aparent(
                    step_output_dir,
                    &bed.to_str().unwrap().to_string(),
                    &twobit,
                    chr,
                );

                inner_jobs.extend(apa_jobs);
            }
        }
    }

    log::info!("INFO [STEP 10a]: Pre-processing completed -> Running APARENT...");

    executor.add_jobs(inner_jobs).and_send(
        config,
        "aparent",
        step_output_dir.clone(),
        1,
        mem as u32,
        None,
        None,
    );

    log::info!("INFO [STEP 10b]: Pre-processing completed -> Polishing...");

    // INFO: looping again to run polish with merged aparent -> per-chr
    for entry in std::fs::read_dir(input_dir)
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", input_dir))
        .flatten()
        .filter(|e| e.path().is_dir())
    {
        let subdir = entry.path(); // INFO: {chr}
        let chr =
            subdir.file_name().unwrap().to_str().unwrap_or_else(|| {
                panic!("ERROR: could not convert chr name to str -> {subdir:?}")
            });

        // INFO: aparent is only ran per-seqs_free-{chr} -> not running isotools on seqs_fusions!
        let apa = merge_aparent(step_output_dir.join(chr), "tmp");

        // INFO: again, we only want to run isotools on {chr}/seqs_free
        let bed = subdir.join("seqs_free").join(format!("{}.reads.bed", chr));
        let outdir = step_output_dir.join(chr);

        // INFO: need to specify cleanup to include it in the cmd!
        let tmp = subdir.join("*/tmp*"); // INFO: remove anything tmp
        let nmd = subdir.join("*/*nmd.bed"); // INFO: move nmds from free/fusions
        let fsn = subdir.join("*/seqs_fusions/*reads.bed"); // INFO: fusions reads
        let preds = subdir.join("*/*predictions*tsv"); // INFO: move nmds from free/fusions
        let reads = subdir.join("*/seqs_free/*reads.bed"); // INFO: mv reads as raw reads

        let cleaning = format!(
            "rm {} && cat {} > {} && mv {} {} && mv {} {} && mv {} {}",
            tmp.display(),
            nmd.display(),
            step_output_dir.join(chr).join("nmd.bed").display(),
            fsn.display(),
            step_output_dir.join(chr).join("fusions.bed").display(),
            preds.display(),
            step_output_dir.join(chr).display(),
            reads.display(),
            step_output_dir.join(chr).join("raw_reads.bed").display(),
        );

        let cmd = format!(
            "{} run --query {} --aparent {} --twobit {} {args} --outdir {} && {}",
            isotools!(ISOTOOLS).display(),
            bed.display(),
            apa.display(),
            twobit,
            outdir.display(),
            cleaning
        );

        log::debug!("DEBUG: executing cmd: {cmd}");
        jobs.push(Job::from(cmd));
    }

    return jobs;
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
fn merge_aparent(outdir: PathBuf, _prefix: &str) -> PathBuf {
    let assets = outdir.join(APARENT_CHUNKS);

    let mut beds = Vec::new();

    for entry in std::fs::read_dir(assets.clone())
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", assets))
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
    let bed_dest = outdir.join(APARENT_OUTPUT);

    log::debug!("DEBUG: will write all chunked .bed to -> {bed_dest:?}");

    write_bed(bed_dest.clone(), bed);

    log::info!("INFO: Merged chunks and cleaning...");
    for entry in std::fs::read_dir(&assets).expect("ERROR: Failed to read assets directory") {
        if let Ok(entry) = entry {
            let path = entry.path();
            if path
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .starts_with("tmp")
            {
                let _ = std::fs::remove_file(path);
            }
        }
    }

    remove_dir(&assets).unwrap_or_else(|e| panic!("ERROR: {e} -> could not remove {assets:?}"));

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
        None,
    );
}
