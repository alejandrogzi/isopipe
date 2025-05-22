use crate::{cat, config::*, consts::*, executor::job::Job, isotools};
use std::{fs::create_dir_all, path::PathBuf};

/// Run polya mod [3 steps]
///
/// # Arguments
///
/// * `step` - The pipeline step to run
/// * `config` - The configuration to use
/// * `input_dir` - The input directory
/// * `output_dir` - The output directory
///
/// # Returns
///
/// A vector of jobs to run
///
/// # Examples
///
/// ```rust, no_run
/// use isopipe::core::polya;
/// use isopipe::config::Config;
/// use std::path::PathBuf;
///
/// let config = Config::default();
/// let input_dir = PathBuf::from("/path/to/input");
/// let output_dir = PathBuf::from("/path/to/output");
///
/// let jobs = polya(&step, &config, &input_dir, &output_dir);
/// ```
pub fn polya(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
) -> Vec<Job> {
    let binary = isotools!(ISO_POLYA);
    let mut jobs = Vec::new();

    let parts = &step_output_dir.join(POLYA_PARTS);
    let _ = create_dir_all(parts);

    let args = config.get_step_args(
        step,
        vec![INPUT_DIR, OUTPUT_DIR, MEMORY, TIME, TOGA, ASSEMBLY],
    );

    // WARN: input_dir needs to be suffixed by /bam
    for entry in std::fs::read_dir(input_dir.join(BAM))
        .expect("Failed to read assets directory")
        .flatten()
        .filter(|entry| {
            entry
                .path()
                .file_name()
                .and_then(|ext| ext.to_str())
                .map(|name| name.ends_with(BAM))
                .unwrap_or(false)
        })
    {
        let bam = entry.path();

        // INFO: for cannonical runs this yields -> chunk*{run}.singletons
        // INFO: and outputs chunk*{run}.singletons.good.bed
        let bind = bam.with_extension("");
        let prefix = bind
            .file_stem()
            .unwrap_or_else(|| panic!("ERROR: could not build prefix for {:?}", bam));

        if std::fs::metadata(&bam).unwrap().len() == 0 {
            log::warn!("WARNING: {} its empty!", bam.display());
            continue;
        }

        let cmd = format!(
            "{} {} --bam {} {args} --prefix {} --outdir {}",
            binary.display(),
            SEGMENT,
            bam.display(),
            prefix.to_string_lossy(),
            &parts.display()
        );

        jobs.push(Job::from(cmd));
    }

    log::info!("INFO [STEP 6]: Pre-processing completed -> Running...");

    return jobs;
}

/// Build command for polya mod
///
/// # Arguments
///
/// * `filename` - The filename to use
/// * `alignment` - The alignment file to use
/// * `args` - The arguments to use
/// * `output_dir` - The output directory to use
/// * `filter` - filter_minimap_qual.py
/// * `correct` - correct_minimap.py
/// * `fields` - Custom fields from args
///
/// # Returns
///
/// A string with the command to run
///
/// # Examples
///
/// ```
/// use isopipe::core::polya::__build_cmd;
/// use std::path::PathBuf;
///
/// let filename = PathBuf::from("test.sam");
/// let alignment = PathBuf::from("alignment.sam");
/// let args = String::from("--arg1 --arg2");
/// let output_dir = PathBuf::from("/path/to/output");
/// let filter = PathBuf::from("filter.py");
/// let correct = PathBuf::from("correct.py");
/// let fields = vec![String::from("field1"), String::from("field2")];
///
/// let cmd = __build_cmd(
///  &filename,
///  &alignment,
///  &args,
///  &output_dir,
///  &filter,
///  &correct,
///  &fields,
/// );
/// ```
#[deprecated]
fn __build_cmd(
    filename: &PathBuf,
    alignment: &PathBuf,
    args: &String,
    output_dir: &PathBuf,
    filter: &PathBuf,
    correct: &PathBuf,
    fields: &Vec<String>,
    step: &PipelineStep,
    config: &Config,
) -> String {
    // INFO: will output all.clustered.aligned.{hq,lq,singletons}.{good,bad}.sam
    // INFO: script.perl {].sam --perID 96 --clip3 50 --polyAReadSuffix 30 --outdir {}/first_pass
    let first_pass = format!(
        "{} {} {} --outdir {}",
        filter.display(),
        alignment.display(),
        args,
        output_dir.join(POLYA_FIRST_PASS).display()
    );

    // INFO: script.py {toga} {}.good.sam {assembly} {].corrected.sam
    let corrected_sam = output_dir.join(filename.with_extension(CORR_MINIMAP_SAM));
    let correct_step = format!(
        "python3 {} {} {} {} {}",
        correct.display(),
        fields
            .get(0)
            .expect(&format!("ERROR: Could not find TOGA -> {:?}", fields)),
        output_dir
            .join(POLYA_FIRST_PASS)
            .join(filename.with_extension(POLYA_GOOD_SAM))
            .display(),
        fields
            .get(1)
            .expect(&format!("ERROR: Could not find assembly -> {:?}", fields)),
        corrected_sam.display()
    );

    // INFO: script.perl {}.corrected.sam --polyAReadSuffix 30 --outdir {}
    // INFO: second run will have +2 --perID -> make it compatible with 98 [quick fix]
    let per_id = config
        .get_step_custom_field(step, PER_ID)
        .parse::<u32>()
        .unwrap_or(96)
        + 2;

    let second_pass = format!(
        "{} {} -polyAReadSuffix 30 --outdir {} --perID {}",
        filter.display(),
        corrected_sam.display(),
        output_dir.display(),
        per_id
    );

    let convert = format!(
        "{} {} -i {} -bed12 > {}",
        BEDTOOLS,
        BAMTOBED,
        output_dir
            .join(filename.with_extension(CORR_MINIMAP_GOOD_SAM))
            .display(),
        output_dir
            .join(filename.with_extension(CORR_MINIMAP_GOOD_BED))
            .display()
    );

    return format!(
        "{} && {} && {} && {}",
        first_pass, correct_step, second_pass, convert
    );
}

/// Builds non-cannonical jobs
///
/// # Arguments
///
/// * `input_dir` - The input directory
/// * `args` - The arguments to use
/// * `output_dir` - The output directory
/// * `filter` - filter_minimap_qual.py
/// * `correct` - correct_minimap.py
/// * `fields` - Custom fields from args
///
/// # Returns
///
/// A vector of jobs to run
///
/// # Examples
///
/// ```rust, no_run
/// use isopipe::core::polya::__build_non_cannonical;
/// use std::path::PathBuf;
///
/// let input_dir = PathBuf::from("/path/to/input");
/// let args = String::from("--arg1 --arg2");
/// let output_dir = PathBuf::from("/path/to/output");
/// let filter = PathBuf::from("filter.py");
/// let correct = PathBuf::from("correct.py");
/// let fields = vec![String::from("field1"), String::from("field2")];
///
/// let jobs = __build_non_cannonical(
///   &input_dir,
///   &args,
///   &output_dir,
///   &filter,
///   &correct,
///   &fields,
/// );
/// ```
#[deprecated]
fn __build_non_cannonical(
    input_dir: &PathBuf,
    args: &String,
    output_dir: &PathBuf,
    filter: &PathBuf,
    correct: &PathBuf,
    fields: &Vec<String>,
    step: &PipelineStep,
    config: &Config,
) -> Vec<Job> {
    log::warn!(
        "WARN: No cannonical jobs found to run for polya in {} -> trying to grab any .sam!",
        input_dir.display()
    );

    let mut jobs = Vec::new();

    for entry in std::fs::read_dir(input_dir)
        .expect("Failed to read assets directory")
        .flatten()
        .filter(|entry| {
            entry
                .path()
                .extension()
                .and_then(|ext| ext.to_str())
                .map(|ext| ext.eq_ignore_ascii_case(SAM))
                .unwrap_or(false)
        })
    {
        let filename = PathBuf::from(entry.path().file_stem().expect(&format!(
            "ERROR: could not get stem from {}",
            entry.path().display()
        )));

        let cmd = __build_cmd(
            &filename,
            &entry.path(),
            &args,
            output_dir,
            &filter,
            &correct,
            &fields,
            step,
            config,
        );

        jobs.push(Job::from(cmd));
    }

    log::info!("INFO [STEP 6]: Pre-processing completed -> Running...");

    jobs
}

/// Merges polya .bed results into a single .bed per category
pub fn merge(input_dir: &PathBuf) {
    log::info!(
        "INFO [MERGE]: Merging polyA segmentation results in {}...",
        &input_dir.join(POLYA_PARTS).display()
    );

    let accepted = input_dir.join("all.aligned.accept.bed");
    let rejected = input_dir.join("all.aligned.reject.bed");
    let singletons = input_dir.join("all.aligned.singletons.bed");

    let files = std::fs::read_dir(input_dir.join(POLYA_PARTS))
        .expect("Failed to read assets directory")
        .flatten()
        .map(|entry| entry.path())
        .filter(|entry| {
            entry
                .file_name()
                .and_then(|ext| ext.to_str())
                .map(|name| name.ends_with(BED))
                .unwrap_or(false)
        })
        .collect::<Vec<_>>();

    if files.is_empty() {
        log::error!("ERROR: could not find any .bed in {}!", input_dir.display());
        std::process::exit(1);
    }

    // INFO: is necessary to separete singletons into a single file
    let mut sgns = Vec::new();
    let mut good = Vec::new();
    let mut bad = Vec::new();

    for path in files {
        if path.ends_with(BED_SGN_ACCEPT) {
            sgns.push(path);
        } else if path.ends_with(BED_ACCEPT) {
            good.push(path);
        } else if path.ends_with(BED_REJECT) {
            bad.push(path);
        }
    }

    cat!(&sgns, singletons);
    cat!(&good, accepted);
    cat!(&bad, rejected);
}
