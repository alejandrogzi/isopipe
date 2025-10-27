use crate::{config::*, consts::*, executor::job::Job, isotools};
use std::path::{Path, PathBuf};

/// Load and concatenate decision files, convert to bigBed format, and optionally upload
///
/// # Arguments
/// * `step` - The pipeline step to run
/// * `config` - The configuration for the pipeline
/// * `input_dir` - The directory containing the input files (from polish step)
/// * `step_output_dir` - The directory to write the output files to
/// * `executor` - Parallel executor for running concatenation jobs
///
/// # Note
///
/// 'input_dir' points to the polish step output, containing chromosome-specific
/// decision directories with categorized BED files. The function concatenates
/// files across all chromosomes for each category, then converts them to bigBed
/// format for visualization.
///
/// Categories processed: pass.bed, trash.bed, rt_reads.bed, retention.bed,
/// intraprimming.bed, truncation.bed, nmd.bed, fusions.bed
///
/// # Returns
/// A vector of jobs to run for bigBed conversion and optional upload
///
/// # Example
/// ```rust, ignore
/// let jobs = load(
///     &step,
///     &config,
///     &input_dir,
///     &step_output_dir,
///     &mut executor,
/// );
/// ```
#[allow(unused_assignments)]
pub fn load(
    step: &PipelineStep,
    config: &Config,
    input_dir: &Path,
    step_output_dir: &PathBuf,
    executor: &mut crate::core::ParallelExecutor,
) -> Vec<Job> {
    let mut jobs = Vec::new();

    let chrom_sizes = config.get_step_custom_field(step, CHROM_SIZES);
    let server = config.get_step_custom_field(step, SERVER);
    let user = config.get_step_custom_field(step, USER);
    let target = config
        .get_step_custom_field(step, TARGET)
        .parse::<PathBuf>()
        .unwrap_or_else(|_| {
            panic!(
                "ERROR: could not parse target directory -> {}",
                config.get_step_custom_field(step, TARGET)
            )
        });
    let web = config
        .get_step_custom_field(step, WEB)
        .parse::<PathBuf>()
        .unwrap_or_else(|_| {
            panic!(
                "ERROR: could not parse web directory -> {}",
                config.get_step_custom_field(step, WEB)
            )
        });
    let upload = config
        .get_step_custom_field(step, UPLOAD_PUBLIC)
        .parse::<bool>()
        .unwrap_or(false);
    let toga = config.get_step_custom_field(step, TOGA);

    let package = config.get_package_from_step(step);

    // INFO: create bed and bb dirs
    for subdir in ["bed", "bb"] {
        let path = step_output_dir.join(subdir);
        std::fs::create_dir_all(&path).unwrap_or_else(|e| {
            panic!(
                "ERROR: could not create {} directory -> {:?}. {e}",
                subdir, path
            )
        });
    }

    let categories = [
        "pass.bed",
        "trash.bed",
        "rt_reads.bed",
        "retention.bed",
        "intraprimming.bed",
        "truncation.bed",
        "nmd.bed",
        "fusions.bed",
    ];

    // INFO: path would look like: {step_polish}/{chr}/decision/ -> looping for each chr
    let inner_jobs = categories
        .iter()
        .map(|target| {
            let pattern = input_dir.join("*").join("decision").join(target);
            let out = target; // INFO: output named same as input basename

            let mut cmd = format!(
                "cat {} >> {}",
                pattern.display(),
                step_output_dir.join("bed").join(out).display()
            );

            // INFO: ignore nmd and pass bc those colors are ok
            if *target != "nmd.bed" && *target != "pass.bed" {
                // INFO: gawk -i inplace -F'\t' 'BEGIN{OFS="\t"} {$9="255,0,0"; print}' rt_reads.bed
                cmd = format!(
                    "{cmd} && gawk -i inplace -F'\\t' 'BEGIN{{OFS=\"\\t\"}} {{$9=\"{}\"; print}}' {}",
                    get_color(target),
                    step_output_dir.join("bed").join(out).display()
                );

            }

            Job::from(cmd)
        })
        .collect();

    executor.add_jobs(inner_jobs).and_send(
        config,
        &step.to_unique_str(),
        step_output_dir.clone(),
        1,
        4,
        Some(package),
        Some("cat"),
    );

    // INFO: each .bed in step11_load should be converted into bigBed
    for bed in std::fs::read_dir(step_output_dir.join("bed"))
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", step_output_dir))
        .flatten()
        .filter(|e| e.path().is_file())
    {
        let input = bed.path();

        // INFO: pass.bed
        let basename = input
            .file_name()
            .unwrap_or_else(|| panic!("ERROR: could not get file name from {bed:?}"))
            .to_string_lossy()
            .to_string();

        let bind = input.with_extension("bb");
        let bb = bind
            .file_name()
            .unwrap_or_else(|| panic!("ERROR: could not get file name from {bed:?}"));

        // INFO: modifying bed file in-place with collapsed -> takes new schema
        // INFO: if bed is passing category, drop orphans first!
        let mut cmd = String::new();

        if basename == "pass.bed" {
            cmd = format!(
                "{} --bed {} --toga {} --all --outdir {} --name pass && {COLLAPSE} run --bed {} --extend --outdir {} --name {} && {BED_TO_BIG_BED} -tab -sort -extraIndex=name -as={SCHEMA} -type={BB_TYPE} {} {} {}",
                // INFO: orphan
                isotools!(ISO_ORPHAN).display(),
                input.display(),
                toga,
                step_output_dir.display(), // INFO: creates orphans/ by default and orphans/pass.orphan_free.bed + orphans.bed
                // INFO: collapse
                step_output_dir.join("orphans").join("pass.orphan_free.bed").display(), // INFO: orphan free output -> pass_orphan_free.bed
                step_output_dir.display(), // INFO: creates collapse/ by default
                basename.clone(), // INFO: e.g pass.bed -> will force output to be collapsed/pass.bed
                // INFO: bigbed
                step_output_dir.join("collapsed").join(&basename).display(), // INFO: collapsed/pass.bed
                chrom_sizes,
                step_output_dir.join("bb").join(bb).display() // INFO: bb/pass.bb
            );

            // TODO: find a way to assert that pass.orphans.bed exists
            let orphans_cmd = format!(
                " && {COLLAPSE} run --bed {} --extend --outdir {} --name {} && {BED_TO_BIG_BED} -tab -sort -extraIndex=name -as={SCHEMA} -type={BB_TYPE} {} {} {}",
                // INFO: collapse
                step_output_dir.join("orphans").join("pass.orphans.bed").display(),
                step_output_dir.display(), // INFO: creates collapse/ by default
                "orphans.bed", // INFO: will force output to be collapsed/orphans.bed
                // INFO: bigbed
                step_output_dir.join("collapsed").join("orphans.bed").display(), // INFO: collapsed/orphans.bed
                chrom_sizes,
                step_output_dir.join("bb").join("orphans.bb").display() // INFO: bb/orphans.bb
            );

            cmd = format!("{cmd}{orphans_cmd}");
        } else {
            cmd = format!(
                "{COLLAPSE} run --bed {} --extend --outdir {} --name {} && {BED_TO_BIG_BED} -tab -sort -extraIndex=name -as={SCHEMA} -type={BB_TYPE} {} {} {}",
                input.display(),
                step_output_dir.display(), // INFO: creates collapse/ by default
                basename.clone(), // INFO: e.g pass.bed
                step_output_dir.join("collapsed").join(&basename).display(), // INFO: collapsed/pass.bed
                chrom_sizes,
                step_output_dir.join("bb").join(bb).display() // INFO: bb/pass.bb
            );
        }

        if upload {
            let bb_name = format!("HLIsoClassAnnot.{}", bb.to_string_lossy()); // INFO: HLIsoClassAnnot.pass.bb

            cmd = format!(
                "{cmd} && ssh {user}@{server} mkdir -p {} && rsync -av {} {user}@{server}:{} && ssh {user}@{server} ln -sf {} {}",
                target.join(ISOPIPE).display(),
                step_output_dir.join("bb").join(bb).display(),
                target.join(ISOPIPE).join(&bb_name).display(),
                target.join(ISOPIPE).join(&bb_name).display(),
                web.join(bb_name).display(),
            );

            if basename == "pass.bed" {
                let orphans_bb_name = "HLIsoClassAnnot.orphans.bb";

                cmd = format!("{cmd} && ssh {user}@{server} mkdir -p {} && rsync -av {} {user}@{server}:{} && ssh {user}@{server} ln -sf {} {}",
                    target.join(ISOPIPE).display(),
                    step_output_dir.join("bb").join("orphans.bb").display(),
                    target.join(ISOPIPE).join(orphans_bb_name).display(),
                    target.join(ISOPIPE).join(orphans_bb_name).display(),
                    web.join(orphans_bb_name).display(),
                );
            }
        }

        log::debug!("DEBUG: executing cmd: {cmd:?}");
        jobs.push(Job::from(cmd));
    }

    jobs
}
