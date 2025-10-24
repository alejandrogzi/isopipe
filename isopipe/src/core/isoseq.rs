use crate::{config::*, consts::*, executor::job::Job};
use std::path::{Path, PathBuf};

/// Run isoseq3 refine
///
/// # Arguments
/// * `step` - The pipeline step being processed.
/// * `config` - The configuration object.
/// * `input_dir` - The directory containing the input BAM files.
/// * `step_output_dir` - The directory where the output files will be written.
/// * `prefix` - The prefix to be used for the output files.
///
/// # Returns
/// A vector of jobs to be executed.
///
/// # Example
/// ```rust, no_run
/// let jobs = __refine(
///     &step,
///     &config,
///     &input_dir,
///     &step_output_dir,
/// );
/// ```
pub fn refine(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &Path,
) -> Vec<Job> {
    let mut jobs = Vec::new();

    let args = config.get_step_args(step, vec![INPUT_DIR, OUTPUT_DIR, MEMORY, TIME, PRIMERS]);
    let fields = config.get_step_custom_fields(step, vec![PRIMERS]);

    // INFO: gzipping .report and .clips
    let report = scan_dir(input_dir, "report");
    let clips = scan_dir(input_dir, "clips");

    crate::gz!(&report);
    crate::gz!(&clips);

    // INFO: format of files: {prefix}.{name}.ccs.merged.fl.{primers}.bam
    for entry in std::fs::read_dir(input_dir)
        .expect("Failed to read assets directory")
        .flatten()
        .filter(|entry| {
            entry
                .path()
                .extension()
                .and_then(|ext| ext.to_str())
                .map(|ext| ext.eq_ignore_ascii_case(BAM))
                .unwrap_or(false)
        })
    {
        let bam = entry.path();
        let basename = bam
            .file_stem()
            .expect("ERROR: failed to get file stem")
            .to_string_lossy();

        let out_bam = step_output_dir.join(format!("{}.flnc.bam", basename));

        let job = Job::new()
            .task(*step)
            .arg(bam.to_str().expect("ERROR: failed to convert path to str"))
            .arg(&fields[0])
            .arg(
                out_bam
                    .to_str()
                    .expect("ERROR: failed to convert path to str"),
            )
            .arg(&args);

        jobs.push(job)
    }

    log::info!("INFO [STEP 3]: Pre-processing completed -> Running...");

    jobs
}

/// Run isoseq3 cluster
///
/// # Arguments
/// * `step` - The pipeline step being processed.
/// * `config` - The configuration for the pipeline.
/// * `input_dir` - The directory containing the input files.
/// * `step_output_dir` - The directory where the output files will be written.
/// * `prefix` - The prefix to use for the output files.
///
/// # Returns
/// A vector of jobs to be executed.
///
/// # Examples
/// ```rust, no_run
/// let jobs = cluster(
///     &step,
///     &config,
///     &input_dir,
///     &step_output_dir,
/// );
/// ```
pub fn cluster(
    step: &PipelineStep,
    config: &Config,
    input_dir: &Path,
    step_output_dir: &Path,
) -> Vec<Job> {
    let refine_fofn = format!("{}/*flnc.{}", input_dir.display(), BAM);
    let all_fofn = format!("{}/{}", step_output_dir.display(), FOFN);

    shell(
        format!("ls {} > {}", refine_fofn, all_fofn),
        "INFO: Grouping flnc reads...",
        CLUSTER,
        false,
    );

    // INFO: gzipping .report.csv
    let report = scan_dir(&input_dir.to_path_buf(), "csv");
    crate::gz!(&report);

    let args = config.get_step_args(step, vec![INPUT_DIR, OUTPUT_DIR, MEMORY, TIME, LOG_FILE]);
    let out_bam = format!("{}/{}", step_output_dir.display(), CLUSTERED_BAM);
    let fields = config.get_step_custom_fields(step, vec![LOG_FILE]);

    let tmp_hq_bam = step_output_dir.join("hq.bam").to_string_lossy().to_string();
    let tmp_singleton_bam = step_output_dir
        .join("singletons.bam")
        .to_string_lossy()
        .to_string();
    let hq_fasta = step_output_dir
        .join("all.clustered.hq.fasta.gz")
        .to_string_lossy()
        .to_string();
    let singleton_fasta = step_output_dir
        .join("all.clustered.singletons.fasta.gz")
        .to_string_lossy()
        .to_string();

    let merging = format!(
        "&& samtools view -h -@ 4 {out_bam} \
        | tee >(awk '{{if($1 ~ /^@/){{print; next}} for(i=12;i<=NF;i++){{if($i ~ /^is:i:1$/){{print; break}}}}}}' \
            | samtools view -@ 4 -bo {tmp_singleton_bam} -) \
        | awk '{{if($1 ~ /^@/){{print; next}} for(i=12;i<=NF;i++){{if($i ~ /^is:i:/){{split($i,a,\":\"); if(a[3]!=1){{print; break}}}}}}}}' \
        | samtools view -@ 4 -bo {tmp_hq_bam} - \
        && samtools fasta -@ 8 -c 9 {tmp_hq_bam} | gzip > {hq_fasta} \
        && samtools fasta -@ 8 -c 9 {tmp_singleton_bam} | gzip > {singleton_fasta} \
        && rm {tmp_hq_bam} {tmp_singleton_bam}"
    );

    let jobs = vec![Job::new()
        .task(*step)
        .arg(&all_fofn)
        .arg(&out_bam)
        .arg(&args)
        .arg(format!("--log-file {}/{}", &step_output_dir.display(), fields[0]).as_str())
        .arg(merging)];

    log::info!("INFO [STEP 4]: Pre-processing completed -> Running...");

    jobs
}
