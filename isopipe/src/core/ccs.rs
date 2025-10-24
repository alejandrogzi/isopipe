use crate::{
    config::*,
    consts::*,
    core::pbindex,
    executor::{job::Job, manager::ParallelExecutor},
    mv,
};

use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Read, Write},
    path::{Path, PathBuf},
};

/// Run ccs
///
/// # Arguments
/// * `step` - The pipeline step being processed.
/// * `config` - The configuration settings for the pipeline.
/// * `input_dir` - The directory containing input files.
/// * `step_output_dir` - The directory where output files will be written.
/// * `prefix` - The prefix to be used for output files.
///
/// # Returns
/// A vector of jobs to be executed.
///
/// # Example
/// ```
/// let jobs = __ccs(
///     &step,
///     &config,
///     &input_dir,
///     &step_output_dir,
///     "sample".to_string(),
/// );
/// ```
pub fn ccs(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &Path,
    executor: &mut ParallelExecutor,
) -> Vec<Job> {
    let prefix = config.get_data_prefix();

    let mut jobs = Vec::new();
    let mut require_pbi = Vec::new();

    let fields = config.get_step_custom_fields(step, vec![CHUNK, REPORT]);
    let args = config.get_step_args(
        step,
        vec![
            INPUT_DIR, PREFIX, OUTPUT_DIR, CHUNK, MEMORY, TIME, REPORT, NUM_CORES,
        ],
    );

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
        let chunk_size = fields[0]
            .parse::<usize>()
            .expect("ERROR: Failed to parse chunk size");
        let bam = entry.path();

        for chunk_idx in 0..chunk_size {
            let chunk_idx = chunk_idx + 1;

            let out_bam = step_output_dir.join(format!(
                "{}.{}.ccs.{}.bam",
                prefix,
                bam.file_stem()
                    .unwrap_or_else(|| {
                        panic!("ERROR: failed to get name from bam: {}", bam.display())
                    })
                    .to_string_lossy(),
                chunk_idx
            ));

            let chunks = format!("--chunk {}/{}", chunk_idx, chunk_size);
            let report = format!(
                "--report-file {}/{}_{}.txt",
                step_output_dir.display(),
                fields[1],
                chunk_idx
            );

            let job = Job::new()
                .task(PipelineStep::Ccs)
                .arg(bam.to_str().expect("ERROR: failed to convert path to str"))
                .arg(
                    out_bam
                        .to_str()
                        .expect("ERROR: failed to convert path to str"),
                )
                .arg(&chunks)
                .arg(&args)
                .arg(&report);

            jobs.push(job)
        }

        // WARN: need to check if bam has a .pbi file -> if not, run pbindex
        let mut pbi = bam.clone();
        pbi.set_extension("bam.pbi");
        if !pbi.exists() {
            log::warn!(
                "WARN: pbi file not found for {}, generating index...",
                bam.display()
            );

            require_pbi.push(bam.clone());
        }
    }

    if !require_pbi.is_empty() {
        pbindex::pbindex(require_pbi, config, executor, step_output_dir);
    }

    log::info!("INFO [STEP 1]: Pre-processing completed -> Running...");

    jobs
}

/// Merges ccs reports into a single file
///
/// # Arguments
/// * `ccs_output_dir` - The directory containing the ccs output files.
/// * `prefix` - The prefix to be used for the output file.
/// * `keep_temp` - Whether to keep the temporary ccs output directory.
///
/// # Example
/// ```
/// let ccs_output_dir = PathBuf::from("/path/to/ccs/output");
/// let prefix = "sample";
/// let keep_temp = true;
///
/// __join_ccs_reports(&ccs_output_dir, &prefix, &keep_temp);
/// ```
pub fn __join_ccs_reports(ccs_output_dir: &Path, prefix: &str, keep_temp: bool) {
    let reports = scan_dir(&ccs_output_dir.to_path_buf(), "txt");
    let mut accumulator: HashMap<String, u64> = HashMap::new();

    let out_reports = ccs_output_dir.join("reports");
    std::fs::create_dir_all(&out_reports).unwrap_or_else(|e| {
        panic!(
            "Failed to create reports directory in ccs output dir {} -> {e}",
            ccs_output_dir.display()
        )
    });

    for report in reports {
        let mut file = File::open(&report).expect("Failed to open report file");
        let mut contents = String::new();
        file.read_to_string(&mut contents)
            .unwrap_or_else(|e| panic!("Failed to read report file {} -> {e}", report.display()));

        for line in contents.lines() {
            if !line.is_empty() && !line.starts_with("Additional") && !line.starts_with("Exclusive")
            {
                let mut fields = line.split(':');

                let key = fields
                    .next()
                    .unwrap_or_else(|| panic!("Failed to get key from line: {line}"))
                    .trim_matches(' ');
                let value = fields
                    .next()
                    .unwrap_or_else(|| panic!("Failed to get value from line: {line}"))
                    .split("(")
                    .next()
                    .unwrap_or_else(|| {
                        panic!("Failed to get value from line while splitting: {line}")
                    })
                    .trim_matches(' ')
                    .parse::<u64>()
                    .unwrap_or_else(|_| panic!("Failed to parse value from line: {line}"));

                // INFO: if key already exists, add to its value
                if let Some(current_value) = accumulator.get_mut(key) {
                    *current_value += value;
                } else {
                    accumulator.insert(key.to_string(), value);
                }
            }
        }

        mv!(&report, &out_reports);
    }

    let mut writer = BufWriter::new(
        File::create(ccs_output_dir.join(format!("{prefix}.merged.ccs.report"))).unwrap_or_else(
            |e| {
                panic!(
                    "Failed to create merged ccs report file {} -> {e}",
                    ccs_output_dir.display()
                )
            },
        ),
    );

    __write_header(&mut writer, &mut accumulator, ccs_output_dir);

    for (key, value) in accumulator {
        writer
            .write_all(format!("{key}: {value}\n").as_bytes())
            .unwrap_or_else(|e| {
                panic!(
                    "Failed to write to merged ccs report file {} -> {e}",
                    ccs_output_dir.display()
                )
            });
    }

    if !keep_temp {
        crate::rm!(ccs_output_dir.join("reports"));
    }
}

/// Writes the header to the ccs report
///
/// # Arguments
/// * `writer` - The writer to write to.
/// * `accumulator` - The accumulator to write the header to.
/// * `ccs_output_dir` - The directory containing the ccs output files.
///
/// # Example
/// ```rust, no_run
/// let mut writer = BufWriter::new(File::create("ccs_report.txt").unwrap());
/// let mut accumulator: HashMap<String, u64> = HashMap::new();
/// let ccs_output_dir = PathBuf::from("/path/to/ccs/output");
///
/// __write_header(&mut writer, &mut accumulator, &ccs_output_dir);
/// ```
pub fn __write_header(
    writer: &mut impl Write,
    accumulator: &mut HashMap<String, u64>,
    ccs_output_dir: &Path,
) {
    let input = accumulator
        .get("ZMWs input")
        .unwrap_or_else(|| panic!("Failed to get input from accumulator"));
    let passes = accumulator
        .get("ZMWs pass filters")
        .unwrap_or_else(|| panic!("Failed to get passes from accumulator"));
    let rejects = accumulator
        .get("ZMWs fail filters")
        .unwrap_or_else(|| panic!("Failed to get rejects from accumulator"));

    writer
        .write_all(
            format!(
                "ZMWs input: {:?}\n\nZMWs pass filters: {:?} ({:.2}%)\nZMWs fail filters: {:?} ({:.2}%)\n\n",
                input, passes, passes * 100 / input, rejects, rejects * 100 / input
            )
            .as_bytes(),
        )
        .unwrap_or_else(|e| {
            panic!(
                "Failed to write to merged ccs report file {:?} -> {e}",
                ccs_output_dir
            )
        });

    accumulator.remove("ZMWs input");
    accumulator.remove("ZMWs pass filters");
    accumulator.remove("ZMWs fail filters");
}
