use crate::{config::*, consts::*, executor::job::Job};
use std::path::PathBuf;

pub fn load(
    step: &PipelineStep,
    config: &Config,
    input_dir: &PathBuf,
    step_output_dir: &PathBuf,
) -> Vec<Job> {
    let mut jobs = Vec::new();
    let chrom_sizes = config.get_step_custom_field(step, CHROM_SIZES);

    // INFO: path would look like: {step_polish}/{chr} -> looping for each chr
    for entry in std::fs::read_dir(input_dir)
        .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", input_dir))
        .flatten()
        .filter(|e| e.path().is_dir())
    {
        let subdir = entry.path(); // INFO: {chr}
        let chr = subdir
            .file_name()
            .unwrap_or_else(|| panic!("ERROR: could not get file naem from {subdir:?}"));
        let outdir = step_output_dir.join(chr);

        // INFO: expects decision/ holding output of pretty_descriptor
        for seqs in std::fs::read_dir(&subdir)
            .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", subdir))
            .flatten()
            .filter(|e| e.path().is_dir())
        {
            let decisions = seqs.path(); // INFO: {chr}/decision

            // INFO: each .bed in decision/ should be converted into bigBed
            for bed in std::fs::read_dir(&decisions)
                .unwrap_or_else(|e| panic!("ERROR: could read directory -> {:?}. {e}", decisions))
                .flatten()
                .filter(|e| e.path().is_file())
            {
                let input = bed.path();

                let bind = input.with_extension("bb");
                let bb = bind
                    .file_name()
                    .unwrap_or_else(|| panic!("ERROR: could not get file name from {bed:?}"));

                let cmd = format!(
                    "{BED_TO_BIG_BED} -as={SCHEMA} -type=12+25 {} {} {}",
                    input.display(),
                    chrom_sizes,
                    outdir.join(bb).display()
                );

                jobs.push(Job::from(cmd));
            }
        }
    }

    jobs
}
