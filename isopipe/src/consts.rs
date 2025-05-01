// formats
pub const BAM: &str = "bam";
pub const SAM: &str = "sam";

// config pub const keys
pub const OUTPUT_DIR: &str = "output_dir";
pub const PRIMERS: &str = "primers";
pub const INPUT_DIR: &str = "input_dir";
pub const MEMORY: &str = "memory";
pub const TIME: &str = "time";
pub const PREFIX: &str = "prefix";
pub const CHUNK: &str = "chunk";
pub const REPORT: &str = "report-file";
pub const RUN_ID: &str = "run_id";

// project-wide pub const
pub const ISOPIPE: &str = "isopipe";
pub const PBINDEX: &str = "pbindex";
pub const ISOTOOLS: &str = "isotools";
pub const OUTPUT: &str = "isopipe_run";
pub const ISOSEQ3: &str = "isoseq3";

// filenames
pub const NF_RUNNER: &str = "execute_joblist.nf";

// managers consts
pub const ASSETS: &str = "assets";
pub const SHORT_QUEUE: &str = "short_queue";
pub const DEFAULT_MEMORY: &str = "default_memory";
pub const DEFAULT_THREADS: &str = "default_threads";

// miscellaneous constants
pub const CHARSET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
pub const RUN_ID_LEN: usize = 4;
