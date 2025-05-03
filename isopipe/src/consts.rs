// formats
pub const BAM: &str = "bam";
pub const SAM: &str = "sam";
pub const FA: &str = "fa";
pub const TWOBIT: &str = ".2bit";
pub const FASTA_GZ: &str = "fasta.gz";

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
pub const LOG_FILE: &str = "log-file";
pub const GENOME: &str = "genome";

// project-wide pub const | names
pub const ISOPIPE: &str = "isopipe";
pub const PBINDEX: &str = "pbindex";
pub const ISOTOOLS: &str = "isotools";
pub const OUTPUT: &str = "isopipe_run";
pub const ISOSEQ: &str = "isoseq";
pub const CLUSTER: &str = "cluster";
pub const SAMTOOLS: &str = "samtools";

// filenames
pub const NF_RUNNER: &str = "execute_joblist.nf";
pub const FOFN: &str = "all.flnc.fofn";
pub const CLUSTERED_BAM: &str = "all.clustered.bam";
pub const CLUSTERED: &str = "all.clustered";
pub const CU_ALN: &str = "all.clustered.aligned";
pub const GENOME_FA: &str = "genome.fa";
pub const MERGED_BAM: &str = "merged.bam";

// manager consts
pub const ASSETS: &str = "assets";
pub const SHORT_QUEUE: &str = "short_queue";
pub const DEFAULT_MEMORY: &str = "default_memory";
pub const DEFAULT_THREADS: &str = "default_threads";

// miscellaneous constants
pub const CHARSET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
pub const RUN_ID_LEN: usize = 4;

// collections
pub const SPECIAL_PARAMETER: &[&str] = &["secondary"];
pub const CLUSTERING_CATEGORIES: &[&str] = &["hq", "lq", "singletons"];
