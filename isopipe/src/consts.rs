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
pub const TOGA: &str = "toga";
pub const ASSEMBLY: &str = "assembly";
pub const NUM_THREADS: &str = "num-threads";
pub const NUM_CORES: &str = "num-cores";

// project-wide pub const | names
pub const ISOPIPE: &str = "isopipe";
pub const PBINDEX: &str = "pbindex";
pub const ISOTOOLS: &str = "isotools";
pub const OUTPUT: &str = "isopipe_run";
pub const ISOSEQ: &str = "isoseq";
pub const CLUSTER: &str = "cluster";
pub const SAMTOOLS: &str = "samtools";
pub const POLYA_FIRST_PASS: &str = "polya_first_pass";
pub const BEDTOOLS: &str = "bedtools";
pub const BAMTOBED: &str = "bamtobed";

// filenames
pub const NF_RUNNER: &str = "execute_joblist.nf";
pub const FOFN: &str = "all.flnc.fofn";
pub const CLUSTERED_BAM: &str = "all.clustered.bam";
pub const CLUSTERED: &str = "all.clustered";
pub const CU_ALN: &str = "all.clustered.aligned";
pub const CU_ALN_HQ_SAM: &str = "all.clustered.aligned.hq.sam";
pub const GENOME_FA: &str = "genome.fa";
pub const MERGED_BAM: &str = "merged.bam";
pub const FILTER_MINIMAP: &str = "filter_minimap_qual.perl";
pub const CORRECT_MINIMAP: &str = "correct_minimap.py";
pub const POLYA_GOOD_SAM: &str = "good.sam";
pub const CORR_MINIMAP_SAM: &str = "corrected.sam";
pub const CORR_MINIMAP_GOOD_SAM: &str = "corrected.good.sam";
pub const CORR_MINIMAP_GOOD_BED: &str = "corrected.good.bed";
pub const TRANSCRIPTS_FA: &str = "transcripts.fa";
pub const ORF: &str = "orf";
pub const ORF_EXE: &str = "orf_tree_pipe.py";

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
pub const FUSION_TYPES: &[&str] = &["free", "fusions", "review", "fakes"];
pub const FUSION_FILES: &[&str] = &[
    "fusions.free.bed",
    "fusions.fusions.bed",
    "fusions.review.bed",
    "fusions.fakes.bed",
];
