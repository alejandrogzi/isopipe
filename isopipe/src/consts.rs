// formats | extensions
pub const BAM: &str = "bam";
pub const SAM: &str = "sam";
pub const BED: &str = "bed";
pub const FA: &str = "fa";
pub const TWOBIT: &str = ".2bit";
pub const FASTA_GZ: &str = "fasta.gz";
pub const FASTQ_GZ: &str = "fastq.gz";
pub const FQ: &str = "fq";
pub const FQ_GZ: &str = "fq.gz";
pub const FASTA: &str = "fasta";
pub const FA_GZ: &str = "fa.gz";
pub const HQ_FA_GZ: &str = "hq.fa.gz";
pub const HQ_FASTA_GZ: &str = "hq.fasta.gz";
pub const SGN_FA_GZ: &str = "singletons.fa.gz";
pub const SGN_FASTA_GZ: &str = "singletons.fasta.gz";
pub const BED_ACCEPT: &str = "good.bed";
pub const BED_REJECT: &str = "bad.bed";
pub const BED_SGN_ACCEPT: &str = "singletons.good.bed";
pub const INDEX: &str = "index";
pub const REDUCED_BED: &str = "reduced.bed";

// config pub const keys
pub const OUTPUT_DIR: &str = "output_dir";
pub const PRIMERS: &str = "primers";
pub const INPUT_DIR: &str = "input_dir";
pub const MEMORY: &str = "memory";
pub const TIME: &str = "time";
pub const PREFIX: &str = "prefix";
pub const CHUNK: &str = "chunk";
pub const CHUNKS: &str = "chunks";
pub const REPORT: &str = "report-file";
pub const RUN_ID: &str = "run_id";
pub const LOG_FILE: &str = "log-file";
pub const GENOME: &str = "genome";
pub const TOGA: &str = "toga";
pub const ASSEMBLY: &str = "assembly";
pub const NUM_THREADS: &str = "num-threads";
pub const NUM_CORES: &str = "num-cores";
pub const PER_ID: &str = "perID";
pub const DATABASE: &str = "database";
pub const ORF_MIN_LEN: &str = "orf_min_len";
pub const ORFIPY: &str = "orfipy";
pub const KEEP_REJECTED: &str = "keep_rejected";
pub const PARALLEL_MODE: &str = "parallel_mode";
pub const KEEP_TEMP: &str = "keep_temp";
pub const CHROM_SIZES: &str = "chrom_sizes";
pub const SERVER: &str = "server";
pub const UPLOAD_PUBLIC: &str = "upload-public";
pub const USER: &str = "user";
pub const TARGET: &str = "target";
pub const WEB: &str = "web";

// project-wide pub const | names
pub const ISOPIPE: &str = "isopipe";
pub const PBINDEX: &str = "pbindex";
pub const ISOTOOLS: &str = "isotools";
pub const OUTPUT: &str = "isopipe_run";
pub const ISOSEQ: &str = "isoseq";
pub const CLUSTER: &str = "cluster";
pub const SAMTOOLS: &str = "samtools";
pub const MINIMAP2: &str = "minimap2";
pub const POLYA_FIRST_PASS: &str = "polya_first_pass";
pub const BEDTOOLS: &str = "bedtools";
pub const BAMTOBED: &str = "bamtobed";
pub const POLYA_APARENT: &str = "aparent";
pub const ISO_POLYA: &str = "iso-polya";
pub const ISO_SPLIT: &str = "iso-split";
pub const ISO_NMD: &str = "iso-nmd";
pub const ISO_FUSION: &str = "iso-fusion";
pub const SEGMENT: &str = "segment";
pub const POLYA_PARTS: &str = "parts";
pub const SEQS: &str = "seqs";
pub const FREE: &str = "free";
pub const FUSIONS: &str = "fusions";
pub const EXTRACT: &str = "extract";
pub const TAI: &str = "tai";

// filenames | suffixes
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
pub const ORF_OUTPUT: &str = "coding.fusions.free.orf.bed";
pub const APARENT_CHUNKS: &str = "TMP_CHUNKS";
pub const APARENT_OUTPUT: &str = "iso_polya_aparent.bed";
pub const SINGLETONS: &str = "SG"; // singleton tag
pub const ALN_POLYA_ACCEPT: &str = "all.aligned.accept.bed";
pub const ALN_POLYA_REJECT: &str = "all.aligned.reject.bed";
pub const ALN_POLYA_SGN: &str = "all.aligned.singletons.bed";

// manager consts
pub const ASSETS: &str = "assets";
pub const SHORT_QUEUE: &str = "short_queue";
pub const DEFAULT_MEMORY: &str = "default_memory";
pub const DEFAULT_THREADS: &str = "default_threads";

// miscellaneous constants
pub const CHARSET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
pub const RUN_ID_LEN: usize = 4;
pub const SGN_COLOR: &str = "51,153,255";
pub const HG_LOAD_BED: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/assets/c/hgLoadBed");
pub const BED_TO_BIG_BED: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/assets/c/bedToBigBed");
pub const SCHEMA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/assets/as/schema.as");
pub const ORF_CHUNKS: usize = 50000; // 50k transcripts

// tool binaries or scripts
pub const PREDICT_PY: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/assets/py/predict/predict.py");
pub const ISOTOOLS_RELEASE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../isotools/isotools/target/release"
);
pub const ORF_RELEASE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/assets/rust/orf/target/release/orf"
);
pub const PRETTY_PY: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/assets/py/pretty/pretty_descriptor.py"
);
pub const EXTRACT_RELEASE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/assets/rust/extract/target/release/extract"
);
pub const TAI_VENV: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/assets/rust/orf/tai/.venv/bin/activate"
);
pub const APARENT_PY: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../isotools/isotools/iso-polya/assets/run_aparent.py"
);

// collections
pub const SPECIAL_PARAMETER: &[&str] = &["secondary"];
pub const CLUSTERING_CATEGORIES: &[&str] = &["hq", "lq", "singletons"];
pub const FUSION_TYPES: &[&str] = &["fusions", "free", "review", "fakes"];
pub const ALN_POLYA_FILES: &[&str] = &[ALN_POLYA_SGN, ALN_POLYA_ACCEPT, ALN_POLYA_REJECT];
pub const FUSION_FILES: &[&str] = &[
    "fusions.free.bed",
    "fusions.fusions.bed",
    "fusions.review.bed",
    "fusions.fakes.bed",
];
pub const GZ_EXTENSIONS: &[&str] = &[
    ".fastq.gz",
    ".fq.gz",
    ".hq.fasta.gz",
    ".hq.fa.gz",
    ".singletons.fasta.gz",
    ".singletons.fa.gz",
];
pub const MUST_FILL: &[&str] = &[
    OUTPUT_DIR,
    INPUT_DIR,
    GENOME,
    PRIMERS,
    MEMORY,
    TIME,
    NUM_THREADS,
    DATABASE,
    ORF_MIN_LEN,
    TOGA,
    "junc-bed",
    "bigwig",
    "user",
    "server",
    "target",
    "web",
    "chrom_sizes",
];
pub const COLOR_SCHEMA: &[(&str, &str)] = &[
    ("trash.bed", "255,41,0"),           // INFO: bright-red
    ("rt_reads.bed", "233,21,237"),      // INFO: pink
    ("retention.bed", "207,125,43"),     // INFO: orange
    ("intraprimming.bed", "173,158,43"), // INFO: dark-yellow
    ("truncation.bed", "102,86,61"),     // INFO: brown
    ("fusions.bed", "128,60,171"),       // INFO: purple
];
