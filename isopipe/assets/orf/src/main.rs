use clap::{self, Parser};
use log::{Level, error, info};
use simple_logger::init_with_level;

use orf::{
    blast::run_blast,
    cli::{Args, Commands},
    merge::merge,
    tai::run_tai,
};

// /beegfs/home/agi/orf/target/release/orf blast -e /beegfs/projects/hillerlab/genome/src/ORFTree/.venv/bin/orfipy --fasta /beegfs/home/agi/tmp_chunk_chr16:2.fa -a /beegfs/home/agi/tmp_chunk_chr16:2.bed --db /projects/hillerlab/genome/data/uniref/TEMP_DMND/swissprot_vertebrates.dmnd --outdir results

fn main() {
    let start = std::time::Instant::now();

    let args = Args::parse();

    init_with_level(args.level).unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();

    match args.command {
        Commands::Blast(args) => run_blast(args),
        Commands::Tai(args) => run_tai(args),
        Commands::Merge(args) => merge(args),
    }

    let elapsed = start.elapsed();
    info!("Elapsed time: {:.3?}", elapsed);
}

// fn main() {
//     use std::io::Write;

//     let dmd = std::path::PathBuf::from(
//         "/Users/alejandrogzi/Documents/projects/isopipe/isopipe/assets/orf/test/orfs.pep.dedup.diamond",
//     );
//     let index = std::path::PathBuf::from(
//         "/Users/alejandrogzi/Documents/projects/isopipe/isopipe/assets/orf/test/orfs.pep.dedup.index",
//     );

//     let index = orf::blast::read_index(&index);
//     let predictions = packbed::reader(&dmd)
//         .unwrap_or_else(|e| panic!("ERROR: failed to read blast predictions file -> {e}"));

//     let mut writer = std::io::BufWriter::new(
//         std::fs::File::create(dmd.with_extension("dmd.result")).unwrap_or_else(|e| {
//             panic!("ERROR: failed to create output file for blast results -> {e}");
//         }),
//     );

//     let accumulator = DashMap::new();

//     // INFO: filtering repeated blast hits by percent_identity -> preserving best
//     // WARN: using a transition collection to retain the best blast record based on % aligned
//     predictions
//         .par_lines()
//         .filter(|line| !line.starts_with('#'))
//         .for_each(|line| {
//             let parts: Vec<&str> = line.split('\t').collect();

//             // INFO: 16 sp|Q9QX47|SON_MOUSE 100 497 0 0 42 538 1089 1585 1.20e-163 515
//             let id = parts[0].parse::<u32>().unwrap_or_else(|_| {
//                 panic!("ERROR: failed to parse ID from line: {}", line);
//             });

//             let data = BlastRecord::from_parts(&parts);

//             // INFO: using a transition collection to retain the best blast record based on % aligned
//             accumulator
//                 .entry(id)
//                 .and_modify(|existing_data: &mut BlastRecord| {
//                     if data.blast_pid > existing_data.blast_pid {
//                         *existing_data = data.clone();
//                     }
//                 })
//                 .or_insert(data);
//         });

//     // INFO: inflate results!
//     accumulator.iter_mut().for_each(|mut record| {
//         let (id, data) = record.pair_mut();

//         // INFO: unpacking index reference -> queries
//         // INFO: for each query all blast records
//         // INFO: { index_id : [(read_id: u16, chr_bytes: [u8; chr_len], orf: u16, subseq_orf: u16, seq_len: usize)] }
//         // INFO: { 0 : [(read_id: 5903, chr_bytes: [16, 32], orf: 1, subseq_orf: 3, seq_len: 350)] }
//         let queries = index.get(id).unwrap_or_else(|| {
//             panic!("ERROR: no queries found for ID: {}", id);
//         });

//         for query in queries.into_iter() {
//             let (read_id, chr, orf, subseq_orf, seq_len, start, end) = query;

//             let query_id = format!(
//                 "R{}_chr{}.p{}@{}",
//                 read_id,
//                 from_utf8(&chr).unwrap(),
//                 orf,
//                 subseq_orf
//             );

//             // INFO: updating blast data with % algined + id
//             data.set_percent_aligned((data.blast_alignment_len as f32 / *seq_len as f32) * 100.0);
//             data.set_id(query_id);

//             writer
//                 .write_all(
//                     format!(
//                         "{}\t{}\t{:.2}\t{:e}\t{}\t{}\t{:2}\t{}\t{}\n",
//                         data.blast_id,
//                         data.blast_idx_id,
//                         data.blast_pid,
//                         data.blast_e_value,
//                         data.blast_offset,
//                         data.blast_alignment_len,
//                         data.percent_aligned,
//                         start,
//                         end
//                     )
//                     .as_bytes(),
//                 )
//                 .unwrap_or_else(|e| {
//                     panic!("ERROR: failed to write blast record to file -> {e}");
//                 });
//         }
//     });
// }
