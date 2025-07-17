use config::{BedColumn, BedColumnValue, bed_to_struct_collection};
use dashmap::{DashMap, DashSet};
use hashbrown::HashMap;
use log::{error, warn};
use packbed::{
    reader as bed_reader,
    record::{Bed6, GenePred},
};
use rayon::prelude::*;
use smol_str::ToSmolStr;

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::PathBuf;

use crate::{cli::TaiArgs, utils::*};

const PREDICTIONS: &str = "predictions";
const RESULT: &str = "tai.result";

const COLUMNS: &[BedColumn] = &[
    BedColumn::Chrom,
    BedColumn::Start,
    BedColumn::End,
    BedColumn::Name,
    BedColumn::Strand,
];

pub fn run_tai(args: TaiArgs) {
    let dir = args.common.outdir.join("tai");
    std::fs::create_dir_all(&dir)
        .unwrap_or_else(|e| panic!("ERROR: could not create directory -> {e}!"));

    let (fasta, index) = refmt(&args.common.fasta, &args.common.alignments, &dir);

    let cmd = format!(
        "translationai -I {} -t {},{} -O {}",
        fasta.display(),
        args.threshold,
        args.threshold,
        fasta.with_extension(PREDICTIONS).display()
    );

    std::process::Command::new("bash")
        .arg("-c")
        .arg(cmd)
        .status()
        .unwrap_or_else(|e| panic!("ERROR: failed to execute orfipy command -> {e}"));

    let idx = read_tai_index(&index);
    unroll_tai(idx, fasta, &args.common.alignments);
}

fn unroll_tai(index: HashMap<String, Vec<String>>, fasta: PathBuf, alignments: &PathBuf) {
    // INFO: predictions in fmt -> id st,sp,st_score,sp_score ...
    let predictions = bed_reader(fasta.with_extension(PREDICTIONS))
        .unwrap_or_else(|e| panic!("ERROR: failed to read predictions file -> {e}"));

    let output = fasta.with_extension(RESULT);
    let mut writer = BufWriter::new(
        File::create(&output)
            .unwrap_or_else(|e| panic!("ERROR: cannot create index from sequences -> {e}")),
    );

    let accumulator = DashSet::new();

    // INFO: inflate results!
    predictions
        .par_lines()
        .filter(|line| !line.starts_with('#'))
        .for_each(|line| {
            let parts: Vec<&str> = line.split('\t').collect();

            // INFO: >chr16:91343975-91360783 +) R9834_chr16__FC37#TC0#PA0#PR0#IY887) 0, 0,)
            let name = parts[0].split("(").collect::<Vec<&str>>();

            let coords = name[0]
                .strip_prefix('>')
                .and_then(|s| s.split(':').nth(1))
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to parse transcript start from header: {}",
                        parts[0]
                    );
                })
                .split('-')
                .collect::<Vec<&str>>()
                .iter()
                .map(|s| s.parse::<u64>().unwrap())
                .collect::<Vec<u64>>();

            let strand = name[1].trim_end_matches(')').to_string();
            let cannonical_id = name[2].trim_end_matches(')').to_string(); // INFO: R9834_chr16__FC37#TC0#PA0#PR0#IY887
            let id = cannonical_id.split("__").collect::<Vec<&str>>()[0]; // INFO: R9834_chr16
            // let chr = id.split('_').collect::<Vec<&str>>()[1]; // INFO: chr16

            // INFO: unpacking index reference -> queries
            // INFO: for each query all orfs in the current record!
            // WARN: leaving option because unique ids are not mapped to an index!
            let queries = index.get(id);

            for (orf_idx, orf) in parts.iter().skip(1).enumerate() {
                // INFO: 13190,13667,0.6413461491465569,0.921319767832756
                let orf_parts: Vec<&str> = orf.split(',').collect();
                if orf_parts.len() < 2 {
                    panic!("ERROR: ORF does not have enough parts to parse: {}", orf);
                }

                let start = orf_parts[0].parse::<u64>().unwrap_or_else(|_| {
                    panic!("ERROR: failed to parse start position from ORF: {}", orf);
                });
                let stop = orf_parts[1].parse::<u64>().unwrap_or_else(|_| {
                    panic!("ERROR: failed to parse stop position from ORF: {}", orf);
                });
                let start_score = orf_parts[2].parse::<f32>().unwrap_or_else(|_| {
                    panic!("ERROR: failed to parse start position from ORF: {}", orf);
                });
                let stop_score = orf_parts[3].parse::<f32>().unwrap_or_else(|_| {
                    panic!("ERROR: failed to parse stop position from ORF: {}", orf);
                });

                if start > stop {
                    panic!(
                        "ERROR: start position is greater than stop position in ORF: {}",
                        orf
                    );
                }

                // INFO: retrieving the reference gene prediction record
                // INFO: since indexing groups exact similar records
                // INFO: we safely assume ref gp record could be applied to all queries
                let ref_id = format!("{}.p{}", id, orf_idx + 1);
                let ref_line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    ref_id,
                    start,
                    stop,
                    start_score,
                    stop_score,
                    strand,
                    coords[0] + start,
                    coords[1] + stop
                );

                accumulator.insert(ref_line);

                if let Some(queries) = queries {
                    // INFO: queries are the orfs in the current record
                    // process_queries(
                    //     queries,
                    //     orf_idx,
                    //     start,
                    //     stop,
                    //     start_score,
                    //     stop_score,
                    //     strand,
                    //     &mut accumulator,
                    // );
                    for query in queries {
                        let query_id = format!("{}.p{}", query, orf_idx + 1);
                        let query_line = format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            query_id,
                            start,
                            stop,
                            start_score,
                            stop_score,
                            strand,
                            coords[0] + start,
                            coords[1] + stop
                        );

                        accumulator.insert(query_line);
                    }
                } else {
                    // here we should do: if log level is debug -> warn
                    warn!("WARN: no queries found for ID: {}", id);
                    continue;
                }
            }
        });

    accumulator.into_iter().for_each(|line| {
        writeln!(writer, "{}", line).unwrap();
    });

    // cleaning the rest of the files
}

fn refmt(fasta: &PathBuf, bed: &PathBuf, outdir: &PathBuf) -> (PathBuf, PathBuf) {
    let records = parse_bed::<Bed6>(bed, COLUMNS.to_vec());
    let seqs =
        parse_fa(fasta).unwrap_or_else(|e| panic!("ERROR: failed to parse FASTA file -> {e}"));

    let fmt = outdir.join(format!(
        "{}.fmt.fa",
        fasta.file_stem().unwrap().to_str().unwrap()
    ));
    let mut writer = create_fasta(&fmt, "fa")
        .unwrap_or_else(|| panic!("ERROR: could not create file {:?}", fmt));
    let mut index = create_index(&fmt);

    let accumulator = DashMap::new();

    records.into_par_iter().for_each(|(chr, rows)| {
        for (header, values) in rows.iter() {
            let fmt_header = get_header_from_values(values, header);
            let seq = seqs.get(&header.clone().to_smolstr()).unwrap_or_else(|| {
                panic!("ERROR: chromosome from {header} not found in sequences -> {chr}!");
            });

            accumulator
                .entry(seq)
                .or_insert_with(Vec::new)
                .push(fmt_header)
        }
    });

    // INFO: will write the first header on collection and will point other headers to it
    accumulator.into_iter().for_each(|(seq, headers)| {
        writeln!(writer, "{}", &headers[0]).unwrap();
        writer.write_all(seq).unwrap();
        writeln!(writer).unwrap();

        if headers.len() > 1 {
            let mut count = 0;
            for header in &headers {
                let (id, chr) = if let Some((read, chr)) = encode_for_tai(header) {
                    (read, chr)
                } else {
                    panic!("ERROR: failed to encode header: {}", header);
                };

                // INFO: first record -> chr byte len + chr bytes
                // INFO: chr-len chr n_ids id[ref] id1 id2
                if count < 1 {
                    index.write_all(&[chr.len() as u8]).unwrap();
                    index.write_all(chr.as_bytes()).unwrap();

                    let n_ids = headers.len() as u16;
                    index.write_all(&n_ids.to_be_bytes()).unwrap();
                }

                index.write_all(&id.to_be_bytes()).unwrap();

                count += 1;
            }
        }
    });

    (fmt.clone(), fmt.with_extension("dedup.index"))
}

// [chr_len: u8]
// [chr_bytes: [u8; chr_len]]
// [n_ids: u16]
// [id_reference: u16]
// [id_1: u16]
// [id_2: u16]
// ...
// [id_n: u16]
fn read_tai_index(index: &PathBuf) -> HashMap<String, Vec<String>> {
    let mut mapper = HashMap::new();

    let mut reader = BufReader::new(
        File::open(index).unwrap_or_else(|e| panic!("ERROR: failed to open index -> {e}")),
    );

    loop {
        let mut chr_len_buf = [0u8; 1];
        if reader.read_exact(&mut chr_len_buf).is_err() {
            break; // EOF reached cleanly
        }
        let chr_len = chr_len_buf[0] as usize;

        let mut chr_buf = vec![0u8; chr_len];
        reader.read_exact(&mut chr_buf).unwrap();
        let chr = String::from_utf8(chr_buf)
            .unwrap_or_else(|e| panic!("ERROR: failed to convert chr bytes to string -> {e}"));

        let mut id_count_buf = [0u8; 2];
        reader.read_exact(&mut id_count_buf).unwrap();
        let n_ids = u16::from_be_bytes(id_count_buf);

        let mut id_buf = [0u8; 2];
        let mut group = Vec::with_capacity(n_ids as usize);
        for _ in 0..n_ids {
            reader.read_exact(&mut id_buf).unwrap();
            let id = u16::from_be_bytes(id_buf);

            // WARN: id fmt -> R{int}_chr{chr}, skipping tags!
            let name = format!("R{}_chr{}", id, chr);
            group.push(name);
        }

        let reference = group[0].clone();
        let queries = group[1..].to_vec();

        mapper
            .entry(reference)
            .or_insert_with(Vec::new)
            .extend(queries);
    }

    mapper
}

fn encode_for_tai(header: &String) -> Option<(u16, &str)> {
    // INFO: ">chr16:45612921-45619040(+)(R6713_chr16__FC48#TC40#PA0#PR0#IY876)(0, 0,)",
    // INFO: 6713 16
    let parts: Vec<&str> = header.split('(').collect();
    if parts.len() < 2 {
        panic!(
            "ERROR: header does not have enough parts to parse: {}",
            header
        );
    }

    let id = parts[2].trim_end_matches(')');
    let id_parts: Vec<&str> = id.split('_').collect();

    let read = id_parts[0]
        .strip_prefix('R')?
        .parse::<u16>()
        .ok()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse read ID from header: {}", header);
        });
    let chr = id_parts[1].strip_prefix("chr").unwrap_or_else(|| {
        panic!("ERROR: failed to parse chromosome from header: {}", header);
    });

    Some((read, chr))
}

/// Reformats a .fa header with specific requirements
fn get_header_from_values(values: &Vec<BedColumnValue>, header: &String) -> String {
    let chr = values[0].as_str().unwrap_or_else(|| {
        panic!(
            "ERROR: start position not found for header: {} - {:?}",
            header, values
        )
    });
    let start = values[1].as_number().unwrap_or_else(|| {
        panic!(
            "ERROR: start position not found for header: {} - {:?}",
            header, values
        )
    });
    let end = values[2].as_number().unwrap_or_else(|| {
        panic!(
            "ERROR: end position not found for header: {} - {:?}",
            header, values
        )
    });
    let name = values[3].as_str().unwrap_or_else(|| {
        panic!(
            "ERROR: name not found for header: {} - {:?}",
            header, values
        )
    });
    let strand = values[4].as_str().unwrap_or_else(|| {
        panic!(
            "ERROR: strand not found for header: {} - {:?}",
            header, values
        )
    });

    let hdr = format!(">{}:{}-{}({})({})(0, 0,)", chr, start, end, strand, name);
    hdr
}
