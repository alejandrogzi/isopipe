use config::bed_to_custom_struct_collection;
use dashmap::DashMap;
use hashbrown::HashMap;
use log::{error, warn};
use packbed::{reader as bed_reader, record::GenePred};
use rayon::prelude::*;

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::str::from_utf8;
use std::sync::Arc;

use crate::{cli::BlastArgs, utils::*};

const ORF_PEP: &str = "orfs.pep.fa";
const ORF_BED: &str = "orfs.bed";
const RESULT: &str = "dmd.result";

pub fn run_blast(args: BlastArgs) {
    let pep = orfipy(&args.common.fasta, &args.common.outdir, &args.orfipy);

    let (dedup, index) = deduplicate(
        &pep,
        true,
        args.orf_min_len,
        args.orf_min_percent,
        "M".as_bytes(),
        SeqType::AminoAcid,
    );

    diamond(&dedup, &args.database, &index, &args.common.alignments);
}

fn orfipy(fasta: &PathBuf, outdir: &PathBuf, executable: &PathBuf) -> PathBuf {
    let dir = outdir.join("orf");

    std::fs::create_dir_all(&dir)
        .unwrap_or_else(|e| panic!("ERROR: could not create directory -> {e}!"));

    let cmd = format!(
        "{} {} --pep {} --bed {} --partial-5 --partial-3 --include-stop --min 100 --ignore-case --outdir {}",
        executable.display(),
        fasta.display(),
        ORF_PEP,
        ORF_BED,
        &dir.display()
    );

    std::process::Command::new("bash")
        .arg("-c")
        .arg(cmd)
        .status()
        .unwrap_or_else(|e| panic!("ERROR: failed to execute orfipy command -> {e}"));

    dir.join(ORF_PEP)
}

fn diamond(dedup: &PathBuf, database: &PathBuf, index: &PathBuf, alignments: &PathBuf) {
    let dmd = dedup.with_extension("diamond");
    let cmd = format!(
        "diamond blastp --query {} --db {} --out {} --outfmt 6 --threads 8 --sensitive -e 1e-10",
        dedup.display(),
        database.display(),
        dmd.display()
    );

    std::process::Command::new("bash")
        .arg("-c")
        .arg(cmd)
        .status()
        .unwrap_or_else(|e| panic!("ERROR: failed to execute diamond command -> {e}"));

    let mut index = read_index(&index);
    let predictions = bed_reader(&dmd)
        .unwrap_or_else(|e| panic!("ERROR: failed to read blast predictions file -> {e}"));

    let mut writer = BufWriter::new(
        File::create(dmd.with_extension(RESULT)).unwrap_or_else(|e| {
            panic!("ERROR: failed to create output file for blast results -> {e}");
        }),
    );

    let accumulator = DashMap::new();

    let records = bed_to_custom_struct_collection::<GenePred>(
        bed_reader(alignments)
            .unwrap_or_else(|e| panic!("ERROR: failed to read BED file -> {e}"))
            .into(),
        config::BedColumn::Name,
        config::BedOperation::SplitName("__", 0), // INFO: R9834_chr16__FC37#TC0#PA0#PR0#IY887
    )
    .unwrap_or_else(|e| panic!("ERROR: failed construct BED to GenePred collection -> {e}"));

    // INFO: filtering repeated blast hits by percent_identity -> preserving best
    predictions
        .par_lines()
        .filter(|line| !line.starts_with('#'))
        .for_each(|line| {
            let parts: Vec<&str> = line.split('\t').collect();

            // INFO: 16 sp|Q9QX47|SON_MOUSE 100 497 0 0 42 538 1089 1585 1.20e-163 515
            let id = parts[0].parse::<u32>().unwrap_or_else(|_| {
                panic!("ERROR: failed to parse ID from line: {}", line);
            });

            let data = BlastRecord::from_parts(&parts);

            // WARN: using a transition collection to retain the best blast record based on % aligned
            accumulator
                .entry(id)
                .and_modify(|existing_data: &mut BlastRecord| {
                    if data.blast_pid > existing_data.blast_pid {
                        *existing_data = data.clone();
                    }
                })
                .or_insert(data);
        });

    // INFO: inflate results!
    accumulator.iter_mut().for_each(|mut record| {
        let (id, data) = record.pair_mut();

        // INFO: unpacking index reference -> queries
        // INFO: for each query all blast records
        // INFO: { index_id : [(read_id: u16, chr_bytes: [u8; chr_len], orf: u16, subseq_orf: u16, seq_len: usize)] }
        // INFO: { 0 : [(read_id: 5903, chr_bytes: [16, 32], orf: 1, subseq_orf: 3, seq_len: 350)] }
        let queries = index.get(id).unwrap_or_else(|| {
            panic!("ERROR: no queries found for ID: {}", id);
        });

        for query in queries.into_iter() {
            let (read_id, chr, orf, subseq_orf, seq_len, start, end) = query;

            let chr = format!("chr{}", from_utf8(&chr).unwrap());
            let cannonical_id = format!("R{}_{}", read_id, chr);
            let query_id = format!("{}.p{}@{}", cannonical_id, orf, subseq_orf);

            // INFO: retrieving the reference gene prediction record
            let (orf_start, orf_end) = records
                .get_mut(&chr)
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: chromosome from {} not found in sequences -> {}!",
                        cannonical_id, chr
                    );
                })
                .get_mut(&cannonical_id)
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: id not found in BED, this is a bug -> {}!",
                        cannonical_id
                    );
                })
                .map_absolute_cds(*start as u64, *end as u64)
                .unwrap_or_default();

            // WARN: skipping unreliable ORFs for the current alignment
            // INFO: none of these will match any other prediction because
            // INFO: the fall off any exonic boundary
            if orf_start == 0 && orf_end == 0 {
                warn!(
                    "WARN: ORF start and end are zero for ID: {}, skipping!",
                    query_id
                );
                continue;
            }

            // INFO: updating blast data with % algined + id
            data.set_percent_aligned((data.blast_alignment_len as f32 / *seq_len as f32) * 100.0);
            data.set_id(query_id);

            writer
                .write_all(
                    format!(
                        "{}\t{}\t{:.2}\t{:e}\t{}\t{}\t{:2}\t{}\t{}\t{}\t{}\n",
                        data.blast_id,
                        data.blast_idx_id,
                        data.blast_pid,
                        data.blast_e_value,
                        data.blast_offset,
                        data.blast_alignment_len,
                        data.percent_aligned,
                        start,
                        end,
                        orf_start,
                        orf_end
                    )
                    .as_bytes(),
                )
                .unwrap_or_else(|e| {
                    panic!("ERROR: failed to write blast record to file -> {e}");
                });
        }

        // INFO: removing the ID from the index to remain with unused IDs
        index.remove(id);
    });

    // INFO: repeating the process for unused ids
    // INFO: add tag DM to the ID -> identify unused ids
    index.iter().for_each(|(id, queries)| {
        for query in queries {
            let (read_id, chr, orf, subseq_orf, _, start, end) = query;

            let chr = format!("chr{}", from_utf8(&chr).unwrap());
            let cannonical_id = format!("R{}_{}", read_id, chr);
            let query_id = format!("{}.p{}@{}#DM", cannonical_id, orf, subseq_orf);

            // INFO: retrieving the reference gene prediction record
            let (orf_start, orf_end) = records
                .get_mut(&chr)
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: chromosome from {} not found in sequences -> {}!",
                        cannonical_id, chr
                    );
                })
                .get_mut(&cannonical_id)
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: id not found in BED, this is a bug -> {}!",
                        cannonical_id
                    );
                })
                .map_absolute_cds(*start as u64, *end as u64)
                .unwrap_or_default();

            // WARN: skipping unreliable ORFs for the current alignment
            // INFO: none of these will match any other prediction because
            // INFO: the fall off any exonic boundary
            if orf_start == 0 && orf_end == 0 {
                warn!(
                    "WARN: ORF start and end are zero for ID: {}, skipping!",
                    query_id
                );
                continue;
            }

            writer
                .write_all(
                    format!(
                        "{}\t{}\t{:.2}\t{:e}\t{}\t{}\t{:2}\t{}\t{}\t{}\t{}\n",
                        query_id, id, 0.0, 1.0, 0, 0, 0.0, start, end, orf_start, orf_end
                    )
                    .as_bytes(),
                )
                .unwrap_or_else(|e| {
                    panic!("ERROR: failed to write unused blast record to file -> {e}");
                });
        }
    });
}

#[derive(Debug, Clone, PartialEq)]
pub struct BlastRecord {
    pub blast_id: String,         // ID of the blast record
    pub blast_idx_id: u32,        // Indexed ID of the blast record
    pub blast_pid: f32,           // Percentage of identical matches
    pub blast_e_value: f64,       // E-value of the match
    pub blast_offset: i32,        // Offset in the query sequence where the match starts
    pub blast_alignment_len: u32, // Length of the alignment
    pub percent_aligned: f32,     // Percentage of the query sequence that is aligned
}

impl BlastRecord {
    pub fn from_parts(parts: &[&str]) -> Self {
        if parts.len() < 10 {
            panic!("ERROR: not enough parts to create BlastRecord -> {parts:?}");
        }

        // INFO: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        // INFO: 16 sp|Q9QX47|SON_MOUSE 100 497 0 0 42 538 1089 1585 1.20e-163 515
        let blast_idx_id = parts[0].parse::<u32>().unwrap_or_else(|_| {
            panic!("ERROR: failed to parse blast ID from parts: {:?}", parts);
        });
        let blast_pid = parts[2].parse::<f32>().unwrap_or_else(|_| {
            panic!("ERROR: failed to parse blast PID from parts: {:?}", parts);
        });
        let blast_e_value = parts[10].parse::<f64>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast E-value from parts: {:?}",
                parts
            );
        });
        let blast_offset = parts[8].parse::<i32>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast offset from parts: {:?}",
                parts
            )
        }) - parts[6].parse::<i32>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast offset from parts: {:?}",
                parts
            );
        });
        let blast_alignment_len = parts[7].parse::<u32>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast offset from parts: {:?}",
                parts
            )
        }) - parts[6].parse::<u32>().unwrap_or_else(|_| {
            panic!(
                "ERROR: failed to parse blast offset from parts: {:?}",
                parts
            )
        }) + 1;

        Self {
            blast_id: String::new(), // INFO: set on the fly
            blast_idx_id,
            blast_pid,
            blast_e_value,
            blast_offset,
            blast_alignment_len,
            percent_aligned: 0.0, // INFO: set on the fly
        }
    }

    pub fn set_percent_aligned(&mut self, percent: f32) {
        self.percent_aligned = percent;
    }

    pub fn set_id(&mut self, id: String) {
        self.blast_id = id;
    }
}

/// Deduplicates a .fa by matching records with repeated sequences,
/// allowing subsequence nesting by any aminoacid ['M' default]
pub fn deduplicate(
    fasta: &PathBuf,
    do_nesting: bool,
    min_len: usize,
    min_percent: f32,
    pattern: &[u8],
    seq_type: SeqType,
) -> (PathBuf, PathBuf) {
    let seqs =
        parse_fa(fasta).unwrap_or_else(|e| panic!("ERROR: failed to parse FASTA file -> {e}"));

    let mut mapper = HashMap::new();
    let mut index = create_index(fasta);
    let mut dedup = create_fasta(fasta, "dedup.fa")
        .unwrap_or_else(|| panic!("ERROR: could not create file {:?}", fasta));

    // INFO: loops through sequences and populates mapper
    for (header, seq) in seqs.iter() {
        let len = seq.len();
        let mut key = Vec::new();
        let header = header.to_string();

        for &b in seq {
            if b != b'\n' {
                key.push(b);
            }
        }

        if do_nesting {
            split_record(
                &header,
                seq,
                len,
                min_len,
                min_percent,
                &mut mapper,
                pattern,
                seq_type,
            )
        }

        let record = Arc::<[u8]>::from(header.into_bytes());
        mapper.entry(key).or_insert(Vec::new()).push(record);
    }

    let _ = make_index(mapper, &mut index, &mut dedup);

    // INFO: cleaning footprint on the fly
    // if !seqs.is_empty() {
    //     std::fs::remove_file(fasta).expect("ERROR: failed to remove original FASTA file");
    // }

    return (
        fasta.with_extension("dedup.fa"),
        fasta.with_extension("dedup.index"),
    );
}

#[derive(Debug, Clone, Copy)]
pub enum SeqType {
    Nucleotide, // search codons like ATG
    AminoAcid,  // search residues like M
}

/// Splits a fasta record based on start signal (ATG or M), respecting codon logic if needed.
fn split_record(
    header: &String,
    seq: &Vec<u8>,
    seq_length: usize,
    min_len: usize,
    min_percent: f32,
    mapper: &mut HashMap<Vec<u8>, Vec<Arc<[u8]>>>,
    needle: &[u8],     // b"ATG" or b"M"
    seq_type: SeqType, // determines scanning logic
) {
    // INFO: always write the original full sequence
    let mut orf_count = 0;

    match seq_type {
        SeqType::Nucleotide => {
            let codon_len = 3;
            let mut pos = codon_len; // INFO: skip first codon
            while pos + codon_len <= seq_length {
                if &seq[pos..pos + codon_len] == needle {
                    let len_remaining = seq_length - pos;
                    let percent = (len_remaining as f32 / 3_f32) / (seq_length as f32 / 3_f32);

                    if (len_remaining as f32 / 3_f32) >= min_len as f32 && percent >= min_percent {
                        orf_count += 1;
                        let sub_seq = &seq[pos..];
                        let sub_id = format!("{}@{}", header, orf_count);

                        let mut inner_key = Vec::with_capacity(sub_seq.len());

                        for &b in sub_seq {
                            if b != b'\n' {
                                inner_key.push(b);
                            }
                        }

                        let record = Arc::<[u8]>::from(sub_id.clone().into_bytes());
                        mapper.entry(inner_key).or_insert(Vec::new()).push(record);
                    }
                }
                pos += codon_len;
            }
        }
        SeqType::AminoAcid => {
            for (pos, &aa) in seq.iter().enumerate().skip(1) {
                if aa == needle[0] {
                    let len_remaining = seq_length - pos;
                    let percent = len_remaining as f32 / seq_length as f32;

                    if len_remaining >= min_len && percent >= min_percent {
                        orf_count += 1;
                        let sub_seq = &seq[pos..];
                        let sub_id = format!("{}@{}", header, orf_count);

                        let mut inner_key = Vec::with_capacity(sub_seq.len());

                        for &b in sub_seq {
                            if b != b'\n' {
                                inner_key.push(b);
                            }
                        }

                        let record = Arc::<[u8]>::from(sub_id.clone().into_bytes());
                        mapper.entry(inner_key).or_insert(Vec::new()).push(record);
                    }
                }
            }
        }
    }
}

fn make_index(
    mapper: HashMap<Vec<u8>, Vec<Arc<[u8]>>>,
    index_writer: &mut BufWriter<File>,
    dedup_writer: &mut BufWriter<File>,
) {
    // INFO: ensuring index is filled up -> will write directly in bytes
    if !mapper.is_empty() {
        let mut count: u32 = 0;

        // INFO: grabs each group and writes 1 sequence per group + an index of all records
        // INFO: pointing to that sequence with the following format:
        // INFO: id_len id seq_len n_headers read_chr [read_ids]
        // INFO: where each [read_id] follows the format: id orf subseq_orf start end
        for (seq, records) in mapper {
            index_writer.write_all(&count.to_be_bytes()).unwrap();

            let seq_len = seq.len() as u32;
            index_writer.write_all(&seq_len.to_be_bytes()).unwrap();

            let n_headers = records.len() as u32;
            index_writer.write_all(&n_headers.to_be_bytes()).unwrap();

            let _ = writeln!(dedup_writer, ">{}\n{}", count, from_utf8(&seq).unwrap());

            // INFO: we do not need to write chr for each record -> assuming same chr for all records
            // INFO: consider the following read names:
            //  >R10589_chr7__FC28#TC0#PA0#PR0#IY896_ORF.6_[4443-4650](+)_type:complete_length:204_frame:1_start:CTG_stop:TAG
            //  >R10589_chr7__FC28#TC0#PA0#PR0#IY896_ORF.6_[4443-4650](+)_type:complete_length:204_frame:1_start:CTG_stop:TAG@1
            // INFO: their indexed names will be -> 10589 7 6 0 4443 4650 and 10589 7 6 1 4443 4650
            let mut chr_count = 0;
            for record in records {
                if let Some((read_id, read_chr, read_orf, read_subseq_orf, start, end)) =
                    get_read_encoding(&record)
                {
                    // INFO: only writing chr for 1st record
                    if chr_count < 1 {
                        index_writer.write_all(&[read_chr.len() as u8]).unwrap();
                        index_writer.write_all(&read_chr).unwrap();
                    }
                    index_writer.write_all(&read_id.to_be_bytes()).unwrap();
                    index_writer.write_all(&read_orf.to_be_bytes()).unwrap();
                    index_writer
                        .write_all(&read_subseq_orf.to_be_bytes())
                        .unwrap();
                    index_writer.write_all(&start.to_be_bytes()).unwrap();
                    index_writer.write_all(&end.to_be_bytes()).unwrap();

                    chr_count += 1;
                } else {
                    panic!("Could not parse header: {:?}", std::str::from_utf8(&record));
                }
            }

            count += 1;
        }
    } else {
        error!("ERROR: No sequences found in the FASTA file!");
        std::process::exit(1);
    }
}

fn get_read_encoding(header: &[u8]) -> Option<(u16, Vec<u8>, u16, u16, u32, u32)> {
    // WARN: cannonical -> R6713_chr16__FC48#TC40#PA0#PR0#IY876
    let header = std::str::from_utf8(header).ok().unwrap_or_else(|| {
        panic!("ERROR: failed to convert header to string");
    });

    // R6456_chr16__FC42#TC48#PA0#PR0#IY907_ORF.247 [44369-44519](+) type:complete length:147 frame:3 start:CTG stop:TAA
    let parts: Vec<&str> = header.split('_').collect();

    if parts.len() < 2 {
        panic!(
            "ERROR: header does not have enough parts to parse: {}",
            &header
        );
    }

    let id = parts[0]
        .strip_prefix('R')?
        .parse::<u16>()
        .ok()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse ID from header: {}", header);
        });
    let chr = parts[1]
        .strip_prefix("chr")
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse chromosome from header: {}", header);
        })
        .as_bytes()
        .to_vec();
    let subseq = parts
        .last()?
        .split('@')
        .nth(1)
        .and_then(|s| s.parse::<u16>().ok())
        .unwrap_or(0);
    let orf = parts
        .get(4)
        .unwrap_or(&"ORF.0")
        .split(" ")
        .next()
        .unwrap_or(&"ORF.0")
        .strip_prefix("ORF.")
        .unwrap_or("0") // WARN: enforcing a default value
        .parse::<u16>()
        .ok()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse ORF from header: {}", header);
        });

    // INFO: ORF.247 [44369-44519](+) type:complete length:147 frame:3 start:CTG stop:TAA
    let coords = parts
        .get(4)
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse coordinates from header: {}", header);
        })
        .split(' ')
        .nth(1)
        .and_then(|s| s.strip_prefix('[')) // 44369-44519](+)
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse coordinates from header: {}", header);
        })
        .split(']')
        .next()
        .unwrap_or_else(|| {
            panic!("ERROR: failed to parse coordinates from header: {}", header);
        })
        .split('-')
        .collect::<Vec<&str>>();

    let start = coords
        .get(0)
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or_else(|| {
            panic!(
                "ERROR: failed to parse start coordinate from header: {}",
                header
            );
        });
    let end = coords
        .get(1)
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or_else(|| {
            panic!(
                "ERROR: failed to parse end coordinate from header: {}",
                header
            );
        });

    return Some((id, chr, orf, subseq, start, end));
}

// [id_bytes: [u8; id_len]]
// [seq_len: u32]
// [n_ids: u32]
// [chr_len: u8]
// [chr_bytes: [u8; chr_len]]
// [read_id: u16]
// ...
//
// { index_id : [(read_id: u16, chr_bytes: [u8; chr_len], orf: u16, subseq_orf: u16, seq_len: usize)] }
// { 0 : [(read_id: 5903, chr_bytes: [16, 32], orf: 1, subseq_orf: 3, seq_len: 350)] }
pub fn read_index(index: &PathBuf) -> HashMap<u32, Vec<(u16, Vec<u8>, u16, u16, usize, u32, u32)>> {
    let mut reader = BufReader::new(
        File::open(index).unwrap_or_else(|e| panic!("ERROR: failed to open index -> {e}")),
    );

    let mut result = HashMap::new();

    loop {
        let mut group_id_buf = [0u8; 4];
        if reader.read_exact(&mut group_id_buf).is_err() {
            break;
        }
        let group_id = u32::from_be_bytes(group_id_buf);

        // 3. Read sequence length
        let mut seq_len_buf = [0u8; 4];
        reader.read_exact(&mut seq_len_buf).unwrap();
        let seq_len = u32::from_be_bytes(seq_len_buf) as usize;

        // 4. Read n_headers
        let mut n_headers_buf = [0u8; 4];
        reader.read_exact(&mut n_headers_buf).unwrap();
        let n_headers = u32::from_be_bytes(n_headers_buf);

        // 5. Read chromosome
        let mut chr_len_buf = [0u8; 1];
        reader.read_exact(&mut chr_len_buf).unwrap();
        let chr_len = chr_len_buf[0] as usize;

        let mut chr_buf = vec![0u8; chr_len];
        reader.read_exact(&mut chr_buf).unwrap();
        let chr = chr_buf;

        // 6. Read n_headers records
        let mut records = Vec::with_capacity(n_headers as usize);
        for _ in 0..n_headers {
            let mut read_id_buf = [0u8; 2];
            reader.read_exact(&mut read_id_buf).unwrap();
            let read_id = u16::from_be_bytes(read_id_buf);

            let mut orf_buf = [0u8; 2];
            reader.read_exact(&mut orf_buf).unwrap();
            let orf = u16::from_be_bytes(orf_buf);

            let mut subseq_buf = [0u8; 2];
            reader.read_exact(&mut subseq_buf).unwrap();
            let subseq = u16::from_be_bytes(subseq_buf);

            let mut start_buf = [0u8; 4];
            reader.read_exact(&mut start_buf).unwrap();
            let start = u32::from_be_bytes(start_buf);

            let mut end_buf = [0u8; 4];
            reader.read_exact(&mut end_buf).unwrap();
            let end = u32::from_be_bytes(end_buf);

            // You now have: group_id, seq_len, n_headers, chr, read_id, orf, subseq, start, end
            records.push((read_id, chr.clone(), orf, subseq, seq_len, start, end));
        }

        result.insert(group_id, records);
    }

    result
}
