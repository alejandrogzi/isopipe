use config::{BedColumn, BedColumnValue, bed_to_nested_collection};
use dashmap::DashMap;
use hashbrown::HashMap;
use log::{error, warn};
use memchr::memchr;
use memmap2::Mmap;
use packbed::reader as bed_reader;
use packbed::record::Bed6;
use smol_str::SmolStr;

use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::str::from_utf8;
use std::sync::Arc;

pub fn parse_fa<F: AsRef<Path>>(
    f: F,
) -> Result<HashMap<SmolStr, Vec<u8>>, Box<dyn std::error::Error>> {
    let file = File::open(f).unwrap();
    let mmap = unsafe { Mmap::map(&file).unwrap() };
    let data = mmap.as_ref();

    let mut acc = HashMap::new();
    let mut pos = 0;

    while let Some(start) = memchr(b'>', &data[pos..]) {
        let start = pos + start;
        let end = memchr(b'>', &data[start + 1..]).map_or(data.len(), |e| start + 1 + e);
        let entry = &data[start + 1..end];
        let header_end = memchr(b'\n', entry).unwrap();
        let header = from_utf8(&entry[..header_end]).unwrap().trim();
        let record = &entry[header_end + 1..];
        let seq = record
            .iter()
            .filter(|&&b| b != b'\n' && b != b'\r') // Remove newlines and carriage returns
            .cloned()
            .collect::<Vec<u8>>();

        acc.insert(SmolStr::new(header), seq);
        pos = end;
    }

    Ok(acc)
}

/// Creates an .index file for a given .fa
pub fn create_index(fasta: &PathBuf) -> BufWriter<File> {
    let index = fasta.with_extension("dedup.index");
    let writer = BufWriter::new(
        File::create(&index)
            .unwrap_or_else(|e| panic!("ERROR: cannot create index from sequences -> {e}")),
    );

    writer
}

pub fn parse_bed<K: config::BedParser>(
    bed: &PathBuf,
    columns: Vec<BedColumn>,
) -> DashMap<String, HashMap<String, Vec<BedColumnValue>>> {
    let rows = Arc::from(
        bed_reader(bed).unwrap_or_else(|e| panic!("ERROR: failed to read BED file -> {e}")),
    );
    let records = bed_to_nested_collection::<K>(rows, BedColumn::Name, columns)
        .unwrap_or_else(|e| panic!("ERROR: failed to convert BED to nested collection -> {e}"));

    records
}

/// Creates a .fa file for a given new extension
pub fn create_fasta(fasta: &PathBuf, extension: &str) -> Option<BufWriter<File>> {
    let file = fasta.with_extension(extension);
    if !file.exists() {
        let writer = BufWriter::new(
            File::create(&file).unwrap_or_else(|e| panic!("ERROR: cannot create file -> {e}")),
        );

        Some(writer)
    } else {
        warn!("WARN: file already exists -> {:?}!", file.display());
        None
    }
}
