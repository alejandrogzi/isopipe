// Copyright (c) 2026 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use bigtools::{utils::reopen::Reopen, BigWigRead};
use dashmap::{DashMap, DashSet};
use log::info;
use rayon::prelude::*;

use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::{Arc, Mutex};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
};

use crate::utils::{build_coordinate_key, reverse_complement_in_place, uppercase_in_place};

/// SpliceAI acceptor BigWig filename for minus strand.
pub const ACCEPTOR_MINUS: &str = "spliceAiAcceptorMinus.bw";
/// SpliceAI acceptor BigWig filename for plus strand.
pub const ACCEPTOR_PLUS: &str = "spliceAiAcceptorPlus.bw";
/// SpliceAI donor BigWig filename for minus strand.
pub const DONOR_MINUS: &str = "spliceAiDonorMinus.bw";
/// SpliceAI donor BigWig filename for plus strand.
pub const DONOR_PLUS: &str = "spliceAiDonorPlus.bw";

/// Type alias: (plus_strand_map, minus_strand_map) for both strands.
pub type SpliceMap = (StrandSpliceMap, StrandSpliceMap);
/// Thread-safe map: chromosome -> set of splice site records.
pub type StrandSpliceMap = DashMap<String, DashSet<Vec<u8>>>;
/// Shared splice scores (unused in current implementation).
pub type SharedSpliceMap = (Option<DashMap<usize, f32>>, Option<DashMap<usize, f32>>);
/// Splice scores: vectors of strand maps for each BigWig type.
pub type SpliceScores = (Vec<StrandSpliceMap>, Vec<StrandSpliceMap>);

/// Fetches and processes splice scores from BigWig files.
///
/// This is a convenience function that wraps `make_splice_map`, handling the case where
/// no splice score files are provided. If files are provided, it calls `make_splice_map` to
/// load the data; otherwise, it returns empty maps.
///
/// # Arguments
///
/// * `splice_scores`: An `Option<T>` with the path to the directory containing splice score BigWigs.
/// * `chrs`: A `Vec<String>` of chromosome names to process.
///
/// # Returns
///
/// * A `SpliceScores` type, which is a tuple of two vectors of `StrandSpliceMap`.
///
/// # Example
///
/// ```rust,ignore
/// let scores = get_splice_scores(splice_scores, chrs);
/// ```
pub fn get_splice_scores<T: AsRef<std::path::Path> + std::fmt::Debug>(
    bigwigs: T,
    chrs: Vec<Vec<u8>>,
    genome: HashMap<Vec<u8>, Vec<u8>>,
) -> (StrandSpliceMap, StrandSpliceMap) {
    // INFO: DashMap<String, DashSet<Vec<u8>>> -> chr -> [ b'pos -> score' ]
    make_splice_map(bigwigs, chrs, genome)
}

/// Creates `StrandSpliceMap`s for both plus and minus strands by parsing BigWig files.
///
/// This function takes a directory containing BigWig files for donor and acceptor splice scores for both
/// strands. It then uses `rayon` to parallelize the parsing of these files into `DashMap`s,
/// which are a thread-safe hash map, and returns the results.
///
/// # Arguments
///
/// * `dir`: The path to the directory containing the BigWig files.
/// * `chrs`: A `Vec<String>` of chromosome names to process.
///
/// # Returns
///
/// * A tuple of two `Vec<StrandSpliceMap>`, where the first vector is for the plus strand
///   (donor and acceptor) and the second is for the minus strand.
///
/// # Example
///
/// ```rust,ignore
/// let (plus_scores, minus_scores) = make_splice_map(dir, chrs);
/// ```
pub fn make_splice_map<T: AsRef<std::path::Path> + std::fmt::Debug>(
    dir: T,
    chrs: Vec<Vec<u8>>,
    genome: HashMap<Vec<u8>, Vec<u8>>,
) -> (StrandSpliceMap, StrandSpliceMap) {
    let plus = vec![
        dir.as_ref().join(DONOR_PLUS),
        dir.as_ref().join(ACCEPTOR_PLUS),
    ];
    let minus = vec![
        dir.as_ref().join(DONOR_MINUS),
        dir.as_ref().join(ACCEPTOR_MINUS),
    ];

    info!("Parsing BigWigs...");
    let (plus, minus) = rayon::join(
        || bigwig_to_map(plus, &chrs, Strand::Forward, &genome),
        || bigwig_to_map(minus, &chrs, Strand::Reverse, &genome),
    );

    let [plus_donor, plus_acceptor] = <[_; 2]>::try_from(plus).unwrap();
    let [minus_donor, minus_acceptor] = <[_; 2]>::try_from(minus).unwrap();

    let donor_scores = merge_splice_maps(plus_donor, minus_donor);
    let acceptor_scores = merge_splice_maps(plus_acceptor, minus_acceptor);

    (donor_scores, acceptor_scores)
}

/// Merges two strand-specific splice maps (e.g., plus and minus strands).
///
/// # Arguments
/// * `primary` - First map (e.g., plus strand)
/// * `secondary` - Second map to merge in (e.g., minus strand)
///
/// # Returns
/// Merged map with all entries
///
/// # Example
/// ```rust,ignore
/// let merged = merge_splice_maps(plus_donor, minus_donor);
/// ```
fn merge_splice_maps(primary: StrandSpliceMap, secondary: StrandSpliceMap) -> StrandSpliceMap {
    let merged = primary;

    secondary.into_iter().for_each(|(chr, splice_sites)| {
        if let Some(existing) = merged.get_mut(&chr) {
            splice_sites.into_iter().for_each(|splice_site| {
                existing.insert(splice_site);
            });
        } else {
            merged.insert(chr, splice_sites);
        }
    });

    merged
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Forward,
    Reverse,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
        }
    }
}

// public enums
/// Splice site type
///
/// This enum is used to store the type of splice site.
///
/// # Example
///
/// ```rust, no_run
/// use bwtoms::spliceai::SpliceSite;
///
/// let donor = SpliceSite::Donor;
/// let acceptor = SpliceSite::Acceptor;
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpliceSite {
    Donor,
    Acceptor,
}

impl std::fmt::Display for SpliceSite {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpliceSite::Donor => write!(f, "D"),
            SpliceSite::Acceptor => write!(f, "A"),
        }
    }
}

/// Converts a vector of BigWig files into a vector of thread-safe maps.
///
/// This function is designed to be run in parallel for plus and minus strands. It iterates through
/// the BigWig files, reads the chromosome data, and populates `DashMap`s with scores that
/// meet a certain significance threshold.
///
/// # Arguments
///
/// * `bigwigs`: A `Vec` of paths to the BigWig files (e.g., donor and acceptor).
/// * `chrs`: A slice of `String` representing the chromosomes to be processed.
///
/// # Returns
///
/// * A `Vec<DashMap<String, DashMap<usize, f32>>>` where the outer vector corresponds to donor/acceptor
///   sites, the middle map keys are chromosome names, and the inner map keys are genomic positions with their scores.
///
/// # Example
///
/// ```rust,ignore
/// let splice_maps = bigwig_to_map(bigwigs, &chrs);
/// ```
fn bigwig_to_map<T: AsRef<std::path::Path> + std::fmt::Debug + Sized + Sync>(
    bigwigs: Vec<T>,
    chrs: &[Vec<u8>],
    strand: Strand,
    genome: &HashMap<Vec<u8>, Vec<u8>>,
) -> Vec<DashMap<String, DashSet<Vec<u8>>>> {
    let total_count = AtomicU32::new(0);
    let rs = Mutex::new(vec![DashMap::new(), DashMap::new()]);

    // [donor, acceptor]
    bigwigs
        .into_par_iter()
        .zip(vec![SpliceSite::Donor, SpliceSite::Acceptor])
        .for_each(|(bigwig, splice_site)| {
            let acc = DashMap::new();

            let bwread = BigWigRead::open_file(bigwig).expect("ERROR: Cannot open BigWig file");
            let chroms: Vec<_> = bwread.chroms().to_vec();
            let splice_site_arc = Arc::new(splice_site);

            chroms.into_par_iter().for_each(|chr| {
                // INFO: per-chromosome map
                let mapper = DashSet::new();
                let local_count = AtomicU32::new(0);

                let mut bwread =
                    BigWigRead::reopen(&bwread).expect("ERROR: Cannot re-open BigWig file");

                if !chrs.contains(&chr.name.as_bytes().to_vec()) {
                    return; // INFO: skip chromosomes not in records
                }

                let name = chr.name.clone();
                let length = chr.length;
                let values = bwread
                    .values(&name, 0, length)
                    .expect("ERROR: Cannot read values from BigWig!");

                values.into_iter().enumerate().for_each(|(i, v)| {
                    if v >= crate::cli::SPLICE_AI_SCORE_RECOVERY_THRESHOLD {
                        let pos = i;
                        let sequence = genome
                            .get(name.as_bytes())
                            .unwrap_or_else(|| panic!("ERROR: cannot chr {} in genome!", name));
                        let Some(dnt) = extract_dinucleotide(sequence, pos, splice_site, strand)
                        else {
                            log::warn!(
                                "Skipping SpliceAI site outside chromosome bounds: {}:{}",
                                name,
                                pos
                            );
                            return;
                        };

                        let line =
                            Minisplice::new(&name, pos, strand, splice_site_arc.clone(), v, dnt);
                        mapper.insert(line.as_bytes());
                        local_count.fetch_add(1, Ordering::Relaxed);
                    }
                });

                acc.insert(name, mapper);
                total_count.fetch_add(local_count.load(Ordering::Relaxed), Ordering::Relaxed);
            });

            let mut guard = rs.lock().expect("ERROR: Cannot lock mutex");
            match splice_site {
                SpliceSite::Donor => guard[0] = acc,
                SpliceSite::Acceptor => guard[1] = acc,
            }
        });

    info!(
        "Parsed and combined {} significant splicing scores from BigWigs!",
        total_count.load(Ordering::Relaxed)
    );

    rs.into_inner()
        .expect("ERROR: Cannot unwrap collection of SpliceAI scores!")
}

/// Extracts the dinucleotide at a splice site from the genome sequence.
///
/// Position depends on splice site type and strand.
///
/// # Arguments
/// * `sequence` - Chromosome sequence
/// * `pos` - Genomic position
/// * `splice_site` - Donor or Acceptor
/// * `strand` - Forward or Reverse
///
/// # Returns
/// Dinucleotide bytes or None if out of bounds
///
/// # Example
/// ```rust,ignore
/// let dnt = extract_dinucleotide(&seq, 100, SpliceSite::Donor, Strand::Forward);
/// ```
fn extract_dinucleotide(
    sequence: &[u8],
    pos: usize,
    splice_site: SpliceSite,
    strand: Strand,
) -> Option<Vec<u8>> {
    let mut dinucleotide = match (splice_site, strand) {
        (SpliceSite::Donor, Strand::Forward) => {
            let end = pos.checked_add(2)?;
            sequence.get(pos..end)?.to_vec()
        }
        (SpliceSite::Donor, Strand::Reverse) => {
            let start = pos.checked_sub(1)?;
            let end = pos.checked_add(1)?;
            sequence.get(start..end)?.to_vec()
        }
        (SpliceSite::Acceptor, Strand::Forward) => {
            let start = pos.checked_sub(1)?;
            let end = pos.checked_add(1)?;
            sequence.get(start..end)?.to_vec()
        }
        (SpliceSite::Acceptor, Strand::Reverse) => {
            let end = pos.checked_add(2)?;
            sequence.get(pos..end)?.to_vec()
        }
    };

    match strand {
        Strand::Forward => {}
        Strand::Reverse => {
            reverse_complement_in_place(&mut dinucleotide);
        }
    }

    uppercase_in_place(&mut dinucleotide);

    Some(dinucleotide)
}

/// Minimal splice record for SpliceAI data.
///
/// Format: chr:pos(strand)\tD/A\tscore\tdinucleotide
///
/// # Fields
/// - `chr`: Chromosome name
/// - `position`: Genomic position
/// - `strand`: Forward or Reverse
/// - `splice_site`: Donor or Acceptor
/// - `score`: SpliceAI score (0-1)
/// - `dinucleotide`: Splice site dinucleotide
#[derive(Debug, Clone, PartialEq)]
pub struct Minisplice {
    pub chr: String,
    pub position: usize,
    pub strand: Strand,
    pub splice_site: Arc<SpliceSite>,
    pub score: f32,
    pub dinucleotide: Vec<u8>,
}

impl Minisplice {
    /// Creates a new Minisplice record.
    ///
    /// # Arguments
    /// * `chr` - Chromosome name
    /// * `position` - Genomic position
    /// * `strand` - Forward or Reverse
    /// * `splice_site` - Donor or Acceptor (as Arc)
    /// * `score` - SpliceAI score
    /// * `dinucleotide` - Splice site dinucleotide
    ///
    /// # Example
    /// ```rust,ignore
    /// let ms = Minisplice::new("chr1", 100, Strand::Forward, Arc::new(SpliceSite::Donor), 0.95, b"GT".to_vec());
    /// ```
    pub fn new(
        chr: &str,
        position: usize,
        strand: Strand,
        splice_site: Arc<SpliceSite>,
        score: f32,
        dinucleotide: Vec<u8>,
    ) -> Self {
        Self {
            chr: chr.to_string(),
            position,
            strand,
            splice_site,
            score,
            dinucleotide,
        }
    }

    /// Converts to tab-separated byte format.
    ///
    /// # Example
    /// ```rust,ignore
    /// let bytes = minisplice.as_bytes();
    /// ```
    pub fn as_bytes(&self) -> Vec<u8> {
        format!(
            "{}:{}({})\t{}\t{}\t{}",
            self.chr,
            self.position,
            self.strand,
            self.splice_site,
            self.score,
            std::str::from_utf8(&self.dinucleotide).unwrap()
        )
        .into_bytes()
    }
}

impl std::fmt::Display for Minisplice {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.chr,
            self.position,
            self.strand,
            self.splice_site,
            self.score,
            std::str::from_utf8(&self.dinucleotide).unwrap()
        )
    }
}

/// Parsed SpliceAI record with extracted fields.
///
/// # Fields
/// - `coordinate_key`: Full coordinate string "chr:pos(strand)"
/// - `chr`: Chromosome bytes
/// - `position`: Genomic position
/// - `strand`: Forward or Reverse
/// - `splice_site`: Donor or Acceptor
/// - `score`: SpliceAI score (0-1)
/// - `dinucleotide`: Splice site dinucleotide
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedSpliceAiRecord {
    pub coordinate_key: Vec<u8>,
    pub chr: Vec<u8>,
    pub position: usize,
    pub strand: Strand,
    pub splice_site: SpliceSite,
    pub score: f32,
    pub dinucleotide: Vec<u8>,
}

/// Parses a tab-separated SpliceAI record line.
///
/// # Arguments
/// * `line` - Tab-separated record: coord\tD/A\tscore\tdinucleotide
///
/// # Returns
/// Parsed record or None if invalid
///
/// # Example
/// ```rust,ignore
/// let parsed = parse_spliceai_record(b"chr1:100(+)\tD\t0.95\tGT")?;
/// ```
pub fn parse_spliceai_record(line: &[u8]) -> Option<ParsedSpliceAiRecord> {
    let mut fields = line.split(|byte| *byte == b'\t');
    let coordinate_field = fields.next()?;
    let splice_site_field = fields.next()?;
    let score = std::str::from_utf8(fields.next()?)
        .ok()?
        .parse::<f32>()
        .ok()?;
    let mut dinucleotide = fields.next()?.to_vec();

    if fields.next().is_some() {
        return None;
    }

    uppercase_in_place(&mut dinucleotide);

    let (chr, position, strand) = parse_coordinate_field(coordinate_field)?;
    let splice_site = parse_splice_site(splice_site_field)?;
    let coordinate_key = build_coordinate_key(&chr, position, strand);

    Some(ParsedSpliceAiRecord {
        coordinate_key,
        chr,
        position,
        strand,
        splice_site,
        score,
        dinucleotide,
    })
}

/// Parses coordinate field "chr:position(strand)" into components.
///
/// # Example
/// ```rust,ignore
/// let (chr, pos, strand) = parse_coordinate_field(b"chr1:12345(+)")?;
/// ```
fn parse_coordinate_field(field: &[u8]) -> Option<(Vec<u8>, usize, Strand)> {
    let strand_open = field.iter().rposition(|byte| *byte == b'(')?;
    let strand_close = field.last().copied()?;
    if strand_close != b')' || strand_open == 0 || strand_open + 2 != field.len() - 1 {
        return None;
    }

    let strand = match field.get(strand_open + 1)? {
        b'+' => Strand::Forward,
        b'-' => Strand::Reverse,
        _ => return None,
    };
    let chr_pos = field.get(..strand_open)?;
    let separator = chr_pos.iter().rposition(|byte| *byte == b':')?;
    let chr = chr_pos.get(..separator)?.to_vec();
    let position = std::str::from_utf8(chr_pos.get(separator + 1..)?)
        .ok()?
        .parse::<usize>()
        .ok()?;

    Some((chr, position, strand))
}

/// Parses splice site type field (D = Donor, A = Acceptor).
///
/// # Example
/// ```rust,ignore
/// let site = parse_splice_site(b"D")?; // Some(SpliceSite::Donor)
/// ```
fn parse_splice_site(field: &[u8]) -> Option<SpliceSite> {
    match field {
        b"D" => Some(SpliceSite::Donor),
        b"A" => Some(SpliceSite::Acceptor),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_spliceai_record_extracts_expected_fields() {
        let parsed = parse_spliceai_record(b"chr1:42(-)\tA\t0.75\tag").unwrap();

        assert_eq!(parsed.coordinate_key, b"chr1:42(-)".to_vec());
        assert_eq!(parsed.chr, b"chr1".to_vec());
        assert_eq!(parsed.position, 42);
        assert_eq!(parsed.strand, Strand::Reverse);
        assert_eq!(parsed.splice_site, SpliceSite::Acceptor);
        assert_eq!(parsed.score, 0.75);
        assert_eq!(parsed.dinucleotide, b"AG".to_vec());
    }
}

/// Writes the results to files.
///
/// This function takes two vectors of `StrandSpliceMap`s and writes the results to files.
/// It first writes the plus strand scores to a file named `prefix.plus.tsv`, and then
/// writes the minus strand scores to a file named `prefix.minus.tsv`.
///
/// # Arguments
///
/// * `plus_scores`: A `Vec<StrandSpliceMap>` containing the plus strand scores.
/// * `minus_scores`: A `Vec<StrandSpliceMap>` containing the minus strand scores.
/// * `prefix`: A `String` representing the prefix for the output files.
///
/// # Example
///
/// ```rust,ignore
/// write_results(plus_scores, minus_scores, prefix);
/// ```
pub fn write_results(
    plus_scores: Vec<StrandSpliceMap>,
    minus_scores: Vec<StrandSpliceMap>,
    prefix: String,
) {
    let mut plus_file = BufWriter::new(
        File::create(format!("{}.plus.tsv", prefix))
            .unwrap_or_else(|e| panic!("ERROR: Cannot create plus.tsv file -> {e}!")),
    );
    let mut minus_file = BufWriter::new(
        File::create(format!("{}.minus.tsv", prefix))
            .unwrap_or_else(|e| panic!("ERROR: Cannot create minus.tsv file -> {e}!")),
    );

    plus_scores.iter().for_each(|entry| {
        entry.iter().for_each(|val| {
            let (chr, scores) = val.pair();
            info!(
                "Writing {} scores for {} in forward strand...",
                scores.len(),
                chr
            );
            scores.iter().for_each(|score| {
                plus_file
                    .write_all(&score)
                    .unwrap_or_else(|e| panic!("ERROR: Cannot write to plus.tsv file -> {e}!"));
                plus_file
                    .write_all(b"\n")
                    .unwrap_or_else(|e| panic!("ERROR: Cannot write to plus.tsv file -> {e}!"));
            });
        });
    });
    plus_file.flush().unwrap();

    minus_scores.iter().for_each(|entry| {
        entry.iter().for_each(|val| {
            let (chr, scores) = val.pair();
            info!(
                "Writing {} scores for {} in reverse strand...",
                scores.len(),
                chr
            );
            scores.iter().for_each(|score| {
                minus_file
                    .write_all(&score)
                    .unwrap_or_else(|e| panic!("ERROR: Cannot write to minus.tsv file -> {e}!"));
                minus_file
                    .write_all(b"\n")
                    .unwrap_or_else(|e| panic!("ERROR: Cannot write to minus.tsv file -> {e}!"));
            });
        });
    });
    minus_file.flush().unwrap();
}
