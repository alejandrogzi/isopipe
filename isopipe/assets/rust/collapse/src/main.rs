use std::collections::HashMap;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::prelude::*;
use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::RwLock;
use std::{fmt::Debug, io::BufReader};

use collapse::cli::{Args, Command};

use clap::Parser;
use log;
use once_cell::sync::Lazy;
use rayon::prelude::*;
use simple_logger::init_with_level;
use xxhash_rust::xxh3::xxh3_64;

static CHROM_INTERNER: Lazy<RwLock<HashMap<String, u32>>> =
    Lazy::new(|| RwLock::new(HashMap::with_capacity(1024)));
static NEXT_CHROM_ID: AtomicU32 = AtomicU32::new(1);

fn main() {
    let start = std::time::Instant::now();

    let args = Args::parse();
    init_with_level(args.level).unwrap_or_else(|e| log::error!("ERROR: Logger init error -> {e}"));

    match args.command {
        Command::Run(args) => {
            let tracks = unpack(args.bed);
        }
        Command::Read(args) => todo!(),
    }

    let elapsed = start.elapsed();
    log::info!("INFO: Elapsed time: {elapsed:?}");
}

fn unpack<P: AsRef<Path> + Debug + Sync + Send>(
    files: Vec<P>,
) -> Result<HashMap<BinKey, Queue>, Box<dyn std::error::Error>> {
    let contents = par_reader(files)?;
    let tracks = par_parse_tracks(&contents)?;

    Ok(tracks)
}

// collapse run --bed 1.bed,2.bed,3.bed --index [binary or .gz?] --write [writes each collpased file always to /collapsed] --outdir
// collapse read --index <INDEX> --read <READ> [print to stdout] --write [writes utf8]
// collapse run --bed 1.bed --extend [appends what would be --index as an additional column + automatically writes + merges all beds] --outdir

fn parse_tracks(contents: &str) -> Result<HashMap<BinKey, Queue>, Box<dyn std::error::Error>> {
    let mut tracks = HashMap::new();

    contents
        .lines()
        .filter(|row| !row.starts_with("#"))
        .filter_map(|line| Record::parse(line).ok())
        .for_each(|record| {
            // INFO: if record not in tracks, create a new queue
            if !tracks.contains_key(&record.key) {
                tracks.insert(
                    record.key,
                    Queue {
                        reads: vec![],
                        count: 0,
                        rep_line: record.line,
                        header: record.read,
                    },
                );
            } else {
                // INFO: if record already in tracks, increment counts and queue
                let queue = tracks.get_mut(&record.key).unwrap_or_else(|| {
                    panic!("ERROR: record key not found in tracks -> {:?}", record.read)
                });

                queue.count += 1;
                queue.reads.push(record.read);
            }
        });

    Ok(tracks)
}

fn par_parse_tracks(contents: &str) -> Result<HashMap<BinKey, Queue>, Box<dyn std::error::Error>> {
    let tracks = contents
        .par_lines()
        .filter(|row| !row.starts_with("#"))
        .filter_map(|line| Record::parse(line).ok())
        .fold(
            || HashMap::new(),
            |mut acc, record| {
                // INFO: if record not in tracks, create a new queue
                acc.entry(record.key)
                    .and_modify(|queue: &mut Queue| {
                        queue.count += 1;
                        queue.reads.push(record.read.clone());
                    })
                    .or_insert(Queue {
                        reads: vec![],
                        count: 0,
                        rep_line: record.line,
                        header: record.read,
                    });

                acc
            },
        )
        .reduce(
            || HashMap::new(),
            |mut left, right| {
                // INFO: if record not in tracks, create a new queue
                for (key, right_queue) in right {
                    left.entry(key)
                        .and_modify(|left_queue| {
                            left_queue.count += right_queue.count + 1; // INFO: accounts for right header

                            left_queue.reads.extend(right_queue.reads.clone());
                            left_queue.reads.push(right_queue.header.clone()); // INFO: accounts for right header
                        })
                        .or_insert(right_queue);
                }

                left
            },
        );

    Ok(tracks)
}

fn reader<P: AsRef<Path> + Debug>(file: P) -> Result<String, Box<dyn std::error::Error>> {
    let mut file = File::open(file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}

fn par_reader<P: AsRef<Path> + Debug + Sync + Send>(
    files: Vec<P>,
) -> Result<String, Box<dyn std::error::Error>> {
    let contents: Vec<String> = files
        .par_iter()
        .map(|path| reader(path).unwrap_or_else(|e| panic!("Error reading file: {:?}", e)))
        .collect();

    Ok(contents.concat())
}

pub struct Record {
    key: BinKey,
    read: Vec<u8>,
    line: Vec<u8>,
}

impl Record {
    pub fn parse(line: &str) -> Result<Self, Box<dyn std::error::Error>> {
        if line.is_empty() {
            return Err("ERROR: Empty line".into());
        }

        let mut fields = line.split('\t');

        let (
            chrom,
            start,
            end,
            name,
            _,
            strand,
            thick_start,
            thick_end,
            _,
            _,
            block_sizes,
            block_starts,
        ) = (
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Could not pasrse chrom field from -> {}", line)),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Could not pasrse start field from -> {}", line))
                .parse::<u32>()
                .unwrap_or_else(|e| {
                    panic!("ERROR: Could not pasrse start field from -> {}. {e}", line)
                }),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Could not pasrse end field from -> {}", line))
                .parse::<u32>()
                .unwrap_or_else(|e| {
                    panic!("ERROR: Could not pasrse end field from -> {}. {e}", line)
                }),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Could not pasrse name field from -> {}", line)),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Could not pasrse score field from -> {}", line)),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Could not pasrse strand field from -> {}", line))
                .parse::<char>()
                .unwrap_or_else(|e| {
                    panic!("ERROR: Could not pasrse score field from -> {}. {e}", line)
                }) as u8,
            fields
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: Could not pasrse thick_start field from -> {}", line)
                })
                .parse::<u32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: Could not pasrse thick_start field from -> {}. {e}",
                        line
                    )
                }),
            fields
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: Could not pasrse thick_end field from -> {}", line)
                })
                .parse::<u32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: Could not pasrse thick_end field from -> {}. {e}",
                        line
                    )
                }),
            fields.next().unwrap_or_else(|| {
                panic!("ERROR: Could not pasrse item_rgb field from -> {}", line)
            }),
            fields.next().unwrap_or_else(|| {
                panic!("ERROR: Could not pasrse block_count field from -> {}", line)
            }),
            fields
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: Could not pasrse block_sizes field from -> {}", line)
                })
                .split(',')
                .map(|s| {
                    s.parse::<u32>().unwrap_or_else(|e| {
                        panic!(
                            "ERROR: Could not pasrse block_sizes field from -> {}. {e}",
                            line
                        )
                    })
                })
                .collect::<Vec<_>>(),
            fields
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: Could not pasrse block_starts field from -> {}",
                        line
                    )
                })
                .split(',')
                .map(|s| {
                    s.parse::<u32>().unwrap_or_else(|e| {
                        panic!(
                            "ERROR: Could not pasrse block_starts field from -> {}. {e}",
                            line
                        )
                    })
                })
                .collect::<Vec<_>>(),
        );

        let gaps = block_sizes
            .iter()
            .zip(block_starts.iter())
            .map(|(block_size, block_start)| {
                (start + block_start, start + block_size + block_start)
            })
            .collect::<Vec<(u32, u32)>>()
            .windows(2)
            .map(|w| (w[0].1 + 1, w[1].0 - 1))
            .collect::<Vec<(u32, u32)>>();

        Ok(Self {
            key: BinKey::from_parts(chrom, start, end, strand, thick_start, thick_end, gaps),
            read: name.as_bytes().to_vec(),
            line: line.as_bytes().to_vec(),
        })
    }
}

#[derive(Clone, Debug)]
pub struct BinKey {
    /// interned chrom id
    pub chrom_id: u32,
    /// strand as byte: b'+' or b'-'
    pub strand: u8,
    /// transcript half-open start (BED convention)
    pub start: u32,
    /// transcript half-open end
    pub end: u32,
    /// CDS thickStart (BED semantics)
    pub thick_start: u32,
    /// CDS thickEnd (BED semantics)
    pub thick_end: u32,
    /// gaps as vector of (intron_start, intron_end)
    pub gaps: Vec<(u32, u32)>,
    /// cached fingerprint (xxh3_64 of to_bytes_canonical())
    pub fingerprint: u64,
}

impl BinKey {
    pub fn from_parts(
        chrom: &str,
        start: u32,
        end: u32,
        strand: u8,
        thick_start: u32,
        thick_end: u32,
        mut gaps: Vec<(u32, u32)>,
    ) -> Self {
        // INFO: clamp gaps into [start, end]
        // INFO: also drop any intron that is invalid (start >= end or outside transcript)
        gaps.retain(|&(s, e)| s < e && s >= start && e <= end);

        // INFO: static indexing of chromosomes -> {chr: 1, chr2: 2, ...}
        let chrom_id = intern_chromosome(chrom);

        let mut key = BinKey {
            chrom_id,
            strand,
            start,
            end,
            thick_start,
            thick_end,
            gaps,
            fingerprint: 0,
        };

        // INFO: deterministic fingerprint and bytes
        let bytes = key.to_bytes_canonical();
        key.fingerprint = xxh3_64(&bytes);

        key
    }

    pub fn to_bytes_canonical(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(64 + self.gaps.len() * 8);

        bytes.extend_from_slice(&self.chrom_id.to_le_bytes());
        bytes.push(self.strand);
        bytes.extend_from_slice(&self.start.to_le_bytes());
        bytes.extend_from_slice(&self.end.to_le_bytes());
        bytes.extend_from_slice(&self.thick_start.to_le_bytes());
        bytes.extend_from_slice(&self.thick_end.to_le_bytes());

        let n_gaps = self.gaps.len() as u16;
        bytes.extend_from_slice(&n_gaps.to_le_bytes());

        for (s, e) in &self.gaps {
            bytes.extend_from_slice(&s.to_le_bytes());

            let len = e - s;
            bytes.extend_from_slice(&len.to_le_bytes());
        }

        bytes
    }

    pub fn from_bytes_canonical(bytes: &[u8]) -> Self {
        let buf = bytes;
        let mut bytes = bytes;

        let chrom_id = u32::from_le_bytes(bytes[..4].try_into().unwrap());
        bytes = &bytes[4..];

        let strand = bytes[0];
        bytes = &bytes[1..];

        let start = u32::from_le_bytes(bytes[..4].try_into().unwrap());
        bytes = &bytes[4..];

        let end = u32::from_le_bytes(bytes[..4].try_into().unwrap());
        bytes = &bytes[4..];

        let thick_start = u32::from_le_bytes(bytes[..4].try_into().unwrap());
        bytes = &bytes[4..];

        let thick_end = u32::from_le_bytes(bytes[..4].try_into().unwrap());
        bytes = &bytes[4..];

        let n_gaps = u16::from_le_bytes(bytes[..2].try_into().unwrap());
        bytes = &bytes[2..];

        let mut gaps = Vec::with_capacity(n_gaps as usize);
        for _ in 0..n_gaps {
            let s = u32::from_le_bytes(bytes[..4].try_into().unwrap());
            bytes = &bytes[4..];

            let len = u32::from_le_bytes(bytes[..4].try_into().unwrap());
            bytes = &bytes[4..];

            gaps.push((s, s + len));
        }

        Self {
            chrom_id,
            strand,
            start,
            end,
            thick_start,
            thick_end,
            gaps,
            fingerprint: xxh3_64(buf),
        }
    }
}

impl PartialEq for BinKey {
    fn eq(&self, other: &Self) -> bool {
        // INFO: fast path, different fingerprints -> not equal
        if self.fingerprint != other.fingerprint {
            return false;
        }

        // INFO: verify by comparing canonical bytes (defends against hash collisions)
        self.to_bytes_canonical() == other.to_bytes_canonical()
    }
}

impl Eq for BinKey {}

impl Hash for BinKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // INFO: we hash the fingerprint (u64). This keeps hashing cheap and consistent with Eq
        state.write_u64(self.fingerprint);
    }
}

/// Intern a chromosome name to a small u32 id (thread-safe).
fn intern_chromosome(chrom: &str) -> u32 {
    // INFO: Fast path (read lock)
    if let Ok(read) = CHROM_INTERNER.read() {
        if let Some(&id) = read.get(chrom) {
            return id;
        }
    }

    // INFO: Slow path (write lock)
    let mut write = CHROM_INTERNER
        .write()
        .unwrap_or_else(|e| panic!("ERROR: poisoned lock: {e} -> {chrom} - {CHROM_INTERNER:?}"));

    if let Some(&id) = write.get(chrom) {
        return id;
    }

    let id = NEXT_CHROM_ID.fetch_add(1, Ordering::Relaxed);
    write.insert(chrom.to_string(), id);
    id
}

#[derive(Debug, Clone)]
pub struct Queue {
    pub reads: Vec<Vec<u8>>,
    pub count: u32,
    pub rep_line: Vec<u8>,
    pub header: Vec<u8>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_binkey_equality() {
        let key = BinKey::from_parts(
            "chr1",
            50,
            500,
            b'+',
            150,
            250,
            vec![(100, 200), (300, 400)],
        );

        assert_eq!(key.chrom_id, 1);
        assert_eq!(key.strand, b'+');
        assert_eq!(key.start, 50);
        assert_eq!(key.end, 500);
        assert_eq!(key.thick_start, 150);
        assert_eq!(key.thick_end, 250);
        assert_eq!(key.gaps, vec![(100, 200), (300, 400)]);

        dbg!(&key.fingerprint);

        let key2 = BinKey::from_parts(
            "chr1",
            50,
            500,
            b'+',
            150,
            250,
            vec![(100, 200), (300, 400)],
        );

        assert_eq!(key.fingerprint, key2.fingerprint);
        assert_eq!(key.to_bytes_canonical(), key2.to_bytes_canonical());
        assert_eq!(key, key2);
    }

    #[test]
    fn test_binkey_to_bytes() {
        let key = BinKey::from_parts(
            "chr1",
            100,
            200,
            b'+',
            150,
            250,
            vec![(100, 200), (300, 400)],
        );

        let bytes = key.to_bytes_canonical();
        assert_eq!(bytes.len(), 31);
    }

    #[test]
    fn test_binkey_from_bytes() {
        let key = BinKey::from_parts(
            "chr1",
            50,
            500,
            b'+',
            150,
            250,
            vec![(100, 200), (300, 400)],
        );

        let bytes = key.to_bytes_canonical();
        let key2 = BinKey::from_bytes_canonical(&bytes);

        assert_eq!(key, key2);
    }
}
