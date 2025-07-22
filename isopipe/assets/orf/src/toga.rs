use config::bed_to_struct_collection;
use dashmap::DashMap;
use flate2::read::MultiGzDecoder;
use hashbrown::HashSet;
use log::{error, warn};
use packbed::{reader, record::Bed6};
use rayon::prelude::*;

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use crate::cli::TogaArgs;

const QUERY_ANNOTATION: &str = "query_annotation.bed";
const TRANSCRIPT_METADATA: &str = "meta/transcript_metadata.tsv.gz";
const SELENOCYSTEINE_CODONS: &str = "meta/selenocysteine_codons.tsv";
const TOGA: &str = "toga";
const TOGA_MERGED: &str = "toga_merged.tsv";

/// Runs the TOGA module of the ORF pipeline.
///
/// This function coordinates the processing of various input files to generate
/// a single, merged output file containing complete TOGA predictions. The pipeline
/// involves three main steps:
/// 1. Reads transcript metadata and selenocysteine codon information from gzipped files.
/// 2. Reads query annotation data from a BED file.
/// 3. Merges the metadata and annotation data by matching records based on their ID.
/// 4. Writes the complete `TogaPrediction` objects to a final merged output file.
///
/// # Arguments
///
/// * `args` - A `TogaArgs` struct containing paths to all necessary input files
///            and the desired output directory.
///
/// # Panics
///
/// This function will panic if:
/// - The output directory cannot be created.
/// - Any of the input files (metadata, selenocysteine codons, or query annotations)
///   cannot be read or parsed correctly.
/// - A record ID from the query annotation BED file does not have a corresponding
///   entry in the metadata map, indicating a data mismatch.
/// - The final output file cannot be created or written to.
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
///
/// struct TogaArgs {
///     path: PathBuf,
///     outdir: PathBuf,
/// }
///
/// let args = TogaArgs {
///     path: PathBuf::from("/path/to/toga_data"),
///     outdir: PathBuf::from("/path/to/output"),
/// };
///
/// run_toga(args);
/// ```
pub fn run_toga(args: TogaArgs) {
    let dir = args.outdir.join(TOGA);
    std::fs::create_dir_all(&dir)
        .unwrap_or_else(|e| panic!("ERROR: could not create directory -> {e}!"));

    let map = metadata(
        args.path.join(TRANSCRIPT_METADATA),
        args.path.join(SELENOCYSTEINE_CODONS),
    );

    let records = bed_to_struct_collection::<Bed6>(
        reader(args.path.join(QUERY_ANNOTATION))
            .unwrap_or_else(|e| panic!("ERROR: failed to read BED file -> {e}"))
            .into(),
        config::BedColumn::Name,
    )
    .unwrap_or_else(|e| panic!("ERROR: failed construct BED to GenePred collection -> {e}"));

    // INFO: process each chromosome in parallel
    records.into_par_iter().for_each(|(_, mut rows)| {
        rows.iter_mut().for_each(|(id, data)| {
            map.get_mut(id)
                .unwrap_or_else(|| {
                    panic!("ERROR: no metadata found for ID {id}!");
                })
                .update_rest(data);
        });
    });

    let file = dir.join(TOGA_MERGED);
    let mut writer = BufWriter::new(
        File::create(&file)
            .unwrap_or_else(|e| panic!("ERROR: could not create file {} -> {e}!", file.display())),
    );

    map.into_iter().for_each(|(id, toga)| {
        if let Err(e) = writeln!(writer, "{}", toga) {
            error!("ERROR: could not write to file {} -> {e}!", file.display());
        } else {
            log::debug!(
                "INFO: wrote TogaPrediction for ID {} to file {}",
                id,
                file.display()
            );
        }
    });
}

/// Reads a gzipped file line by line, parses the first four columns (id, label, pid, blosum),
/// and stores them in a `HashMap<String, TogaPrediction>` where the key is the 'id'.
///
/// # Arguments
///
/// * `file` - A `PathBuf` representing the path to the gzipped input file.
///
/// # Returns
///
/// A `HashMap<String, TogaPrediction>` containing the parsed data, keyed by prediction ID.
///
/// # Example
///
/// ```rust, no_run
/// let path = PathBuf::from("/path/to/file.gz");
/// let data = metadata(file).unwrap();
/// ```
fn metadata(metadata: PathBuf, codons: PathBuf) -> DashMap<String, TogaPrediction> {
    let decoder =
        MultiGzDecoder::new(File::open(&metadata).unwrap_or_else(|e| {
            panic!("ERROR: could not open file {} -> {e}!", metadata.display())
        }));
    let reader = BufReader::with_capacity(64 * 1024, decoder); // 64KB buffer
    let map = DashMap::with_capacity(100_000);

    let codons = add_selenocysteine(codons);

    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                eprintln!("Skipping line due to read error: {e}");
                continue;
            }
        };

        // INFO: only require first 4 fields!
        let mut fields = line.split('\t').take(4);
        let id = match fields.next() {
            Some(f) => f.to_string(),
            None => continue,
        };
        let label = match fields.next() {
            Some(f) => f.to_string(),
            None => continue,
        };
        let pid = match fields.next().and_then(|f| f.parse::<f64>().ok()) {
            Some(p) => p,
            None => continue,
        };
        let blosum = match fields.next().and_then(|f| f.parse::<f64>().ok()) {
            Some(b) => b,
            None => continue,
        };

        // Fill only required fields; others are dummies for now
        let prediction = TogaPrediction {
            id: id.clone(),
            label,
            pid,
            blosum,
            masked: codons.contains(&id),
            start: 0,
            end: 0,
            strand: "".to_string(),
            chrom: "".to_string(),
            key: "".to_string(),
        };

        map.insert(id, prediction);
    }

    map
}

/// Represents a Toga prediction with parsed and default fields.
/// Only `id`, `label`, `pid`, and `blosum` are populated from the input file.
/// The other fields are initialized with default values.
#[derive(Debug, Clone)]
struct TogaPrediction {
    id: String,
    label: String,
    pid: f64,
    blosum: f64,
    masked: bool,
    start: u64,
    end: u64,
    strand: String,
    chrom: String,
    key: String,
}

impl TogaPrediction {
    /// Updates the fields of the TogaPrediction instance with the provided data.
    ///
    /// # Arguments
    ///
    /// * `data` - A reference to a `Bed6` record containing additional data to update.
    fn update_rest(&mut self, data: &Bed6) {
        self.start = data.coord.0;
        self.end = data.coord.1;
        self.strand = data.strand.as_str().to_string();
        self.chrom = data.chrom.clone();
        self.key = format!("{}:{}-{}", data.chrom, data.coord.0, data.coord.1);
    }
}

/// Implements the `Display` trait for `TogaPrediction` to format the output
/// in a tab-separated format suitable for writing to a file
impl std::fmt::Display for TogaPrediction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{:.2}\t{:.2}\t{}\t{}\t{}\t{}\t{}",
            self.id,
            self.chrom,
            self.label,
            self.pid,
            self.blosum,
            self.start,
            self.end,
            self.strand,
            self.key,
            self.masked
        )
    }
}

/// Implements the `Default` trait for `TogaPrediction` to easily initialize
/// fields that are not parsed from the input file.
impl Default for TogaPrediction {
    fn default() -> Self {
        TogaPrediction {
            id: String::new(),     // Default empty string
            label: String::new(),  // Default empty string
            pid: 0.0,              // Default 0.0
            blosum: 0.0,           // Default 0.0
            masked: false,         // Default false
            start: 0,              // Default 0
            end: 0,                // Default 0
            strand: String::new(), // Default empty string
            chrom: String::new(),  // Default empty string
            key: String::new(),    // Default empty string
        }
    }
}

/// Reads a file line by line, extracts the first tab-separated column from each line,
/// and stores these unique identifiers in a `HashSet<String>`.
///
/// # Arguments
///
/// * `codons` - A `PathBuf` representing the path to the input file.
///
/// # Returns
///
/// A `HashSet<String>` containing the unique IDs from the first column of the file.
///
/// # Example
///
/// ```rust, no_run
/// use std::path::PathBuf;
/// use std::collections::HashSet;
///
/// let codons_file = PathBuf::from("path/to/codons.txt");
/// let unique_codons: HashSet<String> = add_selenocysteine(codons_file);
/// ```
fn add_selenocysteine(codons: PathBuf) -> HashSet<String> {
    let reader = BufReader::new(
        File::open(&codons)
            .unwrap_or_else(|e| panic!("ERROR: could not open file {} -> {e}!", codons.display())),
    );
    let mut map = HashSet::with_capacity(500);

    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                eprintln!("WARN: skipping line due to read error: {e}");
                continue;
            }
        };

        if line.trim().is_empty() {
            continue;
        }

        // Split the line by tab and take the first column
        if let Some(id) = line.split('\t').next() {
            map.insert(id.to_string());
        }
    }

    map
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;

    #[test]
    fn test_metadata_function() {
        let test_metadata_file_path = PathBuf::from("test_data_metadata.gz");

        {
            let file = File::create(&test_metadata_file_path)
                .expect("Could not create test metadata file");
            let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
            writeln!(encoder, "gene1\tlabelA\t0.95\t0.80\textra1\textra2").unwrap();
            writeln!(encoder, "gene2\tlabelB\t0.88\t0.75\textra3\textra4").unwrap();
            writeln!(encoder, "gene3\tlabelC\tinvalid_pid\t0.70\textra5").unwrap(); // Malformed PID
            writeln!(encoder, "gene4\tlabelD\t0.90\tinvalid_blosum").unwrap(); // Malformed Blosum
            writeln!(encoder, "gene5\tlabelE\t0.85").unwrap(); // Too few columns
            writeln!(encoder, "gene6\tlabelF\t0.92\t0.82").unwrap(); // Valid with exactly 4 columns
            writeln!(encoder, "").unwrap(); // Empty line
            encoder.finish().unwrap();
        }

        let test_codons_file_path = PathBuf::from("test_codons.txt");
        {
            let mut file =
                File::create(&test_codons_file_path).expect("Could not create test codons file");
            writeln!(file, "gene1\tvalA").unwrap();
            writeln!(file, "gene2\tvalB").unwrap();
            writeln!(file, "gene1\tvalC").unwrap(); // Duplicate
            writeln!(file, "gene3").unwrap(); // Single column
            writeln!(file, "").unwrap(); // Empty line
            writeln!(file, "codon4\tvalD\tvalE").unwrap();
        }

        println!(
            "INFO: Processing metadata file: {}",
            test_metadata_file_path.display()
        );
        let result_map = metadata(
            test_metadata_file_path.clone(),
            test_codons_file_path.clone(),
        );

        dbg!(&result_map);

        // Verify specific entries for metadata
        assert_eq!(result_map.len(), 3); // gene1, gene2, gene6 should be parsed
        assert!(result_map.contains_key("gene1"));
        assert!(result_map.contains_key("gene2"));
        assert!(result_map.contains_key("gene6"));
        assert!(!result_map.contains_key("gene3")); // Should be skipped due to parse error
        assert!(!result_map.contains_key("gene4")); // Should be skipped due to parse error
        assert!(!result_map.contains_key("gene5")); // Should be skipped due to too few columns

        // Clean up the metadata test file
        fs::remove_file(&test_metadata_file_path).expect("Could not remove test metadata file");
        fs::remove_file(&test_codons_file_path).expect("Could not remove test codons file");
    }

    #[test]
    fn test_add_selenocysteine_function() {
        let test_codons_file_path = PathBuf::from("test_codons.txt");
        {
            let mut file =
                File::create(&test_codons_file_path).expect("Could not create test codons file");
            writeln!(file, "codon1\tvalA").unwrap();
            writeln!(file, "codon2\tvalB").unwrap();
            writeln!(file, "codon1\tvalC").unwrap(); // Duplicate
            writeln!(file, "codon3").unwrap(); // Single column
            writeln!(file, "").unwrap(); // Empty line
            writeln!(file, "codon4\tvalD\tvalE").unwrap();
        }

        println!(
            "INFO: Processing codons file: {}",
            test_codons_file_path.display()
        );
        let unique_codons = add_selenocysteine(test_codons_file_path.clone());
        dbg!(&unique_codons);

        // INFO: verify specific entries for add_selenocysteine
        assert_eq!(unique_codons.len(), 4); // codon1, codon2, codon3, codon4
        assert!(unique_codons.contains("codon1"));
        assert!(unique_codons.contains("codon2"));
        assert!(unique_codons.contains("codon3"));
        assert!(unique_codons.contains("codon4"));
        assert!(!unique_codons.contains("valA")); // Should not contain values from other columns

        // INFO: clean up the codons test file
        fs::remove_file(&test_codons_file_path).expect("Could not remove test codons file");
    }
}
