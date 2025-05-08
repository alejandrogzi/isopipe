import logging
from copy import deepcopy
from Bio import SeqIO
from pandas import DataFrame
from bedutils.bed import BedRow

logger = logging.getLogger(__name__)


def multiplex_dataset(input_fasta, output_fasta, output_index, silent=False):
    """
    Maps common sequences in a faste file to a key and writes the relation to an index file.
    :param input_fasta: Input fasta file
    :param output_fasta: Path to output fasta file
    :param output_index: Peth to tab-formatted output index file
    :param silent: Suppress logging
    :return:
    """
    multiplex_dict = dict()
    raw_count = 0

    # Simple FASTA multiplexing procedure:
    # Using sequences as key, collapse record with equal sequences onto one another
    if not silent: logging.info("Reading file")
    with open(input_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if multiplex_dict.get(record.seq):
                targets = multiplex_dict.get(record.seq)
                targets.append(record.id)
                multiplex_dict[record.seq] = targets
            else:
                multiplex_dict[record.seq] = [record.id]
            raw_count += 1

    collapse_ratio = len(multiplex_dict.keys()) / raw_count

    if not silent: logging.info(f"Reduced number of individual sequences by {100 - round(collapse_ratio * 100, 2)}%")

    with open(output_fasta, "w") as fasta:
        with open(output_index, "w", encoding='utf-8') as index:
            for i, key in enumerate(multiplex_dict.keys()):
                index.write(f"sequence_{i}\t{'⍼'.join(multiplex_dict[key])}\n")
                fasta.write(f">sequence_{i}\n{key}\n")

    if not silent:
        logging.info(f"{output_index} \t FASTA: {output_fasta}")


def chunk_bed(input_path, chunk_size, output_path, prefix) -> int:
    """
    Takes an input bed file from disk and splits it into multiple files in output_path/{prefix}_{index}.
    :param input_path: Path to input bed file.
    :param chunk_size: Number of lines per chunk.
    :param output_path: Directory to write chunked files.
    :param prefix: Prefix for chunked files.
    :return: number of chunks produced
    """
    chunk_index = 0
    chunk_buffer = []
    with open(input_path, 'r') as f:
        for i, line in enumerate(f):
            chunk_buffer.append(line)
            if i % chunk_size == 0:
                with open(f"{output_path}/{prefix}_{chunk_index}.bed", 'w') as output_file:
                    output_file.writelines(chunk_buffer)
                    chunk_buffer = []
                    chunk_index += 1

    return chunk_index


def chunk_fasta(input_path, chunk_size, output_path, prefix) -> int:
    """
    Takes an input bed file from disk and splits it into multiple files in output_path/{prefix}_{index}.
    :param input_path: Path to input bed file.
    :param chunk_size: Number of lines per chunk.
    :param output_path: Directory to write chunked files.
    :param prefix: Prefix for chunked files.
    :return: number of chunks produced
    """
    chunk_index = 0
    i = 0
    chunk_buffer = []
    with open(input_path) as f:
        for record in SeqIO.parse(f, "fasta"):
            i += 1
            chunk_buffer.append(record)
            if i % chunk_size == 0:
                with open(f"{output_path}/{prefix}_{chunk_index}.pep", 'w') as output_file:
                    SeqIO.write(chunk_buffer, output_file, "fasta")
                    chunk_buffer = []
                    chunk_index += 1

    if len(chunk_buffer) > 0:
        with open(f"{output_path}/{prefix}_{chunk_index}.pep", 'w') as output_file:
            SeqIO.write(chunk_buffer, output_file, "fasta")
            chunk_index += 1
    return chunk_index


def export_df_orfs(df, input_bed, output_path):
    # populate BED dict
    transcript_dict = dict()
    with open(input_bed, 'r') as f:
        for line in f:
            curr_row = BedRow(line)
            transcript_dict[curr_row.id_str] = curr_row

    dropped = 0
    x = 0
    with open(output_path, 'w', encoding='utf-8') as f:
        for entry in df.itertuples():
            orf_coords = entry.genomic_coords.split("|")
            source_row = deepcopy(transcript_dict.get(entry.canonical_id))

            # Trim BED row to constraints of ORF
            try:
                if orf_coords[-1] == "+":
                    source_row.trim_to_genomic_position_upstream(int(orf_coords[1]))
                    source_row.trim_to_genomic_position_downstream(int(orf_coords[2]))
                else:
                    source_row.trim_to_genomic_position_upstream(int(orf_coords[2]))
                    source_row.trim_to_genomic_position_downstream(int(orf_coords[1]))

            except Exception as e:
                dropped += 1
                continue

            source_row.thick_start = source_row.start
            source_row.thick_stop = source_row.stop
            source_row.id_str = f"{entry.genomic_coords}⍼{entry.canonical_id}"
            f.write(f"{source_row}\n")
            x += 1

    return dropped, x


def write_results(df: DataFrame, threshold, output_path, with_ground_label=False) -> []:
    """

    :param df: results dataframe
    :param threshold: classification threshold
    :param output_path: output folder
    :param with_ground_label: If a ground label is available, let
    :return: Array with ranks of predictions relative to transcript ID
    """
    df = df.reset_index()
    current_id = ""
    current_count = 0
    ranks = []
    with open(f"{output_path}.uncollapsed.bed", 'w') as bed_file:
        if with_ground_label:
            for _, row in df.iterrows():
                if row['canonical_id'] != current_id:
                    current_id = row['canonical_id']
                    current_count = 0
                coords = "\t".join(row['genomic_coords'].split("|"))
                coords_split = row['genomic_coords'].split("|")
                rgb = ""

                if row['ground_label'] and row['class1_probability'] > threshold:
                    rgb = "0,255,0"

                elif not row['ground_label'] and row['class1_probability'] > threshold:
                    rgb = "255,0,0"

                elif row['ground_label'] and row['class1_probability'] < threshold:
                    rgb = "255,165,0"

                elif not row['ground_label'] and row['class1_probability'] < threshold:
                    rgb = "0,0,255"

                bed_file.write(
                    f"{coords[:-1]}\t{row['canonical_id']}_{current_count}_{row['class1_probability']:.4f}\t{row['class1_probability'] * 1000:.2f}\t{coords[-1]}\t{coords_split[1]}\t{coords_split[2]}\t{rgb}\n")
                ranks.append(current_count)
                current_count += 1

        else:
            for _, row in df.iterrows():
                if row['canonical_id'] != current_id:
                    current_id = row['canonical_id']
                    current_count = 0
                coords = "\t".join(row['genomic_coords'].split("|"))
                coords_split = row['genomic_coords'].split("|")

                if row['class1_probability'] > threshold:
                    rgb = "0,0,255"
                else:
                    rgb = "255,0,0"

                bed_file.write(
                    f"{coords[:-1]}\t{row['canonical_id']}_{current_count}_{row['class1_probability']:.4f}\t{row['class1_probability'] * 1000:.2f}\t{coords[-1]}\t{coords_split[1]}\t{coords_split[2]}\t{rgb}\n")
                ranks.append(current_count)
                current_count += 1

    return ranks


def map_bed(df, bed_path, output_path):
    """
    Takes a pandas dataframe and an input BED file and maps to two output files for coding and non-coding transcripts.
    Coding transcripts will have their thickStart / thickEnd fields modified to the ORF candidates.
    :param df: Dataframe after predictions and filtering
    :param bed_path: BED to map (input BED of pipeline)
    :param output_path: Path to write new BED files to
    :return:
    """
    # Assume this has already been filtered
    orfs = dict()
    for _, row in df.iterrows():
        if not orfs.get(row['canonical_id']):
            orfs[row['canonical_id']] = [(row['genomic_coords'], row['class1_probability'], row['toga_overrule'])]
        else:
            orf_list = orfs[row['canonical_id']]
            orf_list.append((row['genomic_coords'], row['class1_probability'], row['toga_overrule']))
            orfs[row['canonical_id']] = orf_list

    logger.info(f"Extraced ORFs for {len(orfs)} transcripts")

    with open(bed_path, 'r') as f:
        with open(f"{output_path}/coding.bed", 'w') as coding:
            with open(f"{output_path}/noncoding.bed", 'w') as noncoding:
                for line in f:
                    vals = line.strip().split("\t")
                    if orfs.get(vals[3]):
                        for orf in orfs[vals[3]]:
                            coords = orf[0].split("|")
                            # Set thickStart / stop
                            vals[6] = coords[1]
                            vals[7] = coords[2]

                            if orf[2]:
                                vals[3] = f"{vals[3]}_OVERRULED"
                            vals[4] = f"{orf[1] * 1000:.2f}"
                            new_row = "\t".join(vals)
                            coding.write(f"{new_row}\n")
                    else:
                        noncoding.write(line)
