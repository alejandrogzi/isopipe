"""
Utilities to clean TranslationAI output data and trim BED tracks to the predicted ORF borders
"""
import copy
import os
import subprocess
import logging
from bedutils.bed import BedRow

logger = logging.getLogger(__name__)
NF_FILE = "/beegfs/projects/hillerlab/genome/src/ORFTree/nextflow/chunked_translationai.nf"
NF_CONFIG = "/beegfs/projects/hillerlab/genome/src/ORFTree/nextflow/translation_ai_nxf.config"

def format_for_translationai(fasta_path: str, bed_path: str, output_path: str):
    """
    Takes a FASTA file whose headers only contain an ID and a bed file which contains features with corresponding IDs
    to reformat the FASTA file for translationAI. This function may be rendered obsolete in the future if matching
    changes to translationAI are made. In this case, please mark this function as obsolete.
    :param fasta_path:
    :param bed_path:
    :return:
    """
    id_to_coords = dict()
    reformatted_ids = []
    with open(bed_path, 'r') as bed_file:
        for line in bed_file:
            cols = line.split("\t")
            # {id} --> {chrom}:{start}-{stop}({strand})({id})(0, 0, )
            id_to_coords[cols[3]] = f">{cols[0]}:{cols[1]}-{cols[2]}({cols[5]})({cols[3]})(0,0, )"
            reformatted_ids.append(f">{cols[0]}:{cols[1]}-{cols[2]}({cols[5]})({cols[3]})(0,0, )")

    with open(fasta_path, 'r') as in_file:
        with open(output_path, 'w') as out_file:
            for line in in_file:
                if not line.startswith(">"):
                    out_file.write(line)
                else:
                    seq_id = line[1:].strip().split(" ")[0].replace("bosTau9_ncbiRefSeq_", "")
                    if id_to_coords.get(seq_id):
                        out_file.write(id_to_coords[seq_id] + "\n")
                    else:
                        logger.critical(f"The sequence ID {seq_id} from {fasta_path} has no corresponding entry in {bed_path}. "
                              f"This means that the file cannot be reformatted correctly. {output_path} is incomplete!")
                        exit(-1)

def clean_translationai_output(input_path, output_path) -> None:
    """
    Reformats the TranslationAI out at input_path into a TSV file with one prediction per row. Writes to output_path.
    Columns are:
    1. ID -> derived from the input ID with the suffix _{n} for the n-th predictions
    2. Start bp position
    3. End bp position
    4. Start probability
    5. End probability
    :param input_path: Raw TranslationAI output
    :param output_path: Path to reformatted TSV file
    """
    with open(input_path, 'r') as raw_file:
        with open(output_path, 'w') as out_file:
            for line in raw_file:
                in_data = line.split("\t")[0]
                my_id = in_data.split("(")[2].split(")")[0]
                variants = line.split("\t")[1:]
                for i, variant in enumerate(variants):
                    variant_cols = variant.strip().split(",")
                    out_file.write(
                        f"{my_id}_{i}\t{variant_cols[0]}\t{variant_cols[1]}\t{variant_cols[2]}\t{variant_cols[3]}\n")


def trim_translationai_bed(input_bed, raw_translation_ai, output_path):
    """
    Trim a BED file sharing IDs with the output of TranslationAI to the ORF predictions
    :param input_bed: Path to input bed file
    :param raw_translation_ai: Path to input TranslationAI file
    :param output_path: Path to output BED file
    :return:
    """
    ref_rows = dict()
    with open(input_bed, 'r') as reference:
        for line in reference:
            row = BedRow(line)
            ref_rows[row.id_str] = row

    variant_dict = dict()
    with open(raw_translation_ai, "r") as results:
        for line in results:
            in_data = line.split("\t")[0]
            my_id = in_data.split("(")[2].split(")")[0]
            variants = line.split("\t")[1:]
            variant_dict[my_id] = [
                {"orf_start": int(variants[0].split(',')[0]), "orf_stop": int(variants[0].split(',')[1]) + 3,
                 "start_score": float(variants[0].split(',')[2]), "stop_score": float(variants[0].split(',')[3])}]
            if len(variants) > 1:
                for variant in variants[1:]:
                    old_list = variant_dict[my_id]
                    old_list.append({"orf_start": int(variant.split(',')[0]), "orf_stop": int(variant.split(',')[1]),
                                     "start_score": float(variant.split(',')[2]),
                                     "stop_score": float(variant.split(',')[3])})
                    variant_dict[my_id] = old_list

    final_rows = []
    dropped = 0
    for key in variant_dict.keys():
        base_row = ref_rows[key]
        variants = variant_dict[key]
        for i, variant in enumerate(variants):
            try:
                row = copy.deepcopy(base_row)
                base_len = row.get_total_block_length()
                row.trim_bp_upstream(variant['orf_start'])
                row.trim_bp_downstream(base_len - variant['orf_stop'])
                avg_score = (variant["start_score"] + variant["stop_score"]) / 2
                row.score = int(avg_score * 1000)
                row.id_str = f"{row.id_str}_{i}"
                final_rows.append(str(row) + '\n')
            except Exception as e:
                dropped += 1

    if dropped > 0:
        logger.warning(f"Dropped {dropped} variants out of {len(variant_dict.keys())} variants")
    with open(output_path, 'w') as out_file:
        out_file.writelines(final_rows)


def translationai_main(transcript_fasta: str, transcript_bed: str, tmp_dir, use_nextflow=False):
    """
    TranslationAI wrapper
    :param transcript_fasta: Transcript FASTA file
    :param transcript_bed: Transcript alignments
    :param tmp_dir: Common temp dir to use
    :return:
    """
    os.mkdir(f"{tmp_dir}/TranslationAI")
    subprocess.run(f"cat {transcript_fasta} | linearize_fasta > {tmp_dir}/TranslationAI/linear.fa", shell=True)
    if use_nextflow:
        format_for_translationai(f"{tmp_dir}/TranslationAI/linear.fa", f"{transcript_bed}", f"{tmp_dir}/TranslationAI/linear.fa.formatted")
        translation_ai_run = subprocess.run(f"cd {tmp_dir}/TranslationAI/ && nextflow {NF_FILE} -c {NF_CONFIG} --query linear.fa.formatted --out ./out.translationai -work-dir TEMP_TRANSLATIONAI", shell=True)
    else:
        translation_ai_run = subprocess.run(f"cd {tmp_dir}/TranslationAI/ && run_translation_ai ../../{transcript_bed} ../../{transcript_fasta}.linear 1000 out.translationai", shell=True)

    if translation_ai_run.returncode != 0:
        logger.critical("TranslationAI failed, dumping output:")
        print(translation_ai_run.stdout ,"=====", translation_ai_run.stdout)

    logger.info("Producing rich tracks for translationai")
    trim_translationai_bed(transcript_bed, f"{tmp_dir}/TranslationAI/out.translationai",f"{tmp_dir}/TranslationAI/translationai.trimmed.bed")

    logger.info("Producing clean output file")
    clean_translationai_output(f"{tmp_dir}/TranslationAI/out.translationai",f"{tmp_dir}/TranslationAI/out.translationai.tsv")
