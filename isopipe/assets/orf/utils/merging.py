"""
Utilities to load data from the three tools into the common dataframe.
"""

import pandas as pd
import logging

logger = logging.getLogger(__name__)


def filter_by_relative_score_strict(group):
    """
    Filter scores by difference to preceeding scores.
    :param group:
    :return:
    """
    sorted_group = group.sort_values(
        by="class1_probability", ascending=True
    ).reset_index(drop=True)
    mask = sorted_group["class1_probability"] >= (
        sorted_group["class1_probability"].iloc[-1] - 0.025
    )
    return sorted_group[mask]


def reformat_id_orfipy(id_str, include_orf=True) -> str:
    """
    Reformat IDs produced by orfipy to a TransDecoder compatible format.
    :param id_str: Original ID
    :param include_orf: Trim the ORF tail appended by the nested ORF finder.
    :return: Reformatted ID
    """
    # base_id = id_str.split(":")[0].replace("_type", "")
    base_id = f"{id_str.split('_[')[0]}"
    base_id = base_id.replace("_ORF.", ".p")

    if id_str.split("_")[-1].startswith("ORF") and include_orf:
        base_id = f"{base_id}_{id_str.split('_')[-1]}"
    return base_id


def load_blast(input_fasta_path, results_path, exported_path) -> pd.DataFrame:
    """
    Load BLAST data for the common dataframe
    :param input_fasta_path: BLAST input .pep file
    :param results_path: BLASTP output
    :param exported_path: Path to exported CDS BED file
    :return:
    """
    blast_dict = dict()
    blast_data_dict = dict()
    dropped_rows = 0

    def reformat_id(id_str):
        # base_id = id_str.split(":")[0].replace("_type", "")
        base_id = f"{id_str.split('_')[0]}_{id_str.split('_')[1]}"

        if id_str.split("_")[-1].startswith("ORF"):
            base_id = f"{base_id}_{id_str.split('_')[-1]}"
        return base_id

    # Load BLAST fmt6 output file as plain text
    with open(results_path, "r") as blast_output_file:
        for line in blast_output_file:
            cols = line.strip().split("\t")
            my_id = reformat_id_orfipy(cols[0])
            if not blast_data_dict.get(my_id):
                pid = float(cols[2])
                e_value = float(cols[-2])
                # offset = int(cols[8]) - int(cols[6])
                offset = int(cols[8]) - int(cols[6])
                alignment_len = int(cols[7]) - int(cols[6]) + 1
                blast_data_dict[my_id] = {
                    "blast_id": my_id,
                    "blast_pid": pid,
                    "blast_e-value": e_value,
                    "blast_offset": offset,
                    "blast_alignment_len": alignment_len,
                }
            else:
                continue

    query_lens = dict()
    with open(input_fasta_path, "r") as pep_file:
        curr_id = ""
        curr_seq_len = -1
        for line in pep_file:
            if line.startswith(">"):
                # vals = line.split(" ")
                # curr_id = reformat_id(vals[0][1:].strip())
                curr_id = reformat_id_orfipy(line[1:].strip())
            else:
                query_lens[curr_id] = len(line.strip()[:-1])

    for key in blast_data_dict.keys():
        prev_data = blast_data_dict[key]
        # prev_data['blast_percentage_aligned'] = blast_data_dict[key]['blast_alignment_len'] / query_lens[blast_data_dict[key]['blast_id']]
        prev_data["blast_percentage_aligned"] = (
            blast_data_dict[key]["blast_alignment_len"] / query_lens[key]
        )
        blast_data_dict[key] = prev_data

    no_hits = []
    blast_data_for_df = []

    # This opens the exported BED tracks -> should have fixed IDs
    logging.info("Loading trimmed alignments")
    with open(exported_path, "r") as blast_file:
        # ID, start, length
        # The top level ORF must be non-nested
        current_top_level_orf = ("", 0, 0)
        for line in blast_file:
            vals = line.strip().split("\t")
            chrom = vals[0]
            my_id = vals[3]
            if blast_data_dict.get(my_id):
                # Thick start / end of BED track
                coding_start = int(vals[6])
                coding_stop = int(vals[7])
                # Strand from bed track
                strand = vals[5]
                my_data = blast_data_dict[my_id]

                # Does the ID contain "ORF" -> not a top level ORF
                if my_data["blast_id"].find("ORF") != -1:
                    my_data["blast_nested"] = int(my_data["blast_id"].split("ORF")[-1])
                    my_data["blast_id"] = "_".join(my_data["blast_id"].split("_")[:-1])
                else:
                    my_data["blast_nested"] = 0
                    if strand == "+":
                        current_top_level_orf = (
                            my_data["blast_id"],
                            coding_start,
                            coding_stop - coding_start,
                        )
                    elif strand == "-":
                        current_top_level_orf = (
                            my_data["blast_id"],
                            coding_stop,
                            coding_stop - coding_start,
                        )
                    else:
                        print(f"Data corrupted: {strand} is not a valid strand")
                    # print(current_top_level_orf)
                if my_data["blast_id"].find(".p") != -1:
                    my_data["blast_id"] = ".".join(my_data["blast_id"].split(".")[:-1])

                my_data["blast_coords"] = (
                    f"{chrom}|{coding_start}|{coding_stop}|{strand}"
                )
                my_data["canonical_id"] = my_data["blast_id"]
                my_data["genomic_coords"] = (
                    f"{chrom}|{coding_start}|{coding_stop}|{strand}"
                )

                try:
                    # Find an offset as the difference between the coding start of the current transcript and the top-
                    # level and divide by the maximum length.
                    if strand == "+":
                        my_data["blast_nested_offset"] = (
                            abs(coding_start - current_top_level_orf[1])
                            / current_top_level_orf[2]
                        )
                        # print(my_id, abs(coding_start - current_top_level_orf[1]), current_top_level_orf[2])
                    else:
                        my_data["blast_nested_offset"] = (
                            abs(coding_stop - current_top_level_orf[1])
                            / current_top_level_orf[2]
                        )
                        # print(my_id, abs(coding_stop - current_top_level_orf[1]), current_top_level_orf[2])
                except ZeroDivisionError:
                    dropped_rows += 1

                blast_data_for_df.append(my_data)
        else:
            no_hits.append(my_id)

    df = pd.DataFrame(
        data=blast_data_for_df,
        columns=[
            "genomic_coords",
            "canonical_id",
            "blast_id",
            "blast_nested",
            "blast_nested_offset",
            "blast_pid",
            "blast_e-value",
            "blast_offset",
            "blast_percentage_aligned",
        ],
    )
    df.set_index(["genomic_coords", "canonical_id"])
    logger.info(f"Finished loading BLAST data")
    logger.warning(f"Dropped {len(no_hits)} | {dropped_rows} rows (safely)")
    return df


def load_translation_ai(results_file, trimmed_bed) -> pd.DataFrame:
    """
    Load TranslationAI data for the common dataframe.
    :param results_file: Cleaned TranslationAI results file
    :param trimmed_bed: Trimmed BED file
    :return:
    """
    translation_ai_data_dict = dict()
    translation_ai_for_df = []
    with open(results_file, "r") as translation_ai_data_file:
        for line in translation_ai_data_file:
            cols = line.strip().split("\t")
            my_id = cols[0]
            orf_start_score = float(cols[3])
            orf_stop_score = float(cols[4].strip())
            translation_ai_data_dict[my_id] = {
                "translationAI_orf_start": orf_start_score,
                "translationAI_orf_stop": orf_stop_score,
            }

    with open(trimmed_bed, "r") as translation_file:
        for line in translation_file:
            vals = line.strip().split("\t")
            chrom = vals[0]
            transcript_id = vals[3]
            coding_start = int(vals[1])
            coding_stop = int(vals[2])
            strand = vals[5]
            my_data = translation_ai_data_dict[transcript_id]
            my_data["translationAI_id"] = f"{'_'.join(transcript_id.split('_')[:-1])}"
            my_data["translationAI_coords"] = (
                f"{chrom}|{coding_start}|{coding_stop}|{strand}"
            )
            translation_ai_for_df.append(my_data)

    df = pd.DataFrame(
        data=translation_ai_for_df,
        columns=[
            "translationAI_coords",
            "translationAI_id",
            "translationAI_orf_start",
            "translationAI_orf_stop",
        ],
    )
    return df


def load_toga(toga_results_file, query_annotation, overlap_file=None) -> pd.DataFrame:
    """
    Load TOGA data for the common dataframe.
    :param toga_results_file: Merged TOGA results file
    :param query_annotation: Query annotation BED file
    :param overlap_file: Overlap dual-bed file produced by bedtools intersect
    :return:
    """
    toga_data = pd.read_csv(toga_results_file, delimiter="\t", header=0, index_col=0)
    toga_data = toga_data.rename(
        columns={
            "pid": "toga_pid",
            "blosum": "toga_blosum",
            "label": "toga_label",
            "masked": "toga_masked",
        }
    )
    toga_id_dict = dict()
    toga_for_df = []
    with open(query_annotation, "r") as toga_file:
        for line in toga_file:
            vals = line.strip().split("\t")
            chrom = vals[0]
            transcript_id = vals[3]
            coding_start = int(vals[6])
            coding_stop = int(vals[7])
            strand = vals[5]
            my_data = dict(toga_data.loc[transcript_id])
            my_data["toga_id"] = transcript_id
            my_data["toga_coords"] = f"{chrom}|{coding_start}|{coding_stop}|{strand}"
            toga_for_df.append(my_data)
            toga_id_dict[transcript_id] = line

    df = pd.DataFrame(
        data=toga_for_df,
        columns=[
            "toga_coords",
            "toga_id",
            "toga_pid",
            "toga_blosum",
            "toga_label",
            "toga_masked",
        ],
    )

    # Load TOGA overlap data if provided
    if overlap_file:
        overlap_list = []
        with open(overlap_file, "r") as overlap_file:
            for line in overlap_file:
                vals = line.strip().split("\t")
                a_id = vals[3]
                b_id = vals[15]
                bp_overlap = int(vals[24])
                overlap_list.append((a_id, b_id, bp_overlap))

        overlap_df = pd.DataFrame(
            overlap_list, columns=["toga_id", "canonical_id", "toga_overlap_bp"]
        )
        df = df.merge(overlap_df, left_on=["toga_id"], right_on=["toga_id"], how="left")

    return df


def load_toga_new(
    toga_results_file, query_annotation, overlap_file=None
) -> pd.DataFrame:
    """
    Load TOGA data for the common dataframe.
    :param toga_results_file: Merged TOGA results file
    :param query_annotation: Query annotation BED file
    :param overlap_file: Overlap dual-bed file produced by bedtools intersect
    :return:
    """
    TOGA_LABELS = ["TOGA_ID", "label", "pid", "blosum"]
    toga_data = pd.read_csv(
        toga_results_file, delimiter="\t", names=TOGA_LABELS, index_col=0
    )

    toga_data = toga_data.rename(
        columns={
            "pid": "toga_pid",
            "blosum": "toga_blosum",
            "label": "toga_label",
            "masked": "toga_masked",
        }
    )
    toga_id_dict = dict()
    toga_for_df = []
    with open(query_annotation, "r") as toga_file:
        for line in toga_file:
            vals = line.strip().split("\t")
            chrom = vals[0]
            transcript_id = vals[3]
            coding_start = int(vals[6])
            coding_stop = int(vals[7])
            strand = vals[5]
            my_data = dict(toga_data.loc[transcript_id])
            my_data["toga_id"] = transcript_id
            my_data["toga_coords"] = f"{chrom}|{coding_start}|{coding_stop}|{strand}"
            toga_for_df.append(my_data)
            toga_id_dict[transcript_id] = line

    df = pd.DataFrame(
        data=toga_for_df,
        columns=[
            "toga_coords",
            "toga_id",
            "toga_pid",
            "toga_blosum",
            "toga_label",
            "toga_masked",
        ],
    )

    # Load TOGA overlap data if provided
    if overlap_file:
        overlap_list = []
        with open(overlap_file, "r", encoding="utf-8") as overlap_file:
            for line in overlap_file:
                vals = line.strip().split("\t")
                a_id = vals[3]
                b_id = vals[15]
                bp_overlap = int(vals[24])

                # If you use ⍼ in your IDs outside of this script, get help
                id_components = b_id.split("⍼")

                overlap_list.append(
                    (a_id, id_components[0], id_components[1], bp_overlap)
                )

        overlap_df = pd.DataFrame(
            overlap_list,
            columns=["toga_id", "genomic_coords", "canonical_id", "toga_overlap_bp"],
        )
        df = df.merge(overlap_df, left_on=["toga_id"], right_on=["toga_id"], how="left")

    return df


def load_toga_two(
    toga_results_file, query_annotation, toga_masked_file, overlap_file=None
) -> pd.DataFrame:
    """
    Load TOGA data for the common dataframe.
    :param toga_results_file: Merged TOGA results file
    :param query_annotation: Query annotation BED file
    :param overlap_file: Overlap dual-bed file produced by bedtools intersect
    :return:
    """
    TOGA_LABELS = ["TOGA_ID", "label", "pid", "blosum"]
    toga_data = pd.read_csv(
        toga_results_file,
        delimiter="\t",
        names=TOGA_LABELS,
        usecols=[0, 1, 2, 3],
        index_col=0,
        compression='gzip'
    )

    toga_data = add_masked_labels(toga_data, toga_masked_file)
    toga_data = toga_data.rename(
        columns={
            "pid": "toga_pid",
            "blosum": "toga_blosum",
            "label": "toga_label",
            "masked": "toga_masked",
        }
    )
    toga_id_dict = dict()
    toga_for_df = []
    logger.info("Loading TOGA coordinates from query annotation")

    with open(query_annotation, "r") as toga_file:
        for line in toga_file:
            vals = line.strip().split("\t")
            chrom = vals[0]
            transcript_id = vals[3]
            coding_start = int(vals[6])
            coding_stop = int(vals[7])
            strand = vals[5]
            my_data = dict(toga_data.loc[transcript_id])
            my_data["toga_id"] = transcript_id
            my_data["toga_coords"] = f"{chrom}|{coding_start}|{coding_stop}|{strand}"
            toga_for_df.append(my_data)
            toga_id_dict[transcript_id] = line

    df = pd.DataFrame(
        data=toga_for_df,
        columns=[
            "toga_coords",
            "toga_id",
            "toga_pid",
            "toga_blosum",
            "toga_label",
            "toga_masked",
        ],
    )

    # Load TOGA overlap data if provided
    if overlap_file:
        logger.info("Loading overlap information")
        overlap_list = []
        with open(overlap_file, "r", encoding="utf-8") as overlap_file:
            for line in overlap_file:
                vals = line.strip().split("\t")
                a_id = vals[3]
                b_id = vals[15]
                bp_overlap = int(vals[24])

                # If you use ⍼ in your IDs outside of this script, get help
                id_components = b_id.split("⍼")

                overlap_list.append(
                    (a_id, id_components[0], id_components[1], bp_overlap)
                )

        overlap_df = pd.DataFrame(
            overlap_list,
            columns=["toga_id", "genomic_coords", "canonical_id", "toga_overlap_bp"],
        )
        df = df.merge(overlap_df, left_on=["toga_id"], right_on=["toga_id"], how="left")

    return df


def add_masked_labels(df, masked_path):
    """
    Load TOGA masked properties from a separate file
    :param df:
    :param masked_path:
    :return:
    """
    mask_data = pd.read_csv(masked_path, delimiter="\t", usecols=[0, 6])
    mask_data.rename(
        columns={"projection": "TOGA_ID", "query_codon": "masked"}, inplace=True
    )

    mask_data["masked"] = mask_data["masked"] == "MASKED"

    mask_data.drop_duplicates(inplace=True)
    mask_data.drop_duplicates(subset="TOGA_ID", inplace=True)
    mask_data.set_index("TOGA_ID", inplace=True)
    df = df.merge(mask_data, how="left", left_index=True, right_index=True)
    df["masked"] = df["masked"].fillna(False).astype(bool)
    return df
