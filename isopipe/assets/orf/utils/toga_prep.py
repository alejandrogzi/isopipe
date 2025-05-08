"""
Prepare TOGA data
"""
import pandas as pd

from bedutils.bed import BedRow
from copy import deepcopy


def export_df_orfs(df, input_bed, output_path):
    # populate BED dict
    transcript_dict = dict()
    with open(input_bed, 'r') as f:
        for line in f:
            curr_row = BedRow(line)
            transcript_dict[curr_row.id_str] = curr_row

    with open(output_path,'w') as f:
        for entry in df.itertuples():
            orf_coords = entry['genomic_coords'].split("|")
            source_row = deepcopy(transcript_dict.get(entry['canonical_id']))

            if orf_coords[-1] == '+':
                source_row.trim_to_genomic_position_upstream(source_row.genomic_position_to_sequence_position(int(orf_coords[1])))
                source_row.trim_to_genomic_position_downstream(source_row.genomic_position_to_sequence_position(int(orf_coords[2])))
            elif orf_coords[-1] == '-':
                source_row.trim_to_genomic_position_upstream(source_row.genomic_position_to_sequence_position(int(orf_coords[2])))
                source_row.trim_to_genomic_position_downstream(source_row.genomic_position_to_sequence_position(int(orf_coords[1])))
            else:
                raise ValueError("Invalid strand symbol")

            # Trim BED row to constraints of ORF
            """if orf_coords[-1] == '+':
                source_row.trim_upstream(source_row.genomic_pos_to_seq_pos(int(orf_coords[1])))
                source_row.trim_downstream(source_row.genomic_pos_to_seq_pos(int(orf_coords[2])))
            else:
                source_row.trim_upstream(source_row.genomic_pos_to_seq_pos(int(orf_coords[2])))
                source_row.trim_downstream(source_row.genomic_pos_to_seq_pos(int(orf_coords[1])))"""

            f.write(str(source_row))



def combine_toga_data(pid_blosum_data, label_data, masking_data, output_path):
    """

    :param pid_blosum_data:
    :param label_data:
    :param masking_data:
    :param output_path:
    :return:
    """
    # Read PID / BLOSUM
    df = pd.read_csv(pid_blosum_data, delimiter="\t", header=0, index_col=0)
    # Read TOGA labels
    label_df = pd.read_csv(label_data, delimiter="\t", index_col=0, names=["projection", "label"])
    df = df.merge(label_df, how="outer", left_index=True, right_index=True)
    df = df.dropna(subset=['pid', 'blosum'])


    masked_dict = dict()
    with open(masking_data, 'r') as masked:
        for line in masked:
            cols = line.split("\t")
            my_id = f"{cols[0][2:]}.{cols[1]}"
            masked_dict[my_id] = {"masked": True}

    # Load TOGA masking
    masked_df = pd.DataFrame.from_dict(masked_dict, orient="index")
    df = df.merge(masked_df, left_index=True, right_index=True, how='left')
    df['masked'] = df['masked'].fillna(False).astype(bool)

    df.to_csv(output_path, sep='\t')
