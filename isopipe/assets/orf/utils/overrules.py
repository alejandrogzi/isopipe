import copy
import logging

stop_codons = ["TAA", "TAG", "TGA"]
logger = logging.getLogger(__name__)


def check_compatible(
    row_a,
    row_b,
    prediction,
    target_coords,
    sequence,
    allowed_stop=None,
    allow_trip=False,
):
    if row_a.strand == "+":
        if not row_a.is_genomic_position_in_block(row_b.start + 1, ignore_thick=True):
            return None, False
        seq_start = row_a.genomic_position_to_sequence_position(row_b.start)

        codon_match = -1
        has_tripped = False
        for i in range(seq_start, len(sequence) - 2, 3):
            codon = sequence[i : i + 3]
            if codon in stop_codons:
                codon_match = row_a.sequence_position_to_genomic_position(i)
                if allowed_stop:
                    if i == allowed_stop:
                        continue
                if not has_tripped:
                    has_tripped = True
                    continue
                break

        my_row = copy.deepcopy(row_a)
        try:
            my_row.trim_bp_upstream(seq_start)
            my_row.trim_bp_downstream(
                abs(my_row.genomic_position_to_sequence_position(codon_match) + 3)
            )
        except ValueError:
            return None, False

        if codon_match + 3 != int(target_coords.split("|")[2]):
            my_row.rgb = "255,0,0"
            compatible = False
        else:
            my_row.rgb = "0,255,0"
            compatible = True

        my_row.id_str = f"{prediction['canonical_id']}_{round(prediction['class1_probability'], 4)}_{row_b.id_str}"

        return my_row, compatible

    else:
        if row_a.is_genomic_position_in_block(row_b.stop, ignore_thick=True):
            seq_start = row_a.genomic_position_to_seq_position(row_b.stop)
            codon_match = -1

            has_tripped = False
            for i in range(seq_start, len(sequence) - 2, 3):
                codon = sequence[i : i + 3]
                if codon in stop_codons:
                    codon_match = row_a.seq_position_to_genomic_position(i)
                    if not has_tripped:
                        has_tripped = True
                        continue
                    if allowed_stop:
                        if i == allowed_stop:
                            continue
                    break

            my_row = copy.deepcopy(row_a)
            try:
                my_row.trim_bp_upstream(seq_start)
                base_len = my_row.get_total_block_length()
                # my_row.trim_downstream(abs(my_row.genomic_position_to_sequence_position(codon_match) + 3))
                my_row.trim_bp_downstream(
                    base_len
                    - abs(my_row.genomic_position_to_sequence_position(codon_match) + 3)
                )
            except ValueError:
                return None, False

            if codon_match - 3 != int(target_coords.split("|")[1]):
                my_row.rgb = "255,0,0"
                compatible = False
            else:
                my_row.rgb = "0,255,0"
                compatible = True

            my_row.id_str = f"{prediction['canonical_id']}_{round(prediction['class1_probability'], 4)}_{row_b.id_str}"
            return my_row, compatible

        else:
            return None, False


def toga_overrule_a(
    unique_canonic, grouped_df, transcript_rows, toga_rows, sequences, threshold=0.03
):
    overrule_rows = []

    for x in unique_canonic:
        predictions = (
            grouped_df.get_group(x)
            .sort_values(by="class1_probability", ascending=False)
            .copy()
        )

        # Skip transcripts with one annotated ORF
        if predictions.iloc[0]["class1_probability"] > threshold:
            continue
        filtered_labels = {
            "FI": predictions[predictions["toga_label"] == "FI"]
            .sort_values(by="toga_pid", ascending=False)
            .drop_duplicates()
            .copy(),
            "I": predictions[predictions["toga_label"] == "I"]
            .sort_values(by="toga_pid", ascending=False)
            .drop_duplicates()
            .copy(),
            "PI": predictions[predictions["toga_label"] == "PI"]
            .sort_values(by="toga_pid", ascending=False)
            .drop_duplicates()
            .copy(),
        }

        overruled = False
        for label in ["FI", "I", "PI"]:
            if overruled:
                break

            if len(filtered_labels[label]) > 0:
                for y in range(len(filtered_labels[label])):
                    candidate = filtered_labels[label].iloc[y].copy()
                    ref_row = transcript_rows[
                        candidate["canonical_id"].replace("mm10_ncbiRefSeq_", "")
                    ]
                    toga_row = toga_rows[candidate["toga_id"]]

                    try:
                        my_row, compatible = check_compatible(
                            ref_row,
                            toga_row,
                            candidate,
                            candidate["toga_coords"],
                            sequences[x],
                            allow_trip=True,
                        )
                        if compatible:
                            overrule_rows.append(candidate)
                            break
                    except Exception as e:
                        pass
    return overrule_rows


def toga_overrule_b(
    unique_canonic, grouped_df, transcript_rows, toga_rows, sequences, threshold=0.03
):
    overrule_rows = []
    overrule_ids = []
    dropped = 0
    for x in unique_canonic:
        try:
            predictions = (
                grouped_df.get_group(x)
                .sort_values(by="class1_probability", ascending=False)
                .copy()
            )
            best_prediction = predictions.iloc[0].copy()

            filtered_labels = {
                "FI": predictions[predictions["toga_label"] == "FI"]
                .sort_values(by="toga_pid", ascending=False)
                .copy(),
                "I": predictions[predictions["toga_label"] == "I"]
                .sort_values(by="toga_pid", ascending=False)
                .copy(),
                "PI": predictions[predictions["toga_label"] == "PI"]
                .sort_values(by="toga_pid", ascending=False)
                .copy(),
            }
            for label in ["FI", "I", "PI"]:
                label_frame = filtered_labels[label]
                if (
                    len(label_frame) > 0
                    and best_prediction["toga_id"] != best_prediction["canonical_id"]
                ):
                    ref_row = transcript_rows[
                        best_prediction["canonical_id"].replace("mm10_ncbiRefSeq_", "")
                    ]
                    toga_row = toga_rows[label_frame.iloc[0]["toga_id"]]

                    # are all TOGAs in a class masked -> correct
                    if (label_frame["toga_masked"] == True).all() and len(
                        label_frame
                    ) > 0:
                        best_coords = best_prediction["genomic_coords"].split("|")

                        # What about different strands?
                        if best_coords[3] == "+":
                            best_stop = int(best_coords[2])
                            toga_stop = toga_row.stop
                        else:
                            best_stop = int(best_coords[1])
                            toga_stop = toga_row.start

                        if ref_row.strand == "+":
                            if toga_stop != best_stop:
                                over_row, is_compatible = check_compatible(
                                    ref_row,
                                    toga_row,
                                    best_prediction,
                                    best_prediction["toga_coords"],
                                    sequences[x],
                                    allowed_stop=ref_row.genomic_position_to_sequence_position(
                                        int(
                                            best_prediction["genomic_coords"].split(
                                                "|"
                                            )[2]
                                        )
                                        - 3
                                    ),
                                )
                            else:
                                over_row, is_compatible = check_compatible(
                                    ref_row,
                                    toga_row,
                                    best_prediction,
                                    best_prediction["toga_coords"],
                                    sequences[x],
                                )
                        else:
                            if toga_stop != best_stop:
                                over_row, is_compatible = check_compatible(
                                    ref_row,
                                    toga_row,
                                    best_prediction,
                                    best_prediction["toga_coords"],
                                    sequences[x],
                                    allowed_stop=ref_row.genomic_position_to_sequence_position(
                                        int(
                                            best_prediction["genomic_coords"].split(
                                                "|"
                                            )[1]
                                        )
                                    )
                                    - 3,
                                )
                            else:
                                over_row, is_compatible = check_compatible(
                                    ref_row,
                                    toga_row,
                                    best_prediction,
                                    sequences[x],
                                    best_prediction["toga_coords"],
                                )

                        if is_compatible:
                            overrule_ids.append(x)
                            best_prediction["genomic_coords"] = best_prediction[
                                "toga_coords"
                            ]
                            overrule_rows.append(best_prediction)

                    elif (label_frame["toga_masked"] == True).any() and len(
                        label_frame
                    ) > 0:
                        best_coords = best_prediction["genomic_coords"].split("|")
                        # What about different strands?
                        if best_coords[3] == "+":
                            best_stop = int(best_coords[2])
                            toga_stop = toga_row.stop
                        else:
                            best_stop = int(best_coords[1])
                            toga_stop = toga_row.start

                        if ref_row.strand == "+":
                            if toga_stop != best_stop:
                                over_row, is_compatible = check_compatible(
                                    ref_row,
                                    toga_row,
                                    best_prediction,
                                    best_prediction["toga_coords"],
                                    sequences[x],
                                    allowed_stop=ref_row.genomic_position_to_sequence_position(
                                        int(
                                            best_prediction["genomic_coords"].split(
                                                "|"
                                            )[2]
                                        )
                                        - 3
                                    ),
                                )
                            else:
                                over_row, is_compatible = check_compatible(
                                    ref_row,
                                    toga_row,
                                    best_prediction,
                                    best_prediction["toga_coords"],
                                    sequences[x],
                                )
                            if is_compatible:
                                overrule_ids.append(best_prediction)
                                best_prediction["genomic_coords"] = best_prediction[
                                    "toga_coords"
                                ]
                                overrule_rows.append(best_prediction)
                        else:
                            if toga_stop != best_stop:
                                over_row, is_compatible = check_compatible(
                                    ref_row,
                                    toga_row,
                                    best_prediction,
                                    best_prediction["toga_coords"],
                                    sequences[x],
                                    allowed_stop=ref_row.genomic_position_to_sequence_position(
                                        int(
                                            best_prediction["genomic_coords"].split(
                                                "|"
                                            )[1]
                                        )
                                    )
                                    - 3,
                                )
                            else:
                                over_row, is_compatible = check_compatible(
                                    ref_row,
                                    toga_row,
                                    best_prediction,
                                    best_prediction["toga_coords"],
                                    sequences[x],
                                )

                            if is_compatible:
                                overrule_ids.append(best_prediction)
                                best_prediction["genomic_coords"] = best_prediction[
                                    "toga_coords"
                                ]
                                overrule_rows.append(best_prediction)
        except Exception as e:
            dropped += 1

    if dropped > 0:
        logger.warning(
            f"Dropped {dropped} overrule candidates for TOGA overrule Type B"
        )
    return overrule_rows
