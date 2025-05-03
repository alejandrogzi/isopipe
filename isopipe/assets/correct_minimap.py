from __future__ import print_function

#!/beegfs/projects/hillerlab/genome/src/correctMinimap/.venv/bin/python3

__author__ = "Lucas Koch"
__maintainer__ = "Lucas Koch"
__credits__ = ["Lucas Koch", "Michael Hiller", "Ariadna Morales"]
__license__ = "MIT"
__version__ = "1.0"


import shutil
import subprocess
from argparse import ArgumentParser
from itertools import groupby


parser = ArgumentParser(
    prog="correctMinimap",
    description="correctMinimap detects aligned transcripts (IsoSeqs) that end at right at an exon boundary and have a small 3' clipped sequence. The script then tests if the 3'clip sequence matches the beginning of a downstream exon. If so, it shifts the alignment end and adjusts the CIGAR.",
)
parser.add_argument(
    "toga",
    type=str,
    metavar="TOGA",
    help="Path to the TOGA annotation (e.g. $genomePath/gbdb-HL/$ref/TOGA/vs_$query/query_annotation.bed)",
)
parser.add_argument(
    "minimap",
    type=str,
    metavar="MINIMAP",
    help="File listing the transcripts mapped with minimap2 -eqx. Input format can be .bam or .sam (will be autodetected)",
)
parser.add_argument(
    "assembly",
    type=str,
    metavar="assembly",
    help="Name of the assembly to which the transcripts are mapped (e.g. hg38)",
)
parser.add_argument(
    "output",
    type=str,
    metavar="output",
    help="Name of output file containing all transcripts, some with corrected alignments (format equals input format. bam or sam)",
)
parser.add_argument(
    "-clip_cutoff",
    type=int,
    default=2,
    metavar="int",
    help="Minimum 3' clip length of a transcript to be considered (default is 2)",
)
parser.add_argument(
    "-wiggle",
    type=int,
    default=2,
    metavar="int",
    help="consider transcripts that end +-$wiggle bp at an exon end. Default is 2 (meaning transcripts ending -1/-2/+1/+2 bp of an exon end will be considered)",
)
parser.add_argument(
    "--debug",
    action="store_const",
    const=True,
    metavar="keep TEMP folder",
    help="If set, output a corrected.sam and discarded.sam that have only the "
    "IsoSeqs that were candidates (end at an exon end) and were or were not corrected",
)
parser.add_argument(
    "--keep_temp",
    action="store_const",
    const=True,
    metavar="keep TEMP folder",
    help="Do not delete the temporary directory",
)


def rev_comp(my_seq):
    """
    Return reverse complement of a given sequence
    :param my_seq:
    :return:
    """
    rev_seq = ""
    rev_char = {"A": "T", "C": "G", "G": "C", "T": "A"}
    for char in my_seq[::-1]:
        rev_seq = f"{rev_seq}{rev_char[char.upper()]}"
    return rev_seq


def filter_internal_exons(bed_path, temp_path):
    """
    Dissect a BED12 file producing a BED file containing only exons which do not start / end the transcript. Exon labels
    are appended a '#' and a consecutive number (counter, per transcript)
    :param bed_path: BED12 filepath
    :param temp_path: temp dir from main()
    :return:
    """
    transcript_dict = dict()
    total_expanded = 0

    # Subdivide BED12 features
    subprocess.check_output(
        f"bed12ToBed6 -i {bed_path} > {temp_path}/query_annotation.bed6", shell=True
    )

    with open(bed_path, "r") as file:
        for line in file:
            cols = line.split("\t")
            chrom = cols[0]
            start = int(cols[1])
            stop = int(cols[2])
            my_id = cols[3]

            transcript_dict[my_id] = (chrom, start, stop)
            total_expanded += 1

    total_discarded = 0
    # WARN: we assume sorted BED file
    # Append index #n to all exons, then direct internal exons to separate file
    with open(f"{temp_path}/query_annotation.bed6", "r") as six_file:
        current_id = ""
        current_count = 0
        with open(f"{temp_path}/indexed_exons.bed", "w") as all_file:
            with open(f"{temp_path}/internal_exons.bed", "w") as internal_file:
                for line in six_file:
                    cols = line.split("\t")
                    my_id = cols[3]

                    # reset counter
                    if my_id != current_id:
                        current_id = my_id
                        current_count = 0

                    # translate to indexed BED file
                    if (
                        int(cols[1]) != transcript_dict[my_id][1]
                        and int(cols[2]) != transcript_dict[my_id][2]
                    ):
                        cols[3] += f"#{current_count}"
                        internal_file.write("\t".join(cols))
                    else:
                        total_discarded += 1
                        cols[3] += f"#{current_count}"
                    current_count += 1
                    all_file.write("\t".join(cols))


def parse_named_fasta(fa_path):
    """
    Load a fasta file into a sequence dict
    :param fa_path: Path to fasta file
    :return: dict, keys are FASTA headers, values are sequences
    """
    exons_seqs_dict = dict()
    with open(fa_path, "r") as file:
        current_seq = ""  # FASTA Seq
        current_id = ""  # FASTA ID
        lines = file.readlines()
        for x, line in enumerate(lines):
            if line.startswith(">") or x == len(lines) - 1:
                if current_id != "":
                    exons_seqs_dict[current_id] = current_seq
                current_id = line[1:].strip()
                current_seq = ""
            elif line.strip() != "":
                current_seq += line.strip()
    return exons_seqs_dict


def parse_sam_into_dicts(original_sam: str, tmp_dir: str, cutoff=0, wiggle=0):
    """
    Parse SAM into Stop pos / clip dict
    :param cutoff:
    :param original_sam:
    :param tmp_dir:
    :return:
    """
    stop_pos_dict = dict()
    clip_dict = dict()
    with open(original_sam, "r") as file:
        for line in file:
            if not line.startswith("@"):
                poly_a = 0
                clip = -1
                cols = line.split("\t")
                my_id = cols[0]
                unwrapped = my_id.split("_")
                for data in unwrapped:
                    if data.startswith("3Clip"):
                        clip = int(data.replace("3Clip", ""))
                    if data.startswith("PolyA") and not data.startswith("PolyARead"):
                        poly_a = int(data.replace("PolyA", ""))
                if clip >= cutoff:
                    seq = cols[9]
                    # (+) strand
                    if cols[1] != "16":
                        clip_dict[my_id] = seq[
                            -(clip + poly_a) - abs(wiggle) : -(clip + poly_a) + clip
                        ]
                    # (-) strand
                    else:
                        clip_dict[my_id] = seq[poly_a : poly_a + clip + abs(wiggle)]

                if clip == -1:
                    print(
                        f"The information abourt clip length is not present in the SAM row '{line}'"
                    )

    # Load the end positions of the SAM entries
    # (from BED for convenience)
    with open(f"{tmp_dir}/minimap.sam.bed8", "r") as file:
        for line in file:
            cols = line.split("\t")
            strand = cols[5]
            chrom = cols[0]
            my_id = cols[3]

            unwrapped = my_id.split("_")
            for data in unwrapped:
                if data.startswith("3Clip"):
                    clip = int(data.replace("3Clip", ""))

            # 3' coord
            if strand == "+":
                end_pos = cols[2]
            else:
                end_pos = cols[1]

            if clip_dict.get(my_id):
                if stop_pos_dict.get("|".join([chrom, end_pos, strand])):
                    old_list = stop_pos_dict.get("|".join([chrom, end_pos, strand]))
                    old_list.append(cols[3])
                else:
                    stop_pos_dict["|".join([chrom, end_pos, strand])] = [cols[3]]

    return stop_pos_dict, clip_dict


def main(args):
    use_bam = False
    # set up temporary directory
    tmp_dir = subprocess.run(
        ["mktemp", "-d", "TEMP_correctMinimap_XXXXXX"], capture_output=True, text=True
    ).stdout.strip()

    # Filter BED for internal exons
    print("→ Filtering internal exons")
    filter_internal_exons(args.toga, tmp_dir)

    # Detect whether the input is SAM or BAM --> use magic number for gzipped files --> 1f 8b
    with open(args.minimap, "rb") as test_f:
        if test_f.read(2) == b"\x1f\x8b":
            print(
                "→ Detected your input to be BAM - final output will be compressed as well"
            )
            subprocess.check_output(
                f"samtools view -h {args.minimap} > {tmp_dir}/in.converted.sam",
                shell=True,
            )
            args.minimap = f"{tmp_dir}/in.converted.sam"
            use_bam = True

    # SAM to bed
    print(
        "→ Converting SAM to BED... depending on the size of your file, this may take a minute"
    )
    # Convert SAM to BAM (handy coordinates)
    subprocess.run(
        f"convert2bed --input=sam < {args.minimap} > {tmp_dir}/minimap.sam.bed",
        shell=True,
    )
    # Trim the file down
    subprocess.check_output(
        f"cut -d$'\t' -f1-8 {tmp_dir}/minimap.sam.bed > {tmp_dir}/minimap.sam.bed8",
        shell=True,
    )
    # read in stop positions and clip sequences
    stop_pos_dict, clip_dict = parse_sam_into_dicts(
        args.minimap, tmp_dir, args.clip_cutoff, wiggle=args.wiggle
    )

    # BED parsing
    matching_transcript_ids = []
    load_by_id = dict()
    exons_by_id = dict()
    matching_tuples = []

    # Load indexed exons
    with open(f"{tmp_dir}/indexed_exons.bed", "r") as all_file:
        for line in all_file:
            cols = line.split("\t")
            exons_by_id[cols[3]] = line

    print("→ Matching reference exons ")
    wiggle_report = []
    # We need to pull the missing bases from the reference genome later --> track this
    matched_with_wiggle = dict()
    with open(f"{tmp_dir}/internal_exons.bed", "r") as exon_file:
        for line in exon_file:
            cols = line.split("\t")
            strand = cols[5].strip()
            chrom = cols[0]
            my_id = cols[3]
            if strand == "+":
                stop = int(cols[2])
            else:
                stop = int(cols[1])

            # we use dicts as a fast and straight-forward coordinate matching approach
            for i in range(-args.wiggle, args.wiggle + 1):
                key = "|".join([chrom, str(stop + i), strand])
                if stop_pos_dict.get(key):
                    if i != 0 and args.debug:
                        wiggle_report.append(
                            "\t".join(
                                [my_id, str(stop_pos_dict.get(key)), str(i) + "\n"]
                            )
                        )
                    if i != 0:
                        for transcript in stop_pos_dict.get(key):
                            matched_with_wiggle[transcript] = True
                    matching_transcript_ids.extend(stop_pos_dict.get(key))
                    # if we matched with a wiggle, adjust the exon coords accordingly
                    matching_tuples.append(
                        (
                            stop_pos_dict.get(key),
                            f"{chrom}\t{int(cols[1])}\t{int(cols[2])}",
                            strand,
                            my_id,
                            i,
                        )
                    )

    # additional wiggle information
    if args.debug:
        with open(f"{tmp_dir}/wiggle_report.tsv", "w") as report_file:
            report_file.writelines(wiggle_report)

    """# if we have a wiggle match, we need to pull additional bases from the reference genome
    # We solve this by creating a BED file and feeding it to twoBitToFa
    with open(f"{tmp_dir}/minimap.wiggled.bed8", "w") as wiggle_file:
        with open(f"{tmp_dir}/minimap.sam.bed8", "r") as file:
            for line in file:
                cols = line.split("\t")
                my_id = cols[3]
                load_by_id[my_id] = line

                if matched_with_wiggle.get(my_id):
                    wiggle_file.write("\t".join(cols[:6]) + "\n")"""

    hits = []
    for my_tuple in matching_tuples:
        for transcript in my_tuple[0]:
            hits.append(
                f"{transcript}\t{my_tuple[1]}\t{my_tuple[2]}\t{my_tuple[3]}\t{clip_dict[transcript]}\t{my_tuple[4]}"
            )

    used_exons = dict()
    with open(f"{tmp_dir}/exon_precursor.bed", "w") as file:
        for my_tuple in matching_tuples:
            for transcript in my_tuple[0]:
                upstream_id = my_tuple[3]
                components = upstream_id.split("#")
                old_count = int(components[-1])
                components = components[:-1]
                wiggle_val = abs(args.wiggle)

                if my_tuple[2] == "+":
                    components.extend(["#", str(old_count + 1)])
                elif my_tuple[2] == "-":
                    components.extend(["#", str(old_count - 1)])
                else:
                    print(f"Expected +/- but got {my_tuple[2]}")
                    exit(1)

                target_id = "".join(components)
                target_exon = exons_by_id[target_id]
                cols = target_exon.split("\t")
                cols[3] = f"{cols[3]}"

                if not used_exons.get(cols[3]):
                    file.write("\t".join(cols))
                    used_exons[f"{cols[3]}_"] = True

    # TOGA exon sequences
    print("→ Fetching Exon sequences")
    subprocess.check_output(
        f"twoBitToFa /projects/hillerlab/genome/gbdb-HL/{args.assembly}/{args.assembly}.2bit -bed={tmp_dir}/exon_precursor.bed {tmp_dir}/exon_seqs.fa",
        shell=True,
    )
    exons_seqs_dict = parse_named_fasta(f"{tmp_dir}/exon_seqs.fa")

    hits_list = []
    for hit in hits:
        cols = hit.split("\t")
        if cols[4].strip() == "+":
            exon_id = cols[5].split("#")
            exon_id[-1] = f"{int(exon_id[-1]) + 1}"
            exon_id = "#".join(exon_id)
            wiggle_shift = args.wiggle - int(cols[7])
            # with a forward wiggle, the last base(s) before the 3'clip actually belong to the next exon
            # therefore, move the frame
            clip_seq = cols[6].strip()[wiggle_shift:]
            wiggle_oriented = int(cols[7])
        else:
            exon_id = cols[5].split("#")
            exon_id[-1] = f"{int(exon_id[-1]) - 1}"
            exon_id = "#".join(exon_id)
            wiggle_shift = args.wiggle + int(cols[7])
            # with a forward wiggle, the last base(s) before the 3'clip actually belong to the next exon
            # therefore, move the frame
            clip_seq = rev_comp(cols[6].strip())[wiggle_shift:]
            wiggle_oriented = -1 * int(cols[7])
        hits_list.append(
            {
                "origin": cols[0],
                "toga_exon": exon_id,
                "clip_seq": clip_seq,
                "toga_seq": exons_seqs_dict[exon_id],
                "strand": cols[4],
                "wiggle_shift": wiggle_oriented,
            }
        )

    # Match starting seqs
    print("→ Matching Clip and Exon Seqs")
    hit_ids = dict()
    hits_counter = 0
    with open(f"{tmp_dir}/hits.filtered.tsv", "w") as file:
        for hit in hits_list:
            # compare sequences and set flags accordingly
            if len(hit["clip_seq"]) > args.clip_cutoff:
                if (
                    hit["clip_seq"].upper()
                    == hit["toga_seq"][: len(hit["clip_seq"])].upper()
                ):
                    file.write(
                        "\t".join([str(hit[key]) for key in hit.keys()])
                        + "\tTRUE"
                        + "\n"
                    )
                    hit_ids[hit["origin"]] = hit["clip_seq"]
                    hits_counter += 1
                    # end iteration early
                    continue
            # if iteration is still going --> mismatch
            file.write(
                "\t".join([str(hit[key]) for key in hit.keys()]) + "\tFALSE" + "\n"
            )
    print(f"→ Matching sequences in {hits_counter} / {len(hits_list)}")

    # classifier SAMs for coordinate match
    print("→ Producing classifier SAMs")

    # get SAM headers
    subprocess.check_output(
        f"restore_sam_seq_headers {args.assembly} > {tmp_dir}/matchExonPos.matchDownstreamSeq.sam",
        shell=True,
    )
    # append to SAM file
    with open(args.minimap, "r") as in_sam:
        with open(f"{tmp_dir}/matchExonPos.matchDownstreamSeq.sam", "a") as out_sam:
            for line in in_sam:
                # exclude headers
                if not line.startswith("@"):
                    cols = line.split("\t")
                    if hit_ids.get(cols[0]):
                        cols[0] = f"{cols[0]}_{hit_ids.get(cols[0])}"
                        out_sam.write("\t".join(cols))

    # SAM CIGAR STRING CORRECTION
    print("→ Loading additional data")
    map_sam_to_exon = dict()
    dict_sam_exons = dict()
    with open(f"{tmp_dir}/hits.filtered.tsv", "r") as hits_file:
        for line in hits_file:
            cols = line.split("\t")
            if cols[-1].strip() == "TRUE":
                if not dict_sam_exons.get(cols[0]):
                    dict_sam_exons[cols[0]] = [cols[1]]
                else:
                    old_list = dict_sam_exons.get(cols[0])
                    old_list.append(cols[1])
                    dict_sam_exons[cols[0]] = old_list
                map_sam_to_exon[f"{cols[0]}&{cols[1]}"] = {
                    "exon_id": cols[1],
                    "clip_seq": cols[2],
                    "exon_seq": cols[3],
                    "wiggle": int(cols[5]),
                }

    indexed_exon_coords = dict()
    with open(f"{tmp_dir}/indexed_exons.bed", "r") as exon_file:
        for line in exon_file:
            cols = line.split("\t")
            indexed_exon_coords[cols[3]] = {"start": int(cols[1]), "stop": int(cols[2])}

    sam_indexed = dict()
    with open(args.minimap, "r") as sam_file:
        for line in sam_file:
            # skip headers
            if not line.startswith("@"):
                cols = line.split("\t")
                sam_indexed[cols[0]] = line

    sam_stops_indexed = dict()
    with open(f"{tmp_dir}/minimap.sam.bed8", "r") as sam_bed:
        for line in sam_bed:
            cols = line.split("\t")
            if dict_sam_exons.get(cols[3]):
                sam_stops_indexed[cols[3]] = {
                    "start": int(cols[1]),
                    "stop": int(cols[2]),
                    "strand": cols[5].strip(),
                }

    # <---- discarded.sam ---->
    discarded_used_ids = dict()
    subprocess.check_output(
        f"restore_sam_seq_headers {args.assembly} > {tmp_dir}/discarded.sam", shell=True
    )
    if args.debug:
        with open(f"{tmp_dir}/hits.filtered.tsv", "r") as hits_file:
            with open(f"{tmp_dir}/discarded.sam", "a") as discarded_file:
                for line in hits_file:
                    cols = line.split("\t")
                    if cols[-1].strip() == "FALSE" and not discarded_used_ids.get(
                        cols[0]
                    ):
                        # get original SAM entry
                        original_sam = sam_indexed[cols[0]]
                        discarded_file.write(original_sam)
                        discarded_used_ids[cols[0]] = True

    # Build new SAM with edited CIGAR strings
    print("→ Building output SAM")
    subprocess.check_output(
        f"restore_sam_seq_headers {args.assembly} > {tmp_dir}/corrected.sam", shell=True
    )
    tracked_cigars = dict()
    discarded_cigars = []
    with open(f"{tmp_dir}/corrected.sam", "a") as out_file:
        # for all eligible SAMs...
        for sam_id in dict_sam_exons.keys():
            # ... walk through all matching exons
            for exon_index, exon in enumerate(dict_sam_exons[sam_id]):
                joined_id = f"{sam_id}&{exon}"
                sam_coords = sam_stops_indexed[sam_id]
                exon_coords = indexed_exon_coords[exon]
                sam_cols = sam_indexed[sam_id].split("\t")
                my_cigar = Cigar(sam_cols[5])
                clip_seq = map_sam_to_exon[joined_id]["clip_seq"]

                # skip clips which are too short
                if len(clip_seq) < args.clip_cutoff:
                    continue

                # (+) strand
                if sam_coords["strand"] == "+":
                    wiggle_val = map_sam_to_exon[joined_id]["wiggle"]
                    dist = exon_coords["start"] - sam_coords["stop"]
                    # dist += map_sam_to_exon[joined_id]['wiggle']
                    # trim upstream exon match
                    my_cigar.features[-2].size = my_cigar.features[-2].size - wiggle_val
                    # insert gap
                    my_cigar.insert_feature(CigarFeature("N", dist + wiggle_val), -1)
                    # insert former clip
                    my_cigar.insert_feature(CigarFeature("=", len(clip_seq)), -1)
                    # shorten softclip
                    my_cigar.features[-1].size = (
                        my_cigar.features[-1].size - len(clip_seq) + wiggle_val
                    )
                    # BED 1-open
                    if map_sam_to_exon[joined_id]["wiggle"] < 0:
                        dist += 1

                    sam_cols[5] = str(my_cigar)
                # (-) strand
                else:
                    # flip to convert to (+) coordinates
                    wiggle_val = -1 * map_sam_to_exon[joined_id]["wiggle"]
                    dist = sam_coords["start"] - exon_coords["stop"]
                    # trim upstream exon match
                    my_cigar.features[1].size = my_cigar.features[1].size - abs(
                        wiggle_val
                    )
                    # insert gap
                    my_cigar.insert_feature(
                        CigarFeature("N", dist + abs(wiggle_val)), 1
                    )
                    # insert former clip
                    my_cigar.insert_feature(CigarFeature("=", len(clip_seq)), 1)
                    # shorten softclip
                    my_cigar.features[0].size = (
                        my_cigar.features[0].size - len(clip_seq) + abs(wiggle_val)
                    )
                    # BED 1-open
                    if map_sam_to_exon[joined_id]["wiggle"] < 0:
                        dist += 1
                    # Exon frame to compare against
                    exon_seq = map_sam_to_exon[joined_id]["exon_seq"][
                        : len(clip_seq)
                    ].upper()
                    # iterate over the clipped seqs and track (mis)matches
                    for i, base in enumerate(clip_seq):
                        if base != exon_seq[i]:
                            print(f"{sam_id} has mismatches")

                    sam_cols[5] = str(my_cigar)
                    sam_cols[3] = str(int(sam_cols[3]) - dist - len(clip_seq))

                # now, update the 3' clip information
                id_vals = sam_cols[0].split("_")
                for i, value in enumerate(id_vals):
                    if value.startswith("3Clip"):
                        id_vals[i] = "3Clip0"
                sam_cols[0] = "_".join(id_vals)

                # create new entry
                if not tracked_cigars.get(sam_id):
                    tracked_cigars[sam_id] = [str(my_cigar)]
                    # sam_cols[0] = f"{sam_cols[0]}#0"
                    out_file.write("\t".join(sam_cols))
                else:
                    # check if an identical cigar string exists for this transcript
                    old_list = tracked_cigars[sam_id]
                    if str(my_cigar) in old_list:
                        discarded_cigars.append("\t".join(sam_cols))
                        continue
                    else:
                        # save next cigar
                        bitwise_flag = 0
                        if strand == "-":
                            bitwise_flag += 16

                        # The CIGAR string is non redudant, but another correction already exists
                        # --> set bitwise flag to indicate secondary alignments
                        # TODO: Add this to the samtypes library
                        bitwise_flag += 256
                        old_list.append(str(my_cigar))
                        tracked_cigars[sam_id] = old_list
                        # sam_cols[0] = f"{sam_cols[0]}#{len(old_list) - 1}"
                        sam_cols[1] = str(bitwise_flag)
                        out_file.write("\t".join(sam_cols))

    print(f"→ Generating collapse report to temp directory")
    total_unique = 0
    for key in tracked_cigars.keys():
        if len(tracked_cigars[key]) > 1:
            total_unique += 1

    print(f"A total of {total_unique} transcripts have unique correction variants")
    if args.debug:
        with open(f"{tmp_dir}/collapsed.sam", "w") as collapse_file:
            for line in discarded_cigars:
                collapse_file.write(line)

    # <--- Generating a combined output file --->
    # Append all uncorrected sam rows from the input file to the corrected file
    # Distinguish by SAM IDs
    print(f"→ Merging corrected and input SAM files")
    with open(f"{tmp_dir}/corrected.sam", "a") as corr_file:
        with open(f"{args.minimap}", "r") as original_file:
            for line in original_file:
                if line.startswith("@"):
                    continue
                else:
                    sam_id = line.split("\t")[0]
                    if sam_id not in dict_sam_exons.keys():
                        corr_file.write(line)

    # Convert to BAM
    if use_bam:
        subprocess.check_output(
            f"samtools view -S -b {tmp_dir}/corrected.sam > {args.output}", shell=True
        )
    else:
        subprocess.check_output(f"cp {tmp_dir}/corrected.sam {args.output}", shell=True)
    # Clean up
    if not args.keep_temp:
        print(f"→ Deleting temporary directory at {tmp_dir}")
        shutil.rmtree(f"{tmp_dir}")
    print(f"→ The final corrected BAM is here: {args.output}")


if __name__ == "__main__":
    # parse args
    arguments = parser.parse_args()
    main(arguments)


class CigarFeature:
    def __init__(self, symbol: str, size: int):
        self.symbol = symbol
        self.size = size

    def set_size(self, size: int):
        self.size = size

    def get_size(self):
        return self.size

    def __str__(self):
        return f"{self.size}{self.symbol}"


class Cigar:
    def __init__(self, value_string):
        new_features = []
        cig_iter = groupby(value_string, lambda c: c.isdigit())
        for g, n in cig_iter:
            next_feature = CigarFeature(
                size=int("".join(n)), symbol="".join(next(cig_iter)[1])
            )
            new_features.append(next_feature)
        self.features = new_features

    def from_string(self, value_string: str):
        new_features = []
        print(value_string)
        if value_string == "*":
            yield (0, None)
            raise StopIteration
        cig_iter = groupby(value_string, lambda c: c.isdigit())
        for g, n in cig_iter:
            new_features.append(CigarFeature(str(n), next(cig_iter)[1]))
        print(new_features)
        self.features = new_features

    def get_features(self):
        return self.features

    def insert_feature(self, feature: CigarFeature, index=-1):
        features = self.features
        features.insert(index, feature)
        self.features = features

    def __str__(self):
        final_str = ""
        for feature in self.features:
            final_str += f"{feature.size}{feature.symbol}"
        return final_str
