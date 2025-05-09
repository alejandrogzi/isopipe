#!/beegfs/projects/hillerlab/genome/src/ORFTree/.venv/bin/python3

import os
import shutil
import subprocess
import logging
import pandas as pd
import argparse
from joblib import load
from Bio import SeqIO
from sklearn.ensemble import RandomForestClassifier
from utils.merging import filter_by_relative_score_strict
from utils import merging
from utils.common import write_results, export_df_orfs, map_bed
from utils.blast_prep import blast_main
from utils.translationai_prep import translationai_main
from bedutils.bed import BedRow
from utils.overrules import toga_overrule_a, toga_overrule_b
from typing import Tuple


MODEL_PATH = "/beegfs/projects/hillerlab/genome/src/ORFTree/model.joblib"
QUERY_ANNOTATION = "query_annotation.bed"
CODONS = "selenocysteine_codons.tsv"
TX_META = "transcripts_meta.tsv.gz"
FEATURES = [
    "blast_pid",
    "blast_e-value",
    "blast_offset",
    "blast_percentage_aligned",
    "toga_pid",
    "toga_blosum",
    "toga_overlap_bp",
    "translationAI_orf_start",
    "translationAI_orf_stop",
]
TOOLS = [
    "TransDecoder.LongOrfs",
    "blastp",
    "diamond",
    "run_translation_ai",
    "/beegfs/projects/hillerlab/genome/src/ORFTree/.venv/bin/orfipy",
]

logger = logging.getLogger(__name__)


def main(args):
    tmp_dir = subprocess.run(
        ["mktemp", "-d", "TEMP_ORF_Tree_XXXXXX"], capture_output=True, text=True
    ).stdout.strip()

    logging.basicConfig(
        level=logging.INFO,
        format="{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M",
    )

    __check_packages(args.use_blast)

    __run_blast(tmp_dir)
    __run_tai(args, tmp_dir)

    df = __load_blast(tmp_dir)

    (tx_meta, bed, codons, df) = __prepare_toga(args, tmp_dir, df)
    df = __run_toga(tx_meta, bed, codons, df, tmp_dir)

    __predict(df, tmp_dir)

    logger.info("Done")


def __check_packages(use_blast: bool) -> None:
    """
    Check if all required packages are installed

    Parameters
    ----------
    use_blast : bool
        If True, check for BLASTP, otherwise check for DIAMOND

    Returns
    -------
    None

    Example
    -------
    >>> __check_packages(use_blast=True)
    """
    logger.info("Looking for tools")
    logger.info("Checking tool dependencies")

    for tool in TOOLS:
        if not shutil.which(tool):
            if (
                tool == "diamond"
                and args.use_blast
                or tool == "diamond"
                and not args.use_blast
            ):
                logger.warning(f"{tool}\tNot found but is not required")
            else:
                logger.error(f"{tool}\tNot found!")
                exit(-1)
        else:
            logger.info(f"Will use {tool} -> {shutil.which(tool)}")


def __run_toga(
    toga_results: str,
    toga_query_annotation: str,
    toga_masked: str,
    df: pd.DataFrame,
    tmp_dir: str,
) -> pd.DataFrame:
    # TOGA
    logger.info("Loading TOGA data")
    inserts = merging.load_toga_two(
        toga_results,
        toga_query_annotation,
        toga_masked,
        f"{tmp_dir}/TOGA/overlap.bed",
    )
    df = df.merge(
        inserts,
        left_on=["canonical_id", "genomic_coords"],
        right_on=["canonical_id", "genomic_coords"],
        how="left",
        right_index=False,
    )
    df.set_index(["genomic_coords", "canonical_id"], inplace=True)
    df["toga_overlap_bp"] = df["toga_overlap_bp"].fillna(-1)
    logger.info(f"Dataset size is now {len(df)} rows")

    # Clean up
    df["toga_masked"] = df["toga_masked"].fillna(False)
    for x in ["blast_id", "toga_id", "translationAI_id", "toga_label"]:
        df[x] = df[x].fillna("X")

    df["blast_nested"] = df["blast_nested"].fillna(-1)
    df = df.fillna(0)
    # clean duplicates
    df.drop_duplicates(keep="first", inplace=True)
    # Use this as a checkpoint file in case the predictions crashes
    # logger.info(f"Writing to {tmp_dir}/mergedData.tsv")
    # df.to_csv(f"{tmp_dir}/mergedData.tsv", sep="\t")

    return df


def __prepare_toga(
    args: argparse.Namespace, tmp_dir: str, df: pd.DataFrame
) -> Tuple[str, str, str, pd.DataFrame]:
    """
    Prepare TOGA data for the analysis

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments
    tmp_dir : str
        Temporary directory to store the results
    df : pd.DataFrame
        DataFrame with the input data

    Returns
    -------
    Tuple[str, str, str]
        Tuple with the paths to the TOGA annotation, codons and tx_meta files

    Example
    -------
    >>> __prepare_toga(args, "/tmp", df)
    ("/tmp/TOGA/query_annotation.bed", "/tmp/TOGA/selenocysteine_codons.tsv", "/tmp/TOGA/transcripts_meta.tsv.gz")
    """
    annotation = f"{args.toga}/{QUERY_ANNOTATION}"
    codons = f"{args.toga}/meta/{CODONS}"
    tx_meta = f"{args.toga}/meta/{TX_META}"

    # We expect pre-calculated TOGA inputs and match them against the alignments to computer reciprocal overlaps of
    # at least 80%
    logger.info("Entering TOGA module")
    os.mkdir(f"{tmp_dir}/TOGA")
    df = df.reset_index()
    ex = export_df_orfs(df, args.alignments, f"{tmp_dir}/TOGA/orf_to_overlap.bed")
    logger.info(
        f"Matching {ex[0]} potential ORFs against TOGA (this can take a moment)"
    )
    subprocess.run(
        f"bedtools intersect -a {annotation} -b {tmp_dir}/TOGA/orf_to_overlap.bed -wo -f 0.8 -r -s > {tmp_dir}/TOGA/overlap.bed",
        shell=True,
    )
    logger.info("Exiting TOGA module")

    return (tx_meta, annotation, codons, df)


def __run_blast(tmp_dir: str) -> None:
    """
    Run BLASTP on the input data

    Parameters
    ----------
    tmp_dir : str
        Temporary directory to store the results

    Returns
    -------
    None

    Example
    -------
    >>> run_blast("/tmp")
    """
    logger.info("INFO: Entering BLAST module")
    blast_main(
        tmp_dir,
        args.fasta,
        args.alignments,
        args.blastdb,
        use_blast=args.use_blast,
        use_nextflow=True,
    )

    logger.info("INFO: Exiting BLAST module")


def __load_blast(tmp_dir: str) -> pd.DataFrame:
    """
    Load BLASTP data and merge it with the input data

    Parameters
    ----------
    tmp_dir : str
        Temporary directory to store the results

    Returns
    -------
    pd.DataFrame
        Merged DataFrame with BLASTP data

    Example
    -------
    >>> df = __load_blast("/tmp")
    """
    logger.info("Loading BLAST data")
    df = merging.load_blast(
        f"{tmp_dir}/ORF/longest_orfs.pep.nestedORFs.fa",
        f"{tmp_dir}/BLASTP/longest_orfs.pep.nestedORFs.fa.fmt6",
        f"{tmp_dir}/BLASTP/trimmed_alignments.bed",
    )
    df.set_index(["genomic_coords", "canonical_id"], inplace=True)
    logger.info(f"Dataset size is now {len(df)} rows")

    # TranslationAI
    logger.info("Loading TranslationAI data")
    inserts = merging.load_translation_ai(
        f"{tmp_dir}/TranslationAI/out.translationai.tsv",
        f"{tmp_dir}/TranslationAI/translationai.trimmed.bed",
    )

    df = df.reset_index().merge(
        inserts,
        left_on=["genomic_coords", "canonical_id"],
        right_on=["translationAI_coords", "translationAI_id"],
        how="outer",
    )
    df.loc[df["canonical_id"].isnull(), "canonical_id"] = df["translationAI_id"]
    df.loc[df["genomic_coords"].isnull(), "genomic_coords"] = df["translationAI_coords"]
    df.set_index(["genomic_coords", "canonical_id"], inplace=True)

    logger.info(f"Dataset size is now {len(df)} rows")

    return df


def __run_tai(args, tmp_dir: str) -> None:
    """
    Run TranslationAI on the input data

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments
    tmp_dir : str
        Temporary directory to store the results

    Returns
    -------
    None

    Example
    -------
    >>> run_tai(args, "/tmp")
    """
    logger.info("INFO: Entering TranslationAI module")
    translationai_main(args.fasta, args.alignments, tmp_dir, use_nextflow=True)
    logger.info("INFO: Exiting TranslationAI module")


def __predict(df: pd.DataFrame, tmp_dir: str) -> None:
    """
    Run the prediction module

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with the input data
    tmp_dir : str
        Temporary directory to store the results

    Returns
    -------
    None
    """
    logger.info("INFO: Entering prediction module")

    # First: load the model
    try:
        with open(MODEL_PATH, "rb") as f:
            clf: RandomForestClassifier = load(f)
    except Exception as err:
        logger.error(f"Failed to load model: {err}")
        exit(1)

    logger.info(f"Finished loading model ({MODEL_PATH})")

    x = df.loc[:, FEATURES]
    probs = clf.predict_proba(x)
    logger.info("Finished predictions")
    df["toga_overrule"] = False
    df["class0_probability"] = probs[:, 0]
    df["class1_probability"] = probs[:, 1]
    df.reset_index(inplace=True)

    df.to_csv(f"{args.output_dir}/predictions_uncollapsed.tsv", sep="\t")

    # Drop extremly redudant rows
    logger.info("Applying BLAST nested cutoff")
    df = df.loc[
        (
            ((1 - df["blast_nested_offset"] >= 0.7) | (df["blast_nested"] <= 4))
            | (
                (df["translationAI_id"] != "X")
                | (df["toga_coords"] == df["genomic_coords"])
            )
            | (df["blast_nested"] == -1)
        )
    ]
    logger.info(f"Applied BLAST nested cutoff, dataset reduced to {len(df)} rows")

    if args.toga_overrule:
        logger.info("Applying TOGA overrules")
        df = df.reset_index()

        sequences = dict()
        with open(args.fasta) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences[record.id] = record.seq

        unique_canonic = df["canonical_id"].unique()
        grouped = df.groupby("canonical_id")

        transcripts = dict()
        with open(args.alignments, "r") as f:
            for line in f:
                row = BedRow(line)
                transcripts[row.id_str] = row

        toga_rows = dict()
        with open(args.toga_query_annotation, "r") as f:
            for line in f:
                row = BedRow(line)
                toga_rows[row.id_str] = row

        overrule_ids = toga_overrule_a(
            unique_canonic, grouped, transcripts, toga_rows, sequences
        )
        overrule_ids.extend(
            toga_overrule_b(unique_canonic, grouped, transcripts, toga_rows, sequences)
        )

        overrule_rows = []
        for x in overrule_ids:
            x["class1_probability"] = x["toga_pid"] / 100
            x["toga_overrule"] = True
            if str(x["genomic_coords"]) != "0":
                overrule_rows.append(x)

    df.to_csv(f"{args.output_dir}/predictions_raw.tsv", sep="\t")

    logger.info("Trimming secondary ORFs (this will take some time)")
    # Apply threshold
    df = df.loc[
        (df["class1_probability"] >= args.threshold) | (df["toga_overrule"] == True)
    ]
    # Apply secondary ORFs
    df = df.groupby("canonical_id", group_keys=False).apply(
        filter_by_relative_score_strict
    )
    logger.info(f"Trimmed secondary ORF candidates, dataset reduced to {len(df)} rows")

    if args.toga_overrule:
        logger.info(f"Appending ORF predictions from TOGA overrules")
        df = pd.concat([df, pd.DataFrame(overrule_rows)], ignore_index=True)

    logger.info(f"Writing un-collapsed output BED to {args.output_dir}")

    # Drop useless columns
    df.drop([0, 2])
    # Drop duplicates -> only keep the better one
    df.sort_values(
        by=["canonical_id", "class1_probability"], inplace=True, ascending=False
    )
    df.drop_duplicates(inplace=True, subset=["canonical_id", "genomic_coords"])
    ranks = write_results(
        df, args.threshold, f"{args.output_dir}/predictions", with_ground_label=False
    )

    # Store the rank of a prediction compared to its peers
    df["rank"] = ranks
    df.to_csv(f"{args.output_dir}/predictions.tsv", sep="\t")

    logger.info(f"Writing (non-)coding BEDs to output dir to {args.output_dir}")
    map_bed(df, args.alignments, args.output_dir)
    logger.info("BED files written")

    if not args.keep_temp:
        logger.info(f"Cleaning temporary directory {tmp_dir}")
        shutil.rmtree(tmp_dir)

    return


def parse() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    -------
    argparse.Namespace

    Example
    -------
    >>> parse()
    """
    parser = argparse.ArgumentParser(description="Hillerlab ORF pipeline!")

    parser.add_argument(
        "--fasta", type=str, metavar="Path to transcript sequences [.fa]", help="Transcript FASTA file"
    )
    parser.add_argument(
        "--alignments",
        type=str,
        metavar="Path to reads aligned to query genome [.bed]",
        help="Transcript alignments in BED format",
    )
    parser.add_argument("--blastdb", type=str, metavar="Path to BLAST database", help="BLAST database")

    # INFO: Unifiying TOGA-dependent flags into one
    parser.add_argument("--toga", type=str, metavar="Path to TOGA directory", help="TOGA results")

    parser.add_argument(
        "--output_dir", type=str, metavar="path", help="Output directory"
    )
    parser.add_argument(
        "--toga_overrule",
        action="store_true",
        help="Overrule with confident TOGA predictions",
    )
    parser.add_argument(
        "--keep_temp", action="store_true", help="Do not clean the temp directory"
    )
    parser.add_argument(
        "--use_blast", action="store_true", help="Use BLASTP rather than diamond"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.03,
        help="Use a non-default threshold for classification",
    )
    parser.add_argument(
        "--no-reduce",
        action="store_true",
        help="Use a non-default threshold for classification",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse()
    main(args)
