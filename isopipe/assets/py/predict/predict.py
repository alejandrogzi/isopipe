#!/usr/bin/env python3

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.3"

import argparse
import logging
from os import PathLike
from pathlib import Path
from typing import Dict, List, Union

import numpy as np
import pandas as pd
import xgboost as xgb

MODEL = "model/orf_classifier_model.json"
MODEL_PATH = Path(__file__).parent / MODEL
UNIQUE_TAI_TAG = "UT"
START_CODONS = ["ATG", "CTG"]
STOP_CODONS = ["TAA", "TAG", "TGA"]
BLAST_COLS: List = [
    "chr",
    "start",
    "end",
    "id",
    "strand",
    "tai_start_score",
    "tai_stop_score",
    "relative_orf_start",
    "relative_orf_end",
    "start_codon",
    "stop_codon",
    "inner_stops",
    "orf_type",
    "blast_pid",
    "blast_evalue",
    "blast_offset",
    "blast_length",
    "blast_percentage_aligned",
]
SAMBA_COLS: List = ["prefix", "rna_score"]
TAI_MASKING_NAN_COLS: List = [
    "tai_start_score",
    "tai_stop_score",
    "tai_mean_score",
    "tai_combined_score",
]
BLAST_MASKING_NAN_COLS: List = [
    "blast_pid",
    "neg_log10_blast_evalue",
    "blast_offset",
    "blast_percentage_aligned",
]
FEATURES: List = [
    "tai_start_score",
    "tai_stop_score",
    "inner_stops",
    "orf_type",
    "blast_pid",
    "neg_log10_blast_evalue",
    "blast_offset",
    "blast_percentage_aligned",
    "rna_score",
    "tai_mean_score",
    "blast_match",
    "log_orf_len",
    "has_canonical_start",
    "has_canonical_stop",
]
ORF_TYPE_MAPPING: Dict = {
    "CO": 1,
    "CN": 2,
    "UN": 3,
    "FP": 4,
    "TP": 5,
    "TN": 6,
    "FN": 7,
}

log = logging.getLogger(__name__)
logging.basicConfig(encoding="utf-8", level=logging.INFO)


def run(args: argparse.Namespace) -> None:
    """
    Executes the main workflow for loading a model, performing predictions,
    and saving the results.

    This function orchestrates the process of loading a RandomForestClassifier
    model, reading and preparing input data from BLAST, TAI, and TOGA files,
    generating predictions, and then saving the final table to a TSV file.
    It also includes a commented-out section for an 'overrule' function,
    indicating potential future functionality.

    Parameters
    ----------
    args : argparse.Namespace
        An object containing command-line arguments, expected to have attributes
        `model` (path to model), `blast` (path to BLAST file),
        `tai` (path to TAI file), `toga` (path to TOGA file), and
        `outdir` (output directory).

    Returns
    -------
    None

    Example
    -------
    >>> import argparse
    >>> # Create a dummy Namespace for demonstration
    >>> # dummy_args = argparse.Namespace(
    >>> #     model="my_model.joblib",
    >>> #     blast="blast.tsv",
    >>> #     tai="tai.tsv",
    >>> #     toga="toga.tsv",
    >>> #     outdir="output_data",
    >>> #     # overrule=False # Include if overrule is an argument
    >>> # )
    >>> # # Assuming necessary files exist and a model can be loaded
    >>> # run(dummy_args)
    >>> # # This would create a 'predictions.tsv' file in 'output_data'
    """
    model = xgb.XGBClassifier()
    model.load_model(args.model)

    log.info("INFO: Model loaded successfully!")
    log.info(f"INFO: Number of features: {model.n_features_in_}")

    table = predict(args.blast, args.samba, model)

    _ = map_to_blocks(
        table,
        args.alignments,
        args.outdir,
        args.prefix,
        args.min_score_max_predictions,
        args.max_predictions,
    )


def map_to_blocks(
    table: pd.DataFrame,
    alignments: Union[str, PathLike],
    outdir: Union[str, PathLike],
    prefix: str,
    min_score_max_predictions: float,
    max_predictions: int,
) -> None:
    """
    Maps prediction results to genomic blocks based on alignment data and saves them.

    This function takes a DataFrame of prediction results, merges it with genomic
    alignment data (BED file format), modifies certain columns to represent block
    coordinates, and then saves both the updated prediction table and the new
    BED-formatted table to the specified output directory.

    Parameters
    ----------
    table : pd.DataFrame
        A pandas DataFrame containing prediction results, expected to have
        'blast_id', 'tai_orf_start', and 'tai_orf_end' columns.
    alignments : Union[str, PathLike]
        Path to the genomic alignments file (BED format).
    outdir : Union[str, PathLike]
        The output directory where the modified prediction table and
        the new BED file will be saved.
    prefix : str
        A string prefix to be added to the output filenames.

    Returns
    -------
    None

    Example
    -------
    >>> import pandas as pd
    >>> # Create dummy dataframes and paths for demonstration
    >>> # dummy_table = pd.DataFrame({
    >>> #     'blast_id': ['gene1__orfA', 'gene2__orfB'],
    >>> #     'tai_orf_start': [100, 500],
    >>> #     'tai_orf_end': [200, 600]
    >>> # })
    >>> # with open("dummy_alignments.bed", "w") as f:
    >>> #     f.write("chr1\\t10\\t1000\\tgene1__orfA\\t0\\t+\\n")
    >>> #     f.write("chr1\\t50\\t1200\\tgene2__orfB\\t0\\t+\\n")
    >>> #
    >>> # map_to_blocks(dummy_table, "dummy_alignments.bed", "./output", "test_")
    >>> # # This would create 'test_predictions.tsv' and 'test_predictions.bed'
    >>> # # in the './output' directory.
    """
    log.info(f"INFO: Initial size of table: {len(table)}")
    table = table.copy()
    # table["id"] = [id.split("__")[0] for id in table.blast_id] -> replace by prefix
    table["tag"] = ["#" + id.split("__")[1] for id in table.id]

    bed = pd.read_csv(alignments, sep="\t", header=None)
    merged = (
        bed.assign(prefix=bed[3].str.split("__").str[0])
        .merge(
            table[["prefix", "start", "end", "tag", "prob_coding"]],
            on="prefix",
            how="left",
        )
        .dropna(subset=["start"])
        .assign(
            start=lambda df: df["start"].astype(int),
            end=lambda df: df["end"].astype(int),
        )
    )

    merged[6] = merged["start"]
    merged[7] = merged["end"]
    merged[3] += merged["tag"]

    table = (
        table[table["prob_coding"] >= min_score_max_predictions]
        .sort_values("prob_coding", ascending=False)
        .groupby("prefix")
        .head(max_predictions)
        .reset_index(drop=True)
    )

    # INFO: add #DU tag to non-best predictions
    table["rank"] = table.groupby("prefix").cumcount() + 1
    table.loc[table["rank"] > 1, "id"] += "#DU"

    log.info(f"INFO: Final size of table: {len(table)}")
    log.info(f"INFO: Writing predictions to {args.outdir}/{prefix}.predictions.tsv")

    table.drop(columns=["prefix", "tag", "rank"]).to_csv(
        f"{args.outdir}/{prefix}.predictions.tsv", index=False, header=True, sep="\t"
    )

    if args.keep_raw:
        merged.drop(columns=["prefix", "start", "end", "tag"]).to_csv(
            f"{outdir}/{prefix}.all.predictions.bed",
            sep="\t",
            header=False,
            index=False,
        )

    merged = (
        merged[merged["prob_coding"] >= min_score_max_predictions]
        .sort_values("prob_coding", ascending=False)
        .groupby("prefix")
        .head(max_predictions)
        .reset_index(drop=True)
    )

    # INFO: add #DU tag to non-best predictions
    merged["rank"] = merged.groupby("prefix").cumcount() + 1
    merged.loc[merged["rank"] > 1, 3] += "#DU"

    merged.drop(
        columns=["prefix", "start", "end", "tag", "prob_coding", "rank"]
    ).to_csv(f"{outdir}/{prefix}.predictions.bed", sep="\t", header=False, index=False)

    return


def predict(
    blast: Union[str, PathLike],
    samba: Union[str, PathLike],
    model: xgb.XGBClassifier,
) -> pd.DataFrame:
    """
    Performs predictions using a loaded RandomForestClassifier model on input data.

    This function reads data from BLAST, TAI, and TOGA files, prepares a query
    table, and then uses the provided model to generate class probabilities.

    Parameters
    ----------
    blast : Union[str, PathLike]
        Path to the BLAST input file.
    tai : Union[str, PathLike]
        Path to the TAI input file.
    toga : Union[str, PathLike]
        Path to the TOGA input file.
    model : RandomForestClassifier
        The pre-trained RandomForestClassifier model used for prediction.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the merged input data along with two
        new columns: 'class0_prob' and 'class1_prob', which
        represent the prediction probabilities for each class.

    Example
    -------
    >>> # Assuming 'blast.tsv', 'tai.tsv', 'toga.tsv' exist and a model is loaded
    >>> # from sklearn.ensemble import RandomForestClassifier
    >>> # dummy_model = RandomForestClassifier() # Replace with your actual loaded model
    >>> # result_df = predict("blast.tsv", "tai.tsv", "toga.tsv", dummy_model)
    >>> # print(result_df.head())
    """
    table = read(blast, samba)
    query = table.loc[:, FEATURES]

    for col in query.columns:
        query[col] = pd.to_numeric(query[col], errors="coerce")

    predictions = model.predict(query)
    probabilities = model.predict_proba(query)
    log.info("INFO: Finished predictions!")

    table["predicted_class"] = predictions
    table["prob_noncoding"] = probabilities[:, 0]
    table["prob_coding"] = probabilities[:, 1]

    log.info("INFO: Predictions summary:")
    log.info(f"INFO:  Predicted as real ORFs (1): {(predictions == 1).sum()}")
    log.info(f"INFO:  Predicted as spurious (0): {(predictions == 0).sum()}")

    return table


def read(
    blast: Union[str, PathLike, Path],
    samba: Union[str, PathLike, Path],
) -> pd.DataFrame:
    """
    Reads and merges data from BLAST, TAI, and TOGA files into a single DataFrame.

    Parameters
    ----------
    blast : Union[str, PathLike, Path]
        Path to the BLAST input file.
    samba : Union[str, PathLike, Path]
        Path to the RNASamba input file.

    Returns
    -------
    pd.DataFrame
        A merged pandas DataFrame containing data from all three input files.

    Example
    -------
    >>> # Assuming 'blast.tsv', 'samba.tsv' exist
    >>> # merged_data = read("blast.tsv", "samba.tsv")
    """
    blast = read_blast(blast)
    samba = read_samba(samba)

    table = merge_tables(blast, samba)
    return table


def read_blast(path: Union[str, PathLike, Path]) -> pd.DataFrame:
    """
    Reads a BLAST tab-separated file and processes it into a pandas DataFrame.

    This function reads the file using predefined column names (BLAST_COLS)
    and generates a unique 'key' column based on 'blast_id', 'blast_orf_start',
    and 'blast_orf_end'.

    Parameters
    ----------
    path : Union[str, PathLike, Path]
        The file path to the BLAST input file.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the BLAST data with an additional 'key' column.

    Example
    -------
    >>> # Assuming 'my_blast.tsv' exists with appropriate tab-separated data
    >>> # blast_df = read_blast("my_blast.tsv")
    >>> # print(blast_df.head())
    """
    blast = pd.read_csv(path, sep="\t", header=None, names=BLAST_COLS)

    # INFO: needs to be the cannonical ID
    # R1_chr1__OR2#NE1 -> R1_chr1 [samba] + R1_chr1:1-10(+) [bed]
    blast["prefix"] = blast["id"].str.split("__").str[0]
    blast["key"] = (
        blast["prefix"]
        + ":"
        + blast["start"].astype(str)
        + "-"
        + blast["end"].astype(str)
        + "("
        + blast["strand"].astype(str)
        + ")"
    )
    return blast


def read_samba(path: Union[str, PathLike, Path]) -> pd.DataFrame:
    """
    Reads a RNASamba tab-separated file and processes it into a pandas DataFrame.

    This function reads the file using predefined column names (SAMBA_COLS)

    Parameters
    ----------
    path : Union[str, PathLike, Path]
        The file path to the RNASamba input file.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the RNASamba data with an additional 'key' column.

    Example
    >>> # Assuming 'my_samba.tsv' exists with appropriate tab-separated data
    >>> # samba_df = read_samba("my_samba.tsv")
    >>> # print(samba_df.head())
    """
    return pd.read_csv(path, sep="\t", header=None, names=SAMBA_COLS)


def merge_tables(
    blast: pd.DataFrame,
    samba: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merges the BLAST and RNASamba DataFrames into a final comprehensive table.

    The merging process involves:
    1. Inner merging BLAST and RNASamba DataFrames on the 'prefix' column.
    2. Extracting chromosome information into 'm_chr' from the merged key.
    3. Creating a 'toga_key' based on chromosome and ORF start/end, considering strand.
    4. Left merging with the TOGA DataFrame on 'toga_key' and TOGA's 'key'.
    5. Adding a binary flag 'toga_overlap_bp' indicating if 'toga_pid' is present.

    Parameters
    ----------
    blast : pd.DataFrame
        The DataFrame containing BLAST data.
    samba : pd.DataFrame
        The DataFrame containing RNASamba data.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame that is the result of merging all three input DataFrames,
        with additional derived columns like 'm_chr', 'toga_key', and 'toga_overlap_bp'.

    Example
    -------
    >>> # Assuming blast_df, samba_df are pre-read DataFrames
    >>> # final_table = merge_tables(blast_df, samba_df)
    >>> # print(final_table.head())
    """
    merged = blast.merge(samba, on="prefix", how="left")

    log.info(f"INFO: Merged df looks like this:\n{merged.head()}")
    log.info(f"INFO: BLAST has {len(blast.prefix.unique())} unique rows!")
    log.info(f"INFO: RNASamba has {len(samba)} rows!")
    log.info(f"INFO: BLAST and RNASamba files have {len(merged)} rows!")

    did_not_get_any_orf = samba[~samba.prefix.isin(blast.prefix)].sort_values(
        "rna_score", ascending=False
    )
    if did_not_get_any_orf.shape[0] > 0:
        log.info(
            f"INFO: These prefixes did not get any ORF prediction:\n{did_not_get_any_orf}"
        )

    not_catched_by_samba = blast[~blast.prefix.isin(samba.prefix)]
    if not_catched_by_samba.shape[0] > 0:
        log.info(
            f"INFO: These prefixes are not caught by samba: {not_catched_by_samba}"
        )

    if len(merged.prefix.unique()) != len(blast.prefix.unique()):
        raise ValueError(
            "ERROR: BLAST and RNASamba files do not match! Some BLAST rows do not have a prediction in RNASamba!"
        )

    merged["tai_mean_score"] = (merged.tai_start_score + merged.tai_stop_score) / 2
    merged["blast_match"] = [1 if x > 0 else 0 for x in merged.blast_percentage_aligned]
    merged["log_orf_len"] = np.log1p(
        merged.relative_orf_end - merged.relative_orf_start
    )
    merged["has_canonical_start"] = merged.start_codon.isin(START_CODONS).astype(int)
    merged["has_canonical_stop"] = merged.stop_codon.isin(STOP_CODONS).astype(int)
    merged["neg_log10_blast_evalue"] = -np.log10(merged.blast_evalue)

    mask = (merged.tai_start_score == 0.0) & (merged.tai_stop_score == 0.0)
    merged.loc[
        mask,
        TAI_MASKING_NAN_COLS,
    ] = np.nan

    mask = merged["blast_evalue"] == 1
    merged.loc[
        mask,
        BLAST_MASKING_NAN_COLS,
    ] = np.nan

    merged["orf_type"] = merged["orf_type"].map(ORF_TYPE_MAPPING)

    return merged


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
    parser = argparse.ArgumentParser(
        description="Run Random Forest classifier to predict open-reading-frames"
    )
    parser.add_argument(
        "-b", "--blast", required=True, help="Path to BLAST results file"
    )
    parser.add_argument(
        "-s", "--samba", required=True, help="Path to RNASamba results file"
    )
    parser.add_argument(
        "-a", "--alignments", required=True, help="Path to aligned .bed file"
    )
    parser.add_argument(
        "-m",
        "--model",
        type=str,
        help="Path to ORF model",
        default=MODEL_PATH,
    )
    parser.add_argument(
        "-T",
        "--threshold",
        type=float,
        default=0.03,
        help="Use a non-default threshold for classification",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="Path to outdir",
        default=".",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        help="Prefix for the output file name",
        default="xgb",
    )
    parser.add_argument(
        "-mm",
        "--min-score-max-predictions",
        type=float,
        default=0.70,
        help="Minimum score to keep a prediction(s), controlled by max_predictions [default: 0.70]",
    )
    parser.add_argument(
        "-mp",
        "--max-predictions",
        type=int,
        default=1,
        help="Maximum number of predictions to keep per query [default: 1]",
    )
    parser.add_argument(
        "-K", "--keep-raw", action="store_true", help="Keep raw predictions"
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"Version: {__version__}",
    )

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = parse()
    run(args)
