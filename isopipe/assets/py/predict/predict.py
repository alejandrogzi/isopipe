#!/usr/bin/env python3

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.1"

import argparse
from os import PathLike
import logging
from joblib import load
from pathlib import Path
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from typing import List, Union


MODEL = "model/model.joblib"
MODEL_PATH = Path(__file__).parent / MODEL
UNIQUE_TAI_TAG = "UT"
FEATURES: List = [
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
BLAST_COLS: List = [
    "blast_id",
    "blast_encoded",
    "blast_pid",
    "blast_e-value",
    "blast_offset",
    "blast_length",
    "blast_percentage_aligned",
    "blast_o_start",
    "blast_o_end",
    "blast_orf_start",
    "blast_orf_end",
    "blast_strand",
    "blast_chr",
]
TAI_COLS = [
    "tai_id",
    "translationAI_orf_start_coord",
    "translationAI_orf_stop_coord",
    "translationAI_orf_start",
    "translationAI_orf_stop",
    "tai_strand",
    "tai_orf_start",
    "tai_orf_end",
    "tai_chr",
]
TOGA_COLS: List = [
    "toga_id",
    "toga_chr",
    "toga_label",
    "toga_pid",
    "toga_blosum",
    "toga_start",
    "toga_end",
    "toga_strand",
    "toga_key",
    "masked",
]
EXTENDED = [
    "blast_id",
    "key",
    "blast_chr",
    "blast_strand",
    "masked",
    "toga_label",
    "tai_orf_start",
    "tai_orf_end",
    "class0_prob",
    "class1_prob",
]

log = logging.getLogger(__name__)


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
    model = load_model(args.model)

    table = predict(args.blast, args.tai, args.toga, model)
    table = table.loc[:, EXTENDED + FEATURES]

    # if args.overrule:
    #     overrule(table)

    _ = map_to_blocks(table, args.alignments, args.outdir, args.prefix)


def map_to_blocks(
    table: pd.DataFrame,
    alignments: Union[str, PathLike],
    outdir: Union[str, PathLike],
    prefix: str,
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
    table.to_csv(
        f"{args.outdir}/{prefix}predictions.tsv", index=False, header=True, sep="\t"
    )

    table = table.copy()
    table["id"] = [id.split("__")[0] for id in table.blast_id]
    table["tag"] = ["#" + id.split("__")[1] for id in table.blast_id]

    bed = pd.read_csv(alignments, sep="\t", header=None)
    merged = (
        bed.assign(id=bed[3].str.split("__").str[0])
        .merge(
            table[["id", "tai_orf_start", "tai_orf_end", "tag"]], on="id", how="left"
        )
        .dropna(subset=["tai_orf_start"])
        .assign(
            tai_orf_start=lambda df: df["tai_orf_start"].astype(int),
            tai_orf_end=lambda df: df["tai_orf_end"].astype(int),
        )
    )

    merged[6] = merged["tai_orf_start"]
    merged[7] = merged["tai_orf_end"]
    merged[3] += merged["tag"]

    merged = merged.drop(columns=["id", "tai_orf_start", "tai_orf_end", "tag"])

    merged.to_csv(
        f"{outdir}/{prefix}predictions.bed", sep="\t", header=False, index=False
    )

    return


def load_model(model: Union[str, PathLike]) -> RandomForestClassifier:
    """
    Loads a pre-trained RandomForestClassifier model from a specified file path.

    Parameters
    ----------
    model : Union[str, PathLike]
        The file path to the pre-trained RandomForestClassifier model.
        This path can be a string or a Path-like object.

    Returns
    -------
    RandomForestClassifier
        The loaded RandomForestClassifier object.

    Raises
    ------
    SystemExit
        Exits the program with an error message if the model fails to load.

    Example
    -------
    >>> from pathlib import Path
    >>> # Assuming 'my_model.joblib' exists and contains a RandomForestClassifier
    >>> loaded_clf = load_model("my_model.joblib")
    >>> type(loaded_clf)
    <class 'sklearn.ensemble.forest.RandomForestClassifier'>
    """
    try:
        with open(model, "rb") as file:
            clf: RandomForestClassifier = load(file)
    except Exception as err:
        log.error(f"ERROR: failed to load model: {err}!")
        exit(1)

    log.info(f"INFO: finished loading model from -> {model}")
    return clf


def predict(
    blast: Union[str, PathLike],
    tai: Union[str, PathLike],
    toga: Union[str, PathLike],
    model: RandomForestClassifier,
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
    table = read(blast, tai, toga)
    query = table.loc[:, FEATURES]

    predictions = model.predict_proba(query)
    log.info("INFO: Finished predictions!")

    table["class0_prob"] = predictions[:, 0]
    table["class1_prob"] = predictions[:, 1]

    return table


def read(
    blast: Union[str, PathLike, Path],
    tai: Union[str, PathLike, Path],
    toga: Union[str, PathLike, Path],
) -> pd.DataFrame:
    """
    Reads and merges data from BLAST, TAI, and TOGA files into a single DataFrame.

    Parameters
    ----------
    blast : Union[str, PathLike, Path]
        Path to the BLAST input file.
    tai : Union[str, PathLike, Path]
        Path to the TAI input file.
    toga : Union[str, PathLike, Path]
        Path to the TOGA input file.

    Returns
    -------
    pd.DataFrame
        A merged pandas DataFrame containing data from all three input files.

    Example
    -------
    >>> # Assuming 'blast.tsv', 'tai.tsv', 'toga.tsv' exist
    >>> # merged_data = read("blast.tsv", "tai.tsv", "toga.tsv")
    >>> # print(merged_data.head())
    """
    blast = read_blast(blast)
    tai = read_tai(tai)
    toga = read_toga(toga)

    table = merge_tables(blast, tai, toga)
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
    blast["key"] = (
        blast["blast_id"].str.split("__").str[0]
        + ":"
        + blast["blast_orf_start"].astype(str)
        + "-"
        + blast["blast_orf_end"].astype(str)
        + "("
        + blast["blast_strand"].astype(str)
        + ")"
    )
    return blast


def read_tai(path: Union[str, PathLike, Path]) -> pd.DataFrame:
    """
    Reads a TAI tab-separated file and processes it into a pandas DataFrame.

    This function reads the file using predefined column names (TAI_COLS)
    and generates a unique 'key' column based on 'tai_id', 'tai_orf_start',
    and 'tai_orf_end'.

    Parameters
    ----------
    path : Union[str, PathLike, Path]
        The file path to the TAI input file.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the TAI data with an additional 'key' column.

    Example
    -------
    >>> # Assuming 'my_tai.tsv' exists with appropriate tab-separated data
    >>> # tai_df = read_tai("my_tai.tsv")
    >>> # print(tai_df.head())
    """
    tai = pd.read_csv(path, sep="\t", header=None, names=TAI_COLS)

    # INFO: needs to be the cannonical ID
    tai["key"] = (
        tai["tai_id"].str.split(".").str[0]
        + ":"
        + tai["tai_orf_start"].astype(str)
        + "-"
        + tai["tai_orf_end"].astype(str)
        + "("
        + tai["tai_strand"].astype(str)
        + ")"
    )
    tai["tai_id"] = [
        id.replace(".p", "__OR") + f"#{UNIQUE_TAI_TAG}" for id in tai.tai_id
    ]  # WARN: introducing unique TAI tag

    return tai


def read_toga(path: Union[str, PathLike, Path]) -> pd.DataFrame:
    """
    Reads a TOGA tab-separated file, processes it, and removes duplicate entries.

    This function reads the file using predefined column names (TOGA_COLS),
    generates a unique 'key' column based on 'toga_chr', 'toga_start', 'toga_end',
    and 'toga_strand'. It then sorts by 'toga_pid' and removes duplicate keys,
    keeping only the first occurrence.

    Parameters
    ----------
    path : Union[str, PathLike, Path]
        The file path to the TOGA input file.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the processed TOGA data with a 'key' column
        and duplicate keys removed.

    Example
    -------
    >>> # Assuming 'my_toga.tsv' exists with appropriate tab-separated data:
    >>> # ENST00000518498.3#TFF3#229 chr17 FI 75.72 67.90 31125482 31129576 - chr17:31125482-31129576(-) false
    >>> # toga_df = read_toga("my_toga.tsv")
    >>> # print(toga_df.head())
    """
    toga = pd.read_csv(path, sep="\t", header=None, names=TOGA_COLS).drop_duplicates(
        subset="toga_key", keep="first"
    )
    toga.sort_values(by="toga_pid", inplace=True, ascending=False)

    return toga


def merge_tables(
    blast: pd.DataFrame, tai: pd.DataFrame, toga: pd.DataFrame
) -> pd.DataFrame:
    """
    Merges the BLAST, TAI, and TOGA DataFrames into a final comprehensive table.

    The merging process involves:
    1. Inner merging TAI and BLAST DataFrames on the 'key' column.
    2. Extracting chromosome information into 'm_chr' from the merged key.
    3. Creating a 'toga_key' based on chromosome and ORF start/end, considering strand.
    4. Left merging with the TOGA DataFrame on 'toga_key' and TOGA's 'key'.
    5. Adding a binary flag 'toga_overlap_bp' indicating if 'toga_pid' is present.

    Parameters
    ----------
    blast : pd.DataFrame
        The DataFrame containing BLAST data.
    tai : pd.DataFrame
        The DataFrame containing TAI data.
    toga : pd.DataFrame
        The DataFrame containing TOGA data.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame that is the result of merging all three input DataFrames,
        with additional derived columns like 'm_chr', 'toga_key', and 'toga_overlap_bp'.

    Example
    -------
    >>> # Assuming blast_df, tai_df, and toga_df are pre-read DataFrames
    >>> # final_table = merge_tables(blast_df, tai_df, toga_df)
    >>> # print(final_table.head())
    """
    merged = tai.merge(blast, on="key", how="outer")  # INFO: preserving ALL predictions

    merged.blast_strand = merged.blast_strand.fillna(merged.tai_strand)
    merged.tai_orf_start = merged.tai_orf_start.fillna(merged.blast_orf_start)
    merged.tai_orf_end = merged.tai_orf_end.fillna(merged.blast_orf_end)
    merged.blast_chr = merged.blast_chr.fillna(merged.tai_chr)
    merged.blast_id = merged.blast_id.fillna(merged.tai_id)

    merged["toga_key"] = [key.split("_")[-1] for key in merged.key]
    table = merged.merge(toga, on="toga_key", how="left")
    table.drop(columns=["toga_key"], inplace=True)

    # INFO: binary flag: 0 if toga_pid NaN else 1
    table["toga_overlap_bp"] = (~table["toga_pid"].isna()).astype(int)

    return table


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
        "-t", "--tai", required=True, help="Path to translationAI results file"
    )
    parser.add_argument(
        "-g", "--toga", required=True, help="Path to TOGA merged results file"
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"Version: {__version__}",
    )
    parser.add_argument(
        "-O",
        "--overrule",
        action="store_const",
        const=True,
        metavar="Flag to activate TOGA overrule on query reads",
    )
    parser.add_argument(
        "-a", "--alignments", required=True, help="Path to aligned .bed file"
    )
    parser.add_argument(
        "-f", "--fasta", required=False, help="Path to alignments .fa file"
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
        default="",
    )

    args = parser.parse_args()

    if args.overrule and (not args.alignments or not args.fasta):
        parser.error(
            "ERROR: Arguments --alignments (-a) and --fasta (-f) are required when --overrule is set."
        )

    return args


if __name__ == "__main__":
    args = parse()
    run(args)
