#!/usr/bin/env python3

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.1"

import pandas as pd
import numpy as np
import ast
import json
from pathlib import Path
import argparse
import logging
import subprocess
import tempfile
from html import escape
from functools import reduce
import os
from typing import Dict, List, Optional, Union, Any

from pretty_consts import (
    # Reference intron columns and configurations
    REF_INTRONS_COLUMNS,
    REFERENCE_INTRON_ALL_COLUMNS,
    REFERENCE_INTRON_ADD_COLUMNS,
    # Intron retention analysis
    DESCRIPTOR_COLUMNS_TO_VECTOR,
    INTRON_METADATA_COLUMNS_FROM_DESCRIPTOR,
    RETENTION_COLUMNS_FROM_DESCRIPTOR,
    RETAINS_RT_COLUMNS_FROM_DESCRIPTOR,
    HAS_RT_COLUMNS_FROM_DESCRIPTOR,
    RETENTION_HTML_COL,
    RETAINS_RT_HTML_COL,
    HAS_RT_HTML_COL,
    INTRON_METADATA_HTML_COL,
    INTRON_READ_STATUS_COL,
    RETENTION_NESTED_HINT_COLUMN,
    RETENTION_COMMON_KEY,
    RETENTION_METADATA_EXCLUDE,
    RETENTION_NESTED_EXCLUDE,
    INTRON_MOD_METADATA_RENAME,
    INTRON_MOD_RETENTION_SUBMOD_RENAME,
    INTRON_MOD_RETAINS_RT_SUBMOD_RENAME,
    INTRON_MOD_HAS_RT_SUBMOD_RENAME,
    HAS_RT_INTRON,
    INTRON_STATUS_COL,
    HAS_RT_BUCKET,
    # Poly-A tail analysis
    POLYA_METADATA_COLUMNS_FROM_DESCRIPTOR,
    POLYA_MOD_METADATA_RENAME,
    POLYA_METADATA_HTML_COL,
    POLYA_METADATA_EXCLUDE,
    POLYA_RESULT_COLS,
    # Truncation analysis
    TRUNCATION_METADATA_COLUMNS_FROM_DESCRIPTOR,
    TRUNCATION_MOD_METADATA_RENAME,
    TRUNCATION_METADATA_HTML_COL,
    TRUNCATION_METADATA_EXCLUDE,
    TRUNCATION_RESULT_COLS,
    # ORF prediction
    ORF_PREDICTOR_COLS,
    ORF_MOD_METADATA_RENAME,
    ORF_METADATA_HTML_COL,
    ORF_METADATA_EXCLUDE,
    ORF_RESULT_COLS,
    # Tag parsing and read metadata
    TAG_MAP,
    TAG_REGEX,
    # Data type mappings and transformations
    DESCRIPTOR_MAPPING_TYPES,
    # Data processing configurations
    DESCRIPTOR_MAPPING_MISSING_VALUES,
    READ_DECIDE_COLUMNS,
    MAP_IGNORE_CATEGORY,
    MAPPING_KEYS_TO_READ_NAME,
    FLOAT_TO_INT_COLS,
    BED12_COLS,
    BB_SCHEMA_COLUMNS,
    SCHEMA_SQL,
    BED_TO_BIG_BED,
)

PathLike = Union[str, os.PathLike]

log = logging.getLogger(__name__)
logging.basicConfig(encoding='utf-8', level=logging.INFO)


def run(args: argparse.Namespace):
    """
    Main pipeline function to process biological sequence data and generate output files.

    Processes BED format reads through descriptor prettification and quality control
    decision making, then outputs results as BED or BigBed files based on configuration.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments containing:
        - outdir: Output directory path
        - reads: Path to BED12 format reads file
        - predictions: Path to ORF predictions file
        - descriptor: Optional path to descriptor file
        - introns: Optional path to introns file
        - flaws: Flaw detection parameters
        - ignore: Parameters to ignore in processing
        - bigbed: Boolean flag for BigBed output format
        - chrom_sizes: Path to chromosome sizes file (required for BigBed)
        - name: Output filename when not using decision factory

    Returns
    -------
    None
        Creates output files in the specified directory

    Example
    -------
    >>> # Command line usage:
    >>> # python script.py --reads data.bed --predictions orfs.tsv --outdir results/
    >>> args = parse_args()
    >>> run(args)
    >>> # Creates files in results/decision/ directory
    """
    log.info("INFO: starting to prettify descriptor and decision factory!")

    outdir = Path(args.outdir) / "decision"
    outdir.mkdir(parents=True, exist_ok=True)
    log.info(f"INFO: output files going to -> {outdir}")

    reads = pd.read_csv(args.reads, sep="\t", header=None, names=BED12_COLS)

    if not args.introns or not args.descriptor:
        schema = get_pretty_descriptor(reads, args.predictions)
    else:
        schema = get_pretty_descriptor(
            reads, args.predictions, args.descriptor, args.introns
        )

    log.info(f"INFO: got schema for -> {args.reads}. Will start deciding...")

    if args.introns and args.descriptor:
        files = decide(schema, args.flaws, args.ignore)
        log.info(f"INFO: Decision made. Found {len(files)} filled categories")

        for name, bed in files.items():
            name = MAPPING_KEYS_TO_READ_NAME[name]
            outfile = outdir / name

            if args.bigbed:
                bed_to_big_bed(bed, name, args.chrom_sizes, outdir)
            else:
                bed.to_csv(outfile, sep="\t", header=None, index=False)
    else:
        outfile = outdir / args.name

        if args.bigbed:
            bed_to_big_bed(schema, outfile, args.chrom_sizes, outdir)
        else:
            schema.to_csv(outfile, sep="\t", header=None, index=False)

    log.info("INFO: Finished deciding!")
    return


def decide(
    schema: pd.DataFrame, flaws: bool, ignore: Optional[List[str]] = None
) -> dict[str, pd.DataFrame]:
    """
    Classify sequence reads into quality control buckets based on pass/fail rules.

    This function performs quality control classification by analyzing various test results
    and grouping reads into different buckets based on their failure patterns. Reads are
    first separated by intron retention status, then further classified by the number
    and type of test failures.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing sequence analysis results with test status columns.
        Must include 'R*read_intron_status' column and columns specified in READ_DECIDE_COLUMNS.

    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary containing classified DataFrames with the following keys:
        - 'pass': Rows that passed all quality tests
        - 'trash': Rows with 2 or more test failures
        - 'R_has_rt': Rows with "HAS_RT_INTRON" status (separated first)
        - <column_name>: Rows with exactly one failure in the specified test column

    Example
    -------
    >>> input_df = pd.DataFrame({
    ...     'id': ['read1', 'read2', 'read3', 'read4'],
    ...     'R*read_intron_status': ['PASS', 'HAS_RT_INTRON', 'PASS', 'PASS'],
    ...     'test1': ['PASS', 'PASS', 'FAIL', 'PASS'],
    ...     'test2': ['PASS', 'PASS', 'PASS', 'FAIL']
    ... })
    >>> buckets = decide(input_df)
    >>> print(list(buckets.keys()))
    ['R_has_rt', 'pass', 'trash', 'test1', 'test2']
    >>> print(len(buckets['pass']))
    1
    """
    buckets = {}

    print(f"INFO: schema has {len(schema)} reads to decide")

    if ignore:
        ignores = [MAP_IGNORE_CATEGORY[category] for category in ignore]
    else:
        ignores = []

    if "HAS_RT_INTRON" in ignores:
        remaining_reads = schema.copy()
    else:
        remaining_reads, df2 = [
            x for _, x in schema.groupby(schema[INTRON_STATUS_COL] == HAS_RT_INTRON)
        ]

        log.info(f"INFO: Found {len(df2)} flawed reads with RT-introns")
        buckets[HAS_RT_BUCKET] = df2

    # INFO: slice to relevant columns
    sliced = remaining_reads[READ_DECIDE_COLUMNS].copy().reset_index(drop=True)

    fail_mask = sliced.drop(columns="id").apply(lambda x: x != "PASS", axis=1)
    fail_counts = fail_mask.sum(axis=1)

    buckets["pass"] = remaining_reads[
        remaining_reads.id.isin(sliced[fail_counts == 0].id)
    ]
    log.info(f"INFO: Found {len(buckets['pass'])} reads without defects!")

    buckets["trash"] = remaining_reads[
        remaining_reads.id.isin(sliced[fail_counts >= flaws].id)
    ]
    log.info(f"INFO: Found {len(buckets['trash'])} reads with more than {flaws} flaws!")

    # INFO: single-fail buckets
    for col in READ_DECIDE_COLUMNS[1:]:  # INFO: skip "id"
        if col in ignores:
            continue
        mask = (fail_counts == 1) & fail_mask[col]
        buckets[col] = remaining_reads[remaining_reads.id.isin(sliced[mask].id)]

        log.info(f"INFO: Found {len(buckets[col])} reads in {col} category!")

    return buckets


def get_pretty_descriptor(
    reads: pd.DataFrame,
    predictions: PathLike,
    global_descriptor: Optional[PathLike] = None,
    ref_introns: Optional[PathLike] = None,
) -> pd.DataFrame:
    """
    Generate a comprehensive descriptor DataFrame by merging multiple data sources.

    This function combines retention data, poly-A data, truncation data, ORF data,
    and tags into a single DataFrame for downstream analysis.

    Parameters
    ----------
    global_descriptor : PathLike
        Path to the main descriptor file containing sequence analysis results
    ref_introns : PathLike
        Path to reference introns file in tab-separated format
    predictions : PathLike
        Path to ORF predictions file

    Returns
    -------
    pd.DataFrame
        Merged DataFrame containing all analysis results with common retention key

    Example
    -------
    >>> descriptor_df = get_pretty_descriptor(
    ...     "data/global_descriptor.tsv",
    ...     "data/reference_introns.tsv",
    ...     "data/orf_predictions.tsv"
    ... )
    >>> print(descriptor_df.columns.tolist())
    ['id', 'retention_data', 'polya_data', 'truncation_data', 'orf_data', 'tags']
    """
    if not global_descriptor or not ref_introns:
        sources = [
            reads,
            get_orf_data(reads, predictions),
            get_tags(reads),
        ]
    else:
        global_descriptor_df = read_descriptor(global_descriptor)

        sources = [
            reads,
            get_retention_data(global_descriptor_df, ref_introns),
            get_polya_data(global_descriptor_df),
            get_truncation_data(global_descriptor_df),
            get_orf_data(global_descriptor_df, predictions),
            get_tags(global_descriptor_df),
        ]

    merged = reduce(
        lambda left, right: pd.merge(left, right, on=RETENTION_COMMON_KEY, how="inner"),
        sources,
    )

    schema = fill_schema(merged)
    return schema.loc[:, ~schema.columns.duplicated(keep="first")]



def bed_to_big_bed(
    bed: pd.DataFrame,
    filename: PathLike,
    chrom_sizes: PathLike,
    outdir: Union[Path, PathLike],
):
    """
    Convert a BED format DataFrame to BigBed format using UCSC bedToBigBed utility.

    Creates a temporary BED file and converts it to BigBed format with a fixed
    bed12+25 schema for genomic feature visualization.

    Parameters
    ----------
    bed : pd.DataFrame
        DataFrame with BED12+25 format columns in correct order
    filename : PathLike
        Output BigBed file path (typically .bb extension)
    chrom_sizes : PathLike
        Path to chromosome sizes file (tab-separated format)
    outdir : Union[Path, PathLike]
        Directory for temporary file creation

    Returns
    -------
    None
        Writes BigBed file to specified filename

    Raises
    ------
    ValueError
        If input DataFrame is empty
    FileNotFoundError
        If chromosome sizes file doesn't exist

    Example
    -------
    >>> bed_df = fill_schema(processed_data)  # Ensure bed12+25 format
    >>> bed_to_big_bed(bed_df, "output.bb", "hg38.chrom.sizes", "/tmp")
    >>> # Creates output.bb file
    """
    if bed.empty:
        raise ValueError("ERROR: input DataFrame is empty")

    if not os.path.exists(chrom_sizes):
        raise FileNotFoundError(
            f"ERROR: chromosome sizes file not found: {chrom_sizes}"
        )

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False, dir=outdir
    ) as temp_bed:
        temp_bed_path = temp_bed.name
        bed.to_csv(temp_bed, sep="\t", header=False, index=False)

        cmd_parts = [
            BED_TO_BIG_BED,
            f"-as={SCHEMA_SQL}",
            "-type=12+24",  # INFO: static
            temp_bed_path,
            str(chrom_sizes),
            str(filename),
        ]

        log.info(f"INFO: will execute -> {' '.join(cmd_parts)}")

        try:
            result = subprocess.run(
                cmd_parts, capture_output=True, text=True, check=True
            )
            log.info(f"Successfully created BigBed file: {filename} -> {result}")
        except subprocess.CalledProcessError as e:
            log.error(f"{e}")
        finally:
            os.unlink(temp_bed_path)


def get_retention_data(
    global_descriptor: pd.DataFrame, ref_introns: PathLike
) -> pd.DataFrame:
    """
    Extract and process intron retention data from descriptor and reference files.

    Parameters
    ----------
    global_descriptor : pd.DataFrame
        Global descriptor DataFrame containing sequence analysis data
    ref_introns : PathLike
        Path to reference introns file

    Returns
    -------
    pd.DataFrame
        Processed DataFrame containing retention analysis results

    Example
    -------
    >>> retention_df = get_retention_data(descriptor_df, "ref_introns.tsv")
    >>> print(retention_df.shape)
    (1000, 6)
    """
    ref_introns_df = read_ref_introns(ref_introns)
    ref_introns_extensible_data = ref_introns_df[REFERENCE_INTRON_ALL_COLUMNS].to_dict(
        orient="index"
    )

    merged = get_retention_submodules(global_descriptor, ref_introns_extensible_data)
    return merged


def fill_schema(merged: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure DataFrame conforms to BED12 and BB schema requirements.

    This function standardizes the DataFrame by adding any missing required columns
    (filled with NaN values) and reordering columns to match the expected schema.
    This is typically used as a final step to prepare data for downstream analysis
    or export in standardized formats.

    Parameters
    ----------
    merged : pd.DataFrame
        Input DataFrame that may be missing some schema columns

    Returns
    -------
    pd.DataFrame
        DataFrame with all required BED12 and BB schema columns in correct order,
        missing columns filled with NaN values

    Example
    -------
    >>> incomplete_df = pd.DataFrame({'chrom': ['chr1'], 'start': [1000]})
    >>> complete_df = fill_schema(incomplete_df)
    >>> print(len(complete_df.columns) == len(BED12_COLS + BB_SCHEMA_COLUMNS))
    True
    >>> print(complete_df['end'].isna().iloc[0])
    True
    """
    required_cols = BED12_COLS + BB_SCHEMA_COLUMNS

    # INFO: add missing columns with NaN
    for col in required_cols:
        if col not in merged.columns:
            merged[col] = np.nan

    # INFO: enforcing schema order
    return merged[required_cols]


def get_metadata(
    global_descriptor: pd.DataFrame,
    columns_from_descriptor: List[str],
    rename_map: Dict[str, str],
    html_col_name: str,
    exclude_fields: List[str],
) -> pd.DataFrame:
    """
    Extract metadata from descriptor DataFrame and generate HTML representation.

    This is a generic function for extracting specific columns, renaming them,
    and creating an HTML column for visualization purposes.

    Parameters
    ----------
    global_descriptor : pd.DataFrame
        Full descriptor DataFrame containing all sequence data
    columns_from_descriptor : List[str]
        List of column names to extract from the descriptor
    rename_map : Dict[str, str]
        Mapping of old column names to new column names
    html_col_name : str
        Name of the column that will contain HTML representation
    exclude_fields : List[str]
        List of field names to exclude from HTML rendering

    Returns
    -------
    pd.DataFrame
        Processed metadata DataFrame with HTML column added

    Example
    -------
    >>> metadata_df = get_metadata(
    ...     descriptor_df,
    ...     ['col1', 'col2', 'col3'],
    ...     {'col1': 'new_col1', 'col2': 'new_col2'},
    ...     'html_representation',
    ...     ['id']
    ... )
    >>> print(metadata_df.columns.tolist())
    ['new_col1', 'new_col2', 'col3', 'html_representation']
    """
    metadata = global_descriptor[columns_from_descriptor].rename(columns=rename_map)

    metadata[html_col_name] = [
        row_to_html(dict(zip(metadata.columns, row)), exclude_fields)
        for row in metadata.itertuples(index=False, name=None)
    ]

    return metadata


def get_retention_metadata(global_descriptor: pd.DataFrame) -> pd.DataFrame:
    """
    Extract intron retention metadata with predefined column mappings.

    Parameters
    ----------
    global_descriptor : pd.DataFrame
        Global descriptor DataFrame

    Returns
    -------
    pd.DataFrame
        Intron retention metadata with HTML representation

    Example
    -------
    >>> retention_meta = get_retention_metadata(descriptor_df)
    >>> print('intron_metadata_html' in retention_meta.columns)
    True
    """
    metadata = get_metadata(
        global_descriptor,
        INTRON_METADATA_COLUMNS_FROM_DESCRIPTOR,
        INTRON_MOD_METADATA_RENAME,
        INTRON_METADATA_HTML_COL,
        RETENTION_METADATA_EXCLUDE,
    )

    metadata[INTRON_READ_STATUS_COL] = [
        "FLAWED"
        if metadata.iloc[idx].R_read_exon_status == "RETAINS_SPLICED_INTRON"
        or metadata.iloc[idx].R_read_intron_status == "HAS_RT_INTRON"
        else "PASS"
        for idx in metadata.index
    ]

    return metadata


def get_polya_metadata(global_descriptor: pd.DataFrame) -> pd.DataFrame:
    """
    Extract poly-A site metadata with predefined column mappings.

    Parameters
    ----------
    global_descriptor : pd.DataFrame
        Global descriptor DataFrame

    Returns
    -------
    pd.DataFrame
        Poly-A metadata with HTML representation

    Example
    -------
    >>> polya_meta = get_polya_metadata(descriptor_df)
    >>> print('polya_metadata_html' in polya_meta.columns)
    True
    """
    return get_metadata(
        global_descriptor,
        POLYA_METADATA_COLUMNS_FROM_DESCRIPTOR,
        POLYA_MOD_METADATA_RENAME,
        POLYA_METADATA_HTML_COL,
        POLYA_METADATA_EXCLUDE,
    )


def get_retention_submodules(
    global_descriptor: pd.DataFrame, extensible_data: Optional[Dict[str, Any]] = None
) -> pd.DataFrame:
    """
    Process retention submodules and determine read status.

    Combines retention metadata with various retention submodules and assigns
    a read status (FLAWED or PASS) based on exon and intron status.

    Parameters
    ----------
    global_descriptor : pd.DataFrame
        Global descriptor DataFrame
    extensible_data : Optional[Dict[str, Any]], default None
        Additional extensible data for enrichment

    Returns
    -------
    pd.DataFrame
        Combined retention data with read status classification

    Example
    -------
    >>> retention_sub = get_retention_submodules(descriptor_df)
    >>> print(retention_sub['read_status'].value_counts())
    PASS      800
    FLAWED    200
    """

    dfs = [
        get_retention_metadata(global_descriptor),
        parse_retention_submodule(
            global_descriptor,
            RETENTION_COLUMNS_FROM_DESCRIPTOR,
            INTRON_MOD_RETENTION_SUBMOD_RENAME,
            RETENTION_HTML_COL,
            extensible_data,
        ),
        parse_retention_submodule(
            global_descriptor,
            RETAINS_RT_COLUMNS_FROM_DESCRIPTOR,
            INTRON_MOD_RETAINS_RT_SUBMOD_RENAME,
            RETAINS_RT_HTML_COL,
            extensible_data,
        ),
        parse_retention_submodule(
            global_descriptor,
            HAS_RT_COLUMNS_FROM_DESCRIPTOR,
            INTRON_MOD_HAS_RT_SUBMOD_RENAME,
            HAS_RT_HTML_COL,
            extensible_data,
        ),
    ]

    merged = reduce(
        lambda left, right: pd.merge(left, right, on=RETENTION_COMMON_KEY, how="inner"),
        dfs,
    )

    return merged[
        [
            RETENTION_COMMON_KEY,
            INTRON_READ_STATUS_COL,
            INTRON_METADATA_HTML_COL,
            RETENTION_HTML_COL,
            RETAINS_RT_HTML_COL,
            HAS_RT_HTML_COL,
            INTRON_STATUS_COL,
        ]
    ]


def read_ref_introns(reference_introns: PathLike) -> pd.DataFrame:
    """
    Read and process reference introns file.

    Reads a tab-separated file containing reference intron coordinates,
    creates a unique key from chromosome and coordinates, and removes duplicates.

    Parameters
    ----------
    reference_introns : PathLike
        Path to reference introns file (tab-separated, no header)

    Returns
    -------
    pd.DataFrame
        Processed DataFrame with intron coordinates indexed by genomic key

    Example
    -------
    >>> ref_introns = read_ref_introns("reference_introns.tsv")
    >>> print(ref_introns.index[0])
    'chr1:1000-2000'
    """
    df = pd.read_csv(
        reference_introns, sep="\t", header=None, names=REF_INTRONS_COLUMNS
    )

    df["key"] = (
        df["chrom"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
    )
    df.drop_duplicates(subset="key", inplace=True)
    df.set_index("key", inplace=True)

    return df


def read_descriptor(path: PathLike) -> pd.DataFrame:
    """
    Read and process the main descriptor file.

    Reads a tab-separated descriptor file and applies various transformations:
    - Converts string representations of lists to actual Python lists
    - Applies column-specific value mappings
    - Handles missing values
    - Converts float columns to integers where appropriate

    Parameters
    ----------
    path : PathLike
        Path to the descriptor file

    Returns
    -------
    pd.DataFrame
        Processed descriptor DataFrame ready for analysis

    Example
    -------
    >>> descriptor = read_descriptor("global_descriptor.tsv")
    >>> print(descriptor.shape)
    (5000, 25)
    """
    descriptor = pd.read_csv(path, sep="\t")

    # Step 1: Convert descriptor columns into proper lists
    for col in DESCRIPTOR_COLUMNS_TO_VECTOR:
        descriptor[col] = descriptor[col].apply(convert_to_proper_list)

    # Step 2: Apply column-specific mappings
    _apply_column_mappings(descriptor)

    # Step 3: Handle missing values
    descriptor.fillna(DESCRIPTOR_MAPPING_MISSING_VALUES, inplace=True)

    # Step 4: Cast float columns to int and round
    descriptor[FLOAT_TO_INT_COLS] = descriptor[FLOAT_TO_INT_COLS].astype(int)
    descriptor = descriptor.round(3)

    return descriptor


def _apply_column_mappings(df: pd.DataFrame) -> None:
    """
    Apply value mappings to specific columns in the descriptor DataFrame.

    This is a private helper function that applies predefined mappings
    to both scalar values and lists within DataFrame columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to modify in place

    Example
    -------
    >>> # This is called internally by read_descriptor()
    >>> _apply_column_mappings(descriptor_df)
    """

    # Generic helper to apply mapping to scalars or lists
    def _map_values(x: Any, mapping: Dict[Any, Any]) -> Any:
        if isinstance(x, list):
            return [mapping.get(item, item) for item in x]
        elif isinstance(x, (str, bool)):
            return mapping.get(x, x)
        return x

    for col, mapping in DESCRIPTOR_MAPPING_TYPES.items():
        df[col] = df[col].apply(lambda x, m=mapping: _map_values(x, m))


def parse_retention_submodule(
    global_descriptor: pd.DataFrame,
    columns: List[str],
    rename_schema: Dict[str, str],
    html_column_name: str,
    extensible_data: Optional[Dict[str, Any]] = None,
) -> pd.DataFrame:
    """
    Parse a specific retention submodule and generate HTML representation.

    Parameters
    ----------
    global_descriptor : pd.DataFrame
        Global descriptor DataFrame
    columns : List[str]
        List of columns to extract for this submodule
    rename_schema : Dict[str, str]
        Column renaming schema
    html_column_name : str
        Name for the HTML output column
    extensible_data : Optional[Dict[str, Any]], default None
        Additional data for extensible columns

    Returns
    -------
    pd.DataFrame
        Processed submodule data with HTML representation

    Example
    -------
    >>> submodule_df = parse_retention_submodule(
    ...     descriptor_df,
    ...     ['col1', 'col2'],
    ...     {'col1': 'new_col1'},
    ...     'html_output'
    ... )
    """
    df = global_descriptor[columns].rename(columns=rename_schema)
    return df_to_nested_html_rows(df, html_column_name, extensible_data)


def df_to_nested_html_rows(
    df: pd.DataFrame,
    output_col_name: str,
    extensible_data: Optional[Dict[str, Any]] = None,
) -> pd.DataFrame:
    """
    Convert DataFrame rows to nested HTML table representations.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame to process
    output_col_name : str
        Name of the column to store HTML output
    extensible_data : Optional[Dict[str, Any]], default None
        Additional data for extensible columns

    Returns
    -------
    pd.DataFrame
        DataFrame with added HTML column containing nested table representations

    Example
    -------
    >>> html_df = df_to_nested_html_rows(data_df, 'html_table')
    >>> print(html_df['html_table'].iloc[0][:50])
    '<table border="1" style="border-collapse:collapse;">'
    """
    df[output_col_name] = [
        nested_row_to_html(
            dict(zip(df.columns, row)),
            RETENTION_NESTED_HINT_COLUMN,
            RETENTION_NESTED_EXCLUDE,
            extensible_data,
            REFERENCE_INTRON_ADD_COLUMNS,
        )
        for row in df.itertuples(index=False, name=None)
    ]

    return df


def convert_to_proper_list(value: Any) -> List[Any]:
    """
    Convert various input types to Python lists.

    Handles string representations of lists (JSON format or Python literal format),
    existing lists, and scalar values. Supports JavaScript-style boolean conversion.

    Parameters
    ----------
    value : Any
        Input value to convert to list

    Returns
    -------
    List[Any]
        Converted list, empty list if conversion fails

    Example
    -------
    >>> convert_to_proper_list('[1, 2, 3]')
    [1, 2, 3]
    >>> convert_to_proper_list('[true, false]')
    [True, False]
    >>> convert_to_proper_list('single_value')
    ['single_value']
    """
    if pd.isna(value):
        return []
    elif isinstance(value, list):
        return value
    elif isinstance(value, str):
        try:
            if value.strip() == "[]":
                return []

            # First try json.loads (handles true/false natively)
            if value.startswith("[") and value.endswith("]"):
                try:
                    return json.loads(value)
                except json.JSONDecodeError:
                    # If json fails, try ast.literal_eval with JavaScript conversion
                    json_compatible = value.replace("true", "True").replace(
                        "false", "False"
                    )
                    return ast.literal_eval(json_compatible)
            else:
                return [value]
        except (ValueError, SyntaxError, json.JSONDecodeError):
            return []
    else:
        return []


def row_to_html(
    row: Dict[str, Any],
    exclude: Optional[List[str]] = None,
    extensible_lookup: Optional[Dict[str, Any]] = None,
    extensible_cols: Optional[List[str]] = None,
) -> str:
    """
    Convert a dictionary row to an HTML table representation.

    Parameters
    ----------
    row : Dict[str, Any]
        Dictionary containing row data
    exclude : Optional[List[str]], default None
        List of keys to exclude from the table
    extensible_lookup : Optional[Dict[str, Any]], default None
        Lookup dictionary for extensible data
    extensible_cols : Optional[List[str]], default None
        List of extensible column names

    Returns
    -------
    str
        HTML table string representation of the row

    Example
    -------
    >>> html_table = row_to_html(
    ...     {'name': 'Gene1', 'score': 0.95, 'status': 'PASS'},
    ...     exclude=['id']
    ... )
    >>> print('<table' in html_table)
    True
    """
    if exclude is None:
        exclude = []

    # Columns to include (skip 'id')
    cols = [c for c in row.keys() if c not in exclude]

    # Build header
    header = "".join(f"<th>{escape(col)}</th>" for col in cols)

    # Build body
    body = "".join(f"<td>{escape(str(row.get(col)))}</td>" for col in cols)

    if extensible_cols and extensible_lookup is not None:
        for col in extensible_cols:
            val = extensible_lookup.get(col)
            header += f"<th>{escape(col)}</th>"
            body += f"<td>{escape(str(val))}</td>"

    header_html = "<tr>" + header + "</tr>"
    body_html = "<tr>" + body + "</tr>"

    # Assemble table
    table_html = f"""<table border="1" style="border-collapse:collapse;"><thead>{header_html}</thead><tbody>{body_html}</tbody></table>"""
    return table_html.strip()


def nested_row_to_html(
    row: Dict[str, Any],
    length_hint_col: str,
    exclude: Optional[List[str]] = None,
    extensible_lookup: Optional[Dict[str, Any]] = None,
    extensible_cols: Optional[List[str]] = None,
) -> str:
    """
    Convert a row containing list data into a nested HTML table.

    Creates an HTML table where each row represents an "event" from the list data,
    with columns showing the corresponding values for each event.

    Parameters
    ----------
    row : Dict[str, Any]
        Dictionary containing row data with list-like columns
    length_hint_col : str
        Column name to use for determining the number of events
    exclude : Optional[List[str]], default None
        List of keys to exclude from the table
    extensible_lookup : Optional[Dict[str, Any]], default None
        Lookup dictionary for additional data
    extensible_cols : Optional[List[str]], default None
        List of extensible column names

    Returns
    -------
    str
        HTML table string with nested event structure

    Example
    -------
    >>> nested_html = nested_row_to_html(
    ...     {'scores': [0.8, 0.9, 0.7], 'positions': [100, 200, 300]},
    ...     'scores'
    ... )
    >>> print('EVENT 1' in nested_html)
    True
    """
    if exclude is None:
        exclude = []

    # Columns to include (skip 'id')
    cols = [c for c in row.keys() if c not in exclude]
    cols.insert(0, "EVENT")

    # Assume all list-like columns have equal length
    vals = row.get(length_hint_col, [])
    if isinstance(vals, pd.Series):  # unwrap if it's a Series of lists
        vals = vals.iloc[0]
    if not isinstance(vals, (list, tuple)):
        vals = []

    events = len(vals)

    if events < 1:
        return "<table border='1' style='border-collapse:collapse;'></table>"

    # Build header
    header = "".join(f"<th>{escape(col)}</th>" for col in cols)
    if extensible_cols:
        for col in extensible_cols:
            header += f"<th>{escape(col)}</th>"
    header_html = f"<tr>{header}</tr>"

    # Build body
    body_rows = []
    for event in range(events):
        cells = []
        for col in cols:
            if col == "EVENT":
                cells.append(f"<td>{escape(f'EVENT {event+1}')}</td>")
            else:
                cells.append(f"<td>{escape(str(row[col][event]))}</td>")

        if extensible_cols and extensible_lookup is not None:
            key = vals[event]
            for col in extensible_cols:
                val = extensible_lookup.get(key, {}).get(col, "")
                cells.append(f"<td>{escape(str(val))}</td>")

        body_rows.append("<tr>" + "".join(cells) + "</tr>")

    body_html = "".join(body_rows)

    # Assemble table
    return f"""<table border="1" style="border-collapse:collapse;">
<thead>{header_html}</thead><tbody>{body_html}</tbody></table>""".strip()


def get_polya_data(global_descriptor: pd.DataFrame) -> pd.DataFrame:
    """
    Extract poly-A site analysis data.

    Parameters
    ----------
    global_descriptor : pd.DataFrame
        Global descriptor DataFrame

    Returns
    -------
    pd.DataFrame
        Poly-A data subset with relevant columns

    Example
    -------
    >>> polya_data = get_polya_data(descriptor_df)
    >>> print(polya_data.shape[1] <= len(POLYA_RESULT_COLS))
    True
    """
    metadata = get_polya_metadata(global_descriptor)
    return metadata[POLYA_RESULT_COLS]


def get_truncation_metadata(global_descriptor: pd.DataFrame) -> pd.DataFrame:
    """
    Extract truncation analysis metadata with predefined mappings.

    Parameters
    ----------
    global_descriptor : pd.DataFrame
        Global descriptor DataFrame

    Returns
    -------
    pd.DataFrame
        Truncation metadata with HTML representation

    Example
    -------
    >>> trunc_meta = get_truncation_metadata(descriptor_df)
    >>> print('truncation_metadata_html' in trunc_meta.columns)
    True
    """
    return get_metadata(
        global_descriptor,
        TRUNCATION_METADATA_COLUMNS_FROM_DESCRIPTOR,
        TRUNCATION_MOD_METADATA_RENAME,
        TRUNCATION_METADATA_HTML_COL,
        TRUNCATION_METADATA_EXCLUDE,
    )


def get_truncation_data(global_descriptor: pd.DataFrame) -> pd.DataFrame:
    """
    Extract truncation analysis data.

    Parameters
    ----------
    global_descriptor : pd.DataFrame
        Global descriptor DataFrame

    Returns
    -------
    pd.DataFrame
        Truncation data subset with relevant columns

    Example
    -------
    >>> trunc_data = get_truncation_data(descriptor_df)
    >>> print(trunc_data.shape[1] <= len(TRUNCATION_RESULT_COLS))
    True
    """
    metadata = get_truncation_metadata(global_descriptor)
    return metadata[TRUNCATION_RESULT_COLS]


def get_tags(global_descriptor: pd.DataFrame) -> pd.DataFrame:
    """
    Extract and parse tags from sequence identifiers.

    Parses tags embedded in sequence names (after "__") and converts them
    to structured data columns.

    Parameters
    ----------
    global_descriptor : pd.DataFrame
        Global descriptor DataFrame with 'id' column

    Returns
    -------
    pd.DataFrame
        DataFrame with parsed tag columns

    Example
    -------
    >>> tags_df = get_tags(descriptor_df)
    >>> print('id' in tags_df.columns)
    True
    """
    return pd.DataFrame(
        [get_tags_from_read(r) for r in global_descriptor.id],
        index=global_descriptor.id,
    ).reset_index()


def get_tags_from_read(name: str) -> Dict[str, Union[int, bool]]:
    """
    Parse tags from a single sequence identifier.

    Extracts structured tag information from sequence names formatted as:
    "sequence_name__tag1#tag2:value#tag3" etc.

    Parameters
    ----------
    name : str
        Sequence identifier containing embedded tags

    Returns
    -------
    Dict[str, Union[int, bool]]
        Dictionary of parsed tags with their values

    Example
    -------
    >>> tags = get_tags_from_read("seq1__quality:95#filtered#length:1000")
    >>> print(tags.get('quality'))
    95
    """
    tags = {}
    _, _, tag_str = name.partition("__")  # everything after "__"
    for match in tag_str.split("#"):
        m = TAG_REGEX.fullmatch(match)
        if not m:
            continue
        key, val = m.groups()
        col = TAG_MAP.get(key)
        if col:
            if val is not None:
                tags[col] = int(val)
            else:
                tags[col] = True
    return tags


def get_orf_metadata(path: PathLike) -> pd.DataFrame:
    """
    Read and process ORF (Open Reading Frame) prediction metadata.

    Reads ORF prediction results, applies column renaming, rounds numeric values,
    and generates HTML representation for visualization.

    Parameters
    ----------
    path : PathLike
        Path to ORF predictions file (tab-separated)

    Returns
    -------
    pd.DataFrame
        Processed ORF metadata with HTML representation

    Example
    -------
    >>> orf_meta = get_orf_metadata("orf_predictions.tsv")
    >>> print('O_metadata_html' in orf_meta.columns)
    True
    """
    data = (
        pd.read_csv(path, sep="\t")
        .rename(columns=ORF_MOD_METADATA_RENAME)[ORF_PREDICTOR_COLS]
        .round({"O_blast_pid": 3, "O_class_0_prob": 3, "O_blast_percentage_aligned": 3})
    )

    data[ORF_METADATA_HTML_COL] = [
        row_to_html(dict(zip(data.columns, row)), ORF_METADATA_EXCLUDE)
        for row in data.itertuples(index=False, name=None)
    ]

    return data[ORF_RESULT_COLS]


def get_orf_data(
    global_descriptor: pd.DataFrame, predictions: PathLike
) -> pd.DataFrame:
    """
    Extract and merge ORF prediction data with sequence identifiers.

    Matches ORF predictions to sequence identifiers by parsing the sequence names
    and merging on the prefix portion before "__".

    Parameters
    ----------
    global_descriptor : pd.DataFrame
        Global descriptor DataFrame with sequence identifiers | Optionally bed file [4th col -> id name]
    predictions : PathLike
        Path to ORF predictions file

    Returns
    -------
    pd.DataFrame
        Merged DataFrame containing ORF scores and metadata for matched sequences

    Example
    -------
    >>> orf_data = get_orf_data(descriptor_df, "orf_predictions.tsv")
    >>> print(orf_data.columns.tolist())
    ['id', 'O_read_orf_score', 'O_metadata_html']
    """
    ids_df = global_descriptor.id.to_frame()
    ids_df[["prefix", "suffix_main"]] = ids_df["id"].str.split("__", n=1, expand=True)

    orf = get_orf_metadata(predictions)
    orf[["prefix", "suffix"]] = orf["O_blast_id"].str.split("__", n=1, expand=True)

    merged = pd.merge(
        ids_df, orf, on="prefix", how="left"
    )  # INFO: orf contains also NMDs
    merged = merged[
        merged.apply(lambda r: r["suffix_main"].endswith(r["suffix"]), axis=1)
    ]

    return merged[["id", "O_read_orf_score", "O_metadata_html"]]


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
        description="Isopipe descriptor prettifier that produces bigBed inputs"
    )
    parser.add_argument(
        "-r", "--reads", required=True, help="Path to aligned .bed file"
    )
    parser.add_argument(
        "-d", "--descriptor", required=False, help="Path to global_descriptor.tsv file"
    )
    parser.add_argument(
        "-p", "--predictions", required=True, help="Path to ORF predictions.tsv file"
    )
    parser.add_argument(
        "-i",
        "--introns",
        required=False,
        help="Path to reference_introns.tsv, output of iso-classify",
    )
    parser.add_argument(
        "-bb",
        "--bigbed",
        action="store_const",
        const=True,
        metavar="Flag to activate bigBed conversion",
    )
    parser.add_argument(
        "-cs",
        "--chrom-sizes",
        type=str,
        required=False,
        metavar="Path to chrom.sizes file",
    )
    parser.add_argument(
        "-f",
        "--flaws",
        type=int,
        help="Minimum number of flaws to be send to trash",
        default=2,
    )
    parser.add_argument(
        "-ig",
        "--ignore",
        choices=["retention", "truncation", "intraprimming", "rtintron"],
        nargs="+",
        help="Categories to exclude from decision. Options: retention, truncation, intraprimming",
        default=None,
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="Path to outdir",
        default=".",
    )
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        required=False,
        help="Optinal user-desired output file name [Only for non-cannonical decisions]",
        default="output",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"Version: {__version__}",
    )

    args = parser.parse_args()

    if not args.descriptor and not args.introns:
        log.info(
            f"INFO: running in non-cannonical mode for: {args.reads}. Will fill with NaNs"
        )

    if args.bigbed:
        if not args.chrom_sizes:
            raise FileNotFoundError(
                "ERROR: chrom.sizes file is necessary to run bigBed converter!"
            )
        log.info(
            f"INFO: bigBed conversion activated, will prepare files and schema for {args.reads} using {args.chrmom_sizes}"
        )

    return args


if __name__ == "__main__":
    args = parse()
    run(args)
