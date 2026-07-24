#!/usr/bin/env python3
"""Rename ORF predictions in a BED12 file and its companion TSV records.

This script takes a BED12 file of predicted open reading frames (ORFs) and a
tab-separated (TSV) file holding additional metadata for each prediction, and
rewrites the shared ``id`` column of both files with human-readable,
information-rich names.

Input identifiers must have the shape ``ROOT_ORF.N[@M][#DU]`` or
``ROOT.pN[@M][#DU]``. The terminal marker separates the opaque root from ORF
metadata; dots and underscores inside the root are preserved. Metadata is
converted into explicit tags in the renamed identifier. The result has the
general shape::

    [PREFIX#...]ROOT[#TAG]...[#SC{SCORE}][#DU][#OR{N}][#TI{N}][#IN{M}][@PROTEIN]

With ``--rebase``, the root is replaced by a short, collision-free hash while
the ORF metadata remains visible as tags. Predictions sharing a root therefore
share a hash and remain distinguishable by their metadata.

Reserved characters ``#`` and ``@`` are not permitted inside user-supplied
prefixes or tags because they are used as the structural separators of the
generated names.

Example
-------
    rename_predictions.py \\
        --bed predictions.bed \\
        --tsv predictions.tsv \\
        --unmask --tag LIVER --append-orf-score \\
        --output-bed renamed.bed --output-tsv renamed.tsv

See :func:`parse_args` for the full list of command-line options.
"""

import argparse
import hashlib
import re
import sys
from itertools import repeat
from typing import Iterable, NamedTuple, Optional

import pandas as pd

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.4"

BED_COLS = [
    "chrom",
    "chrom_start",
    "chrom_end",
    "id",
    "score",
    "strand",
    "thick_start",
    "thick_end",
    "item_rgb",
    "block_count",
    "block_sizes",
    "block_starts",
]
RESERVED = {"#", "@"}
ORF_DELIMITER = "_ORF"
ORF_METADATA_PATTERN = re.compile(
    r"^\.(?P<orf_number>\d+)(?:@(?P<inner_number>\d+))?(?P<duplicate>#DU)?$"
)
TRANSLATION_ID_PATTERN = re.compile(
    r"^(?P<root>.+)\.p(?P<orf_number>\d+)"
    r"(?:@(?P<inner_number>\d+))?(?P<duplicate>#DU)?$"
)


class PredictionId(NamedTuple):
    """Parsed components of an input prediction identifier."""

    root: str
    orf_number: str
    translation_index: Optional[str]
    inner_number: Optional[str]
    is_duplicate: bool


def parse_prediction_id(value: str) -> PredictionId:
    """Parse either supported prediction ID without interpreting its root."""
    if ORF_DELIMITER in value:
        parts = value.split(ORF_DELIMITER)
        if len(parts) != 2 or not parts[0]:
            raise ValueError(
                "Prediction ID must contain exactly one '_ORF' delimiter after "
                f"a non-empty root: {value!r}"
            )

        root, metadata = parts
        match = ORF_METADATA_PATTERN.fullmatch(metadata)
        translation_index = None
    else:
        match = TRANSLATION_ID_PATTERN.fullmatch(value)
        root = match.group("root") if match is not None else ""
        translation_index = match.group("orf_number") if match is not None else None

    if match is None:
        raise ValueError(
            "Prediction ID must end in '_ORF.N' or '.pN', optionally followed "
            f"by '@M' and '#DU': {value!r}"
        )

    return PredictionId(
        root=root,
        orf_number=match.group("orf_number"),
        translation_index=translation_index,
        inner_number=match.group("inner_number"),
        is_duplicate=match.group("duplicate") is not None,
    )


def validate_token(value: str, label: str) -> str:
    """Reject tokens containing reserved separator characters.

    Prefixes and tags must not include ``#`` or ``@`` because those characters
    are used to separate the components of a generated name. ``label`` is used
    to produce a helpful error message identifying which argument was invalid.
    """
    if any(char in value for char in RESERVED):
        raise ValueError(f"{label} cannot contain '#' or '@': {value!r}")
    return value


def short_hashes(values: Iterable[str]) -> dict[str, str]:
    """Map each unique root to the shortest collision-free 64-bit hash prefix.

    A 64-bit BLAKE2s digest is computed once per distinct root. Sorted integer
    digests make it possible to derive the required prefix width from adjacent
    values in one pass, avoiding repeated large temporary sets. The digest map
    is then converted to its final string values in place to limit peak memory.

    Raises
    ------
    ValueError
        If two distinct roots have the same complete 64-bit digest.
    """
    digests = {}
    for value in values:
        if value not in digests:
            digest = hashlib.blake2s(value.encode(), digest_size=8).digest()
            digests[value] = int.from_bytes(digest, byteorder="big")

    sorted_digests = sorted(digests.values())
    width = 4
    for previous, current in zip(sorted_digests, sorted_digests[1:]):
        differing_bits = previous ^ current
        if differing_bits == 0:
            raise ValueError(
                "Two distinct prediction roots produced the same 64-bit hash"
            )
        shared_nibbles = (64 - differing_bits.bit_length()) // 4
        width = max(width, shared_nibbles + 1)

    del sorted_digests

    for root, digest in digests.items():
        digests[root] = f"h{digest:016x}"[: width + 1]
    return digests


def format_score(value: str) -> str:
    """Format an ORF score as a compact integer string.

    The value is parsed as a integer (0-100) and rendered with two significant digits
    (``%.2g``) so scores from the TSV are normalized consistently in the
    generated names.
    """
    try:
        return f"{float(value) * 100:.0f}"
    except ValueError as error:
        raise ValueError(f"Invalid ORF score: {value!r}") from error


def parse_args() -> argparse.Namespace:
    """Define and parse the command-line interface.

    Returns
    -------
    argparse.Namespace
        Parsed arguments. Available attributes:

        ``bed``
            Path to the input BED12 file.
        ``tsv``
            Path to the input TSV metadata file.
        ``output_bed``
            Destination path for the renamed BED12 file.
        ``output_tsv``
            Destination path for the renamed TSV file.
        ``deactivate``
            When set, copy the inputs unchanged to the outputs and exit.
        ``rebase``
            Replace the parsed root with a short content hash while preserving
            ORF metadata as tags.
        ``append_orf_score``
            Append the ORF score (from ``score_column``) to each name.
        ``score_column``
            TSV column used as the ORF score source.
        ``append_protein_name``
            Append the BLAST reference id as the protein component.
        ``custom_prefix``
            Extra prefix prepended to every generated name.
        ``unmask``
            Prepend the ``UNMASK`` prefix to every generated name.
        ``tag``
            Zero or more tags appended to every generated name.
    """
    parser = argparse.ArgumentParser(
        description="Rename BED predictions and their TSV records."
    )
    parser.add_argument("--bed", required=True)
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--output-bed", default="renamed.bed")
    parser.add_argument("--output-tsv", default="renamed.tsv")
    parser.add_argument("--deactivate", action="store_true")
    parser.add_argument("--rebase", action="store_true")
    parser.add_argument("--append-orf-score", action="store_true")
    parser.add_argument("--score-column", default="prob_coding")
    parser.add_argument("--append-protein-name", action="store_true")
    parser.add_argument("--custom-prefix")
    parser.add_argument("--unmask", action="store_true")
    parser.add_argument("--tag", action="append", default=[])
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    bed = pd.read_csv(args.bed, sep="\t", header=None, dtype=str, keep_default_na=False)
    if bed.shape[1] < len(BED_COLS):
        raise ValueError(f"Expected BED12, found {bed.shape[1]} columns")

    if bed.shape[1] > len(BED_COLS):
        print(f"WARN: BED12 has {bed.shape[1]} columns; dropping extra columns")
        bed.drop(columns=bed.columns[len(BED_COLS) :], inplace=True)

    bed.columns = BED_COLS

    tsv = pd.read_csv(args.tsv, sep="\t", dtype=str, keep_default_na=False)
    if "id" not in tsv.columns:
        raise ValueError("TSV is missing the 'id' column")
    if tsv["id"].duplicated().any():
        raise ValueError("TSV contains duplicate IDs")

    bed_ids = set(bed["id"])
    tsv_ids = set(tsv["id"])
    if bed_ids != tsv_ids:
        missing_tsv = sorted(bed_ids - tsv_ids)[:5]
        missing_bed = sorted(tsv_ids - bed_ids)[:5]
        raise ValueError(
            f"BED/TSV IDs differ; missing from TSV: {missing_tsv}; "
            f"missing from BED: {missing_bed}"
        )

    if args.deactivate:
        bed.to_csv(args.output_bed, sep="\t", header=False, index=False)
        tsv.to_csv(args.output_tsv, sep="\t", index=False)
        return

    prefix_parts = []
    if args.unmask:
        prefix_parts.append("UNMASK")
    if args.custom_prefix:
        prefix_parts.append(validate_token(args.custom_prefix, "custom prefix"))

    tags = [validate_token(tag, "tag").upper() for tag in args.tag]

    if args.append_orf_score and args.score_column not in tsv.columns:
        raise ValueError(f"TSV is missing score column {args.score_column!r}")
    if args.append_protein_name and "blast_reference_id" not in tsv.columns:
        raise ValueError("TSV is missing the 'blast_reference_id' column")

    hash_map = (
        short_hashes(parse_prediction_id(old_id).root for old_id in tsv["id"])
        if args.rebase
        else {}
    )

    scores = tsv[args.score_column] if args.append_orf_score else repeat("")
    proteins = tsv["blast_reference_id"] if args.append_protein_name else repeat("")

    rename_map = {}
    generated_ids = set()
    for old_id, score, protein in zip(tsv["id"], scores, proteins):
        parsed = parse_prediction_id(old_id)
        name_parts = list(prefix_parts)

        if args.rebase:
            name_parts.extend(["xORF", hash_map[parsed.root]])
        else:
            name_parts.append(parsed.root)

        name_parts.extend(tags)
        if args.append_orf_score and score != "":
            name_parts.append(f"SC{format_score(score)}")
        if parsed.is_duplicate:
            name_parts.append("DU")
        name_parts.append(f"OR{parsed.orf_number}")
        if parsed.translation_index is not None:
            name_parts.append(f"TI")
        if parsed.inner_number is not None:
            name_parts.append(f"IN{parsed.inner_number}")

        name = "#".join(name_parts)
        if protein:
            name += f"@{protein}"

        if name in generated_ids:
            raise ValueError(f"Renaming produced a duplicate ID: {name!r}")
        generated_ids.add(name)
        rename_map[old_id] = name
    del generated_ids

    bed["id"] = bed["id"].map(rename_map)
    tsv["id"] = tsv["id"].map(rename_map)

    bed.to_csv(args.output_bed, sep="\t", header=False, index=False)
    tsv.to_csv(args.output_tsv, sep="\t", index=False)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, pd.errors.ParserError) as error:
        sys.exit(f"ERROR: {error}")
