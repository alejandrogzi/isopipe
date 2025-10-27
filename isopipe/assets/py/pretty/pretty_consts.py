#!/usr/bin/env python3

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.2"

import re
from typing import Dict, List, Union

# =============================================================================
# REFERENCE INTRON COLUMNS AND CONFIGURATIONS
# =============================================================================

REF_INTRONS_COLUMNS: List[str] = [
    "chrom",  # Chromosome identifier
    "start",  # Start position of the intron (1-based)
    "end",  # End position of the intron (1-based)
    "strand",  # Genomic strand (+ or -)
    "seen",  # Frequency of how many reads contain this intron
    "spanned",  # Frequency of how many reads span this intron
    "splice_ai_donor",  # SpliceAI score for the donor site
    "splice_ai_acceptor",  # SpliceAI score for the acceptor site
    "max_ent_donor",  # MaxEntScan score for the donor site
    "max_ent_acceptor",  # MaxEntScan score for the acceptor site
    "donor_sequence",  # Nucleotide sequence around the donor site
    "acceptor_sequence",  # Nucleotide sequence around the acceptor site
    "donor_context",  # MaxEntScan 9-mer donor context sequence
    "acceptor_context",  # MaxEntScan 23-mer acceptor context sequence
    "intron_position",  # Classification of the intron's position according to TOGA
    "is_toga_supported",  # Boolean indicating if the intron is supported by TOGA
    "is_in_frame",  # Boolean indicating if the intron maintains the reading frame
    "donor_rt_context",  # RT-switch context sequence for the donor site
    "acceptor_rt_context",  # RT-switch context sequence for the acceptor site
    "is_rt_intron",  # Boolean indicating if the intron is an RT-switch intron
    "is_nag_intron",  # Boolean indicating if the intron is a TOGA-nag intron
    "support",  # Classification of the intron's support type
]

REFERENCE_INTRON_ALL_COLUMNS: List[str] = [
    "strand",
    "seen",
    "spanned",
    "splice_ai_donor",
    "splice_ai_acceptor",
    "max_ent_donor",
    "max_ent_acceptor",
    "donor_sequence",
    "acceptor_sequence",
    "donor_context",
    "acceptor_context",
    "intron_position",
    "is_toga_supported",
    "is_in_frame",
    "donor_rt_context",
    "acceptor_rt_context",
    "is_rt_intron",
    "is_nag_intron",
    "support",
]

REFERENCE_INTRON_ADD_COLUMNS_WITH_SPLICEAI: List[str] = [
    "strand",
    "seen",
    "spanned",
    "splice_ai_donor",
    "splice_ai_acceptor",
    "max_ent_donor",
    "max_ent_acceptor",
    "donor_sequence",
    "acceptor_sequence",
    "donor_context",
    "acceptor_context",
    "is_toga_supported",
    "is_nag_intron",
]

REFERENCE_INTRON_ADD_COLUMNS: List[str] = [
    "strand",
    "seen",
    "spanned",
    "max_ent_donor",
    "max_ent_acceptor",
    "donor_sequence",
    "acceptor_sequence",
    "donor_context",
    "acceptor_context",
    "is_toga_supported",
    "is_nag_intron",
]

# =============================================================================
# INTRON RETENTION ANALYSIS
# =============================================================================

DESCRIPTOR_COLUMNS_TO_VECTOR: List[str] = [
    "coords_of_retention",
    "is_intron_retained_in_frame",
    "location_of_retention",
    "retention_donor_score",
    "retention_acceptor_score",
    "retention_support_type",
    "retains_rt_intron_map",
    "has_rt_intron_map",
]

INTRON_METADATA_COLUMNS_FROM_DESCRIPTOR: List[str] = [
    "id",
    "exonic_status",
    "intronic_status",
    "intron_retention",
    "retains_rt_intron",
    "number_of_retentions",
    "component_size",
    "component_retention_ratio",
    "is_dirty_intron_component",
    "query_intron_component_size",
    "ref_introns_component_size",
]

RETENTION_COLUMNS_FROM_DESCRIPTOR: List[str] = [
    "id",
    "coords_of_retention",
    "is_intron_retained_in_frame",
    "location_of_retention",
    "retention_donor_score",
    "retention_acceptor_score",
    "retention_support_type",
]

RETAINS_RT_COLUMNS_FROM_DESCRIPTOR: List[str] = ["id", "retains_rt_intron_map"]
HAS_RT_COLUMNS_FROM_DESCRIPTOR: List[str] = ["id", "has_rt_intron_map"]

RETENTION_HTML_COL: str = "R_retentions_html"
RETAINS_RT_HTML_COL: str = "R_retains_rt_intron_html"
HAS_RT_HTML_COL: str = "R_rt_html"
INTRON_METADATA_HTML_COL: str = "R_metadata_html"
INTRON_READ_STATUS_COL: str = "R_read_status"
RETENTION_NESTED_HINT_COLUMN: str = "R_genomic_coords"
RETENTION_COMMON_KEY: str = "id"
INTRON_STATUS_COL: str = "R_read_intron_status"
HAS_RT_INTRON: str = "HAS_RT_INTRON"
HAS_RT_BUCKET: str = "R_has_rt"

RETENTION_METADATA_EXCLUDE: List[str] = ["id", "R_read_status"]
RETENTION_NESTED_EXCLUDE: List[str] = ["id"]

INTRON_MOD_METADATA_RENAME: Dict[str, str] = {
    "exonic_status": "R_read_exon_status",
    "intronic_status": "R_read_intron_status",
    "retains_rt_intron": "R_read_rt_status",
    "intron_retention": "R_read_retention_status",
    "number_of_retentions": "R_number_of_retentions",
    "component_size": "R_component_size",
    "component_retention_ratio": "R_component_intron_retention_ratio",
    "is_dirty_intron_component": "R_was_component_reviewed",
    "query_intron_component_size": "R_number_of_query_introns",
    "ref_introns_component_size": "R_number_of_reference_introns",
}

INTRON_MOD_RETENTION_SUBMOD_RENAME: Dict[str, str] = {
    "coords_of_retention": "R_genomic_coords",
    "is_intron_retained_in_frame": "R_intron_in_frame",
    "location_of_retention": "R_gene_location",
    "retention_acceptor_score": "R_acceptor_score",
    "retention_donor_score": "R_donor_score",
    "retention_support_type": "R_support_type",
}

INTRON_MOD_RETAINS_RT_SUBMOD_RENAME: Dict[str, str] = {
    "retains_rt_intron_map": "R_genomic_coords"
}
INTRON_MOD_HAS_RT_SUBMOD_RENAME: Dict[str, str] = {
    "has_rt_intron_map": "R_genomic_coords"
}

# =============================================================================
# POLY-A TAIL ANALYSIS
# =============================================================================

POLYA_METADATA_COLUMNS_FROM_DESCRIPTOR: List[str] = [
    "id",
    "forced_poly_a",
    "genomic_poly_a",
    "intrapriming_comp_ratio",
    "is_dirty_polya_component",
    "is_intrapriming",
    "is_poly_a_supported",
    "poly_a_location",
    "poly_a_score",
    "truncation_support_ratio",
    "whole_poly_a_length",
]

POLYA_MOD_METADATA_RENAME: Dict[str, str] = {
    "is_intrapriming": "P_read_status",
    "forced_poly_a": "P_is_forced_poly_a",
    "genomic_poly_a": "P_genomic_poly_a_len",
    "intrapriming_comp_ratio": "P_intrapriming_comp_ratio",
    "is_dirty_polya_component": "P_was_component_reviewed",
    "is_poly_a_supported": "P_is_poly_a_supported",
    "poly_a_location": "P_poly_a_location_in_toga",
    "poly_a_score": "P_poly_a_score",
    "truncation_support_ratio": "P_truncation_support_ratio",
    "whole_poly_a_length": "P_whole_poly_a_length",
}

POLYA_METADATA_HTML_COL: str = "P_metadata_html"
POLYA_METADATA_EXCLUDE: List[str] = ["id", "P_read_status"]
POLYA_RESULT_COLS: List[str] = ["id", "P_read_status", "P_metadata_html"]

# =============================================================================
# TRUNCATION ANALYSIS
# =============================================================================

TRUNCATION_METADATA_COLUMNS_FROM_DESCRIPTOR: List[str] = [
    "id",
    "component_truncation_ratio",
    "is_dirty_utr_component",
    "is_novel_start",
    "is_read_truncated",
    "is_truncation_supported",
    "query_utr_component_size",
    "ref_utr_component_size",
    "truncation_support_ratio",
]

TRUNCATION_MOD_METADATA_RENAME: Dict[str, str] = {
    "component_truncation_ratio": "T_component_truncation_ratio",
    "is_dirty_utr_component": "T_was_component_reviewed",
    "is_novel_start": "T_is_novel_start",
    "is_read_truncated": "T_read_status",
    "is_truncation_supported": "T_is_truncation_supported",
    "query_utr_component_size": "T_query_utr_component_size",
    "ref_utr_component_size": "T_ref_utr_component_size",
    "truncation_support_ratio": "T_truncation_support_ratio",
}

TRUNCATION_METADATA_HTML_COL: str = "T_metadata_html"
TRUNCATION_METADATA_EXCLUDE: List[str] = ["id", "T_read_status"]
TRUNCATION_RESULT_COLS: List[str] = ["id", "T_read_status", "T_metadata_html"]

# =============================================================================
# ORF (OPEN READING FRAME) PREDICTION
# =============================================================================

ORF_PREDICTOR_COLS: List[str] = [
    "O_id",
    "O_tai_start_score",
    "O_tai_stop_score",
    "O_relative_orf_start",
    "O_relative_orf_end",
    "O_start_codon",
    "O_stop_codon",
    "O_inner_stops",
    "O_orf_type",
    "O_blast_pid",
    "O_blast_evalue",
    "O_blast_offset",
    "O_blast_length",
    "O_blast_percentage_aligned",
    "O_rna_score",
    "O_tai_mean_score",
    "O_log_orf_len",
    "O_predicted_class",
    "O_read_orf_score",
]


ORF_MOD_METADATA_RENAME: Dict[str, str] = {
    "chr": "O_chr",
    "start": "O_start",
    "end": "O_end",
    "id": "O_id",
    "strand": "O_strand",
    "tai_start_score": "O_tai_start_score",
    "tai_stop_score": "O_tai_stop_score",
    "relative_orf_start": "O_relative_orf_start",
    "relative_orf_end": "O_relative_orf_end",
    "start_codon": "O_start_codon",
    "stop_codon": "O_stop_codon",
    "inner_stops": "O_inner_stops",
    "orf_type": "O_orf_type",
    "blast_pid": "O_blast_pid",
    "blast_evalue": "O_blast_evalue",
    "blast_offset": "O_blast_offset",
    "blast_length": "O_blast_length",
    "blast_percentage_aligned": "O_blast_percentage_aligned",
    "key": "O_key",
    "rna_score": "O_rna_score",
    "tai_mean_score": "O_tai_mean_score",
    "blast_match": "O_blast_match",
    "log_orf_len": "O_log_orf_len",
    "has_canonical_start": "O_has_canonical_start",
    "has_canonical_stop": "O_has_canonical_stop",
    "neg_log10_blast_evalue": "O_neg_log10_blast_evalue",
    "tai_combined_score": "O_tai_combined_score",
    "predicted_class": "O_predicted_class",
    "prob_noncoding": "O_prob_noncoding",
    "prob_coding": "O_read_orf_score",
}

ORF_METADATA_HTML_COL: str = "O_metadata_html"
ORF_METADATA_EXCLUDE: List[str] = ["O_id", "O_read_orf_score"]
ORF_RESULT_COLS: List[str] = ["O_id", "O_read_orf_score", "O_metadata_html"]

ORF_TYPE_MAPPING: Dict[Union[str, int], str] = {
    1: "COMPLETE",
    2: "COMPLETE_NESTED",
    3: "TAI_UNIQUE",
    4: "FIVE_PRIME_PARTIAL",
    5: "THREE_PRIME_PARTIAL",
    6: "THREE_PRIME_PARTIAL_NESTED",
    7: "FIVE_PRIME_PARTIAL_NESTED",
    "1": "COMPLETE",
    "2": "COMPLETE_NESTED",
    "3": "TAI_UNIQUE",
    "4": "FIVE_PRIME_PARTIAL",
    "5": "THREE_PRIME_PARTIAL",
    "6": "THREE_PRIME_PARTIAL_NESTED",
    "7": "FIVE_PRIME_PARTIAL_NESTED",
}


# =============================================================================
# TAG PARSING AND READ METADATA
# =============================================================================

TAG_MAP: Dict[str, str] = {
    "FC": "TAG_five_clip_len",
    "TC": "TAG_three_clip_len",
    "PA": "TAG_polya_len",
    "PR": "TAG_polya_read_len",
    "IY": "TAG_mapping_identity",
    "SG": "TAG_singleton",
    "OR": "TAG_orf_number",
    "NE": "TAG_nested_orf_number",
    "UT": "TAG_unique_tai",
    "FK": "TAG_fake_fusion",
    "RW": "TAG_review_fusion",
    "SN": "TAG_strong_nmd",
    "WN": "TAG_weak_nmd",
    "DU": "TAG_duplicated_non_best_orf",
}

TAG_REGEX: re.Pattern = re.compile(r"([A-Z]{2})(\d+)?")

# =============================================================================
# DATA TYPE MAPPINGS AND TRANSFORMATIONS
# =============================================================================

FRAME_MAPPING: Dict[Union[bool, str], str] = {
    True: "IN_FRAME",
    False: "OUT_OF_FRAME",
    "true": "IN_FRAME",
    "false": "OUT_OF_FRAME",
}

LOCATION_MAPPING: Dict[Union[bool, str], str] = {
    "Mixed": "UTR-CDS",
    "CDS": "CDS",
    "UTR": "UTR",
}

EXON_STATUS_MAPPING: Dict[Union[bool, str], str] = {
    "KEEP": "PASS",
    "DISCARD": "RETAINS_SPLICED_INTRON",
    "UNCLEAR": "UNCLEAR",
}

INTRON_STATUS_MAPPING: Dict[Union[bool, str], str] = {
    "KEEP": "PASS",
    "DISCARD": "HAS_RT_INTRON",
}

RETAINS_RT_MAPPING: Dict[Union[bool, str], str] = {
    True: "RETAINS_RT_INTRON",
    False: "NULL",
    "True": "RETAINS_RT_INTRON",
    "False": "NULL",
}

RETAINS_ANY_MAPPING: Dict[Union[bool, str], str] = {
    True: "RETAINS_INTRON",
    False: "NO_RETENTION",
    "True": "RETAINS_INTRON",
    "False": "NO_RETENTION",
}

IS_INTRAPRIMMING_MAPPING: Dict[Union[bool, str], str] = {
    True: "INTRAPRIMMING",
    False: "PASS",
    "true": "INTRAPRIMMING",
    "false": "PASS",
}

IS_POLYA_SUPPORTED_MAPPING: Dict[Union[bool, str], str] = {
    True: "SUPPORTED",
    False: "NOT_SUPPORTED",
    "true": "SUPPORTED",
    "false": "NOT_SUPPORTED",
}

IS_TRUNCATED_MAPPING: Dict[Union[bool, str], str] = {
    True: "TRUNCATED",
    False: "PASS",
    "true": "TRUNCATED",
    "false": "PASS",
}

IS_TRUNCATION_SUPPORTED_MAPPING: Dict[Union[bool, str], str] = {
    True: "SUPPORTED",
    False: "NOT_SUPPORTED",
    "true": "SUPPORTED",
    "false": "NOT_SUPPORTED",
}

DESCRIPTOR_MAPPING_TYPES: Dict[str, Dict[Union[bool, str], str]] = {
    "is_intron_retained_in_frame": FRAME_MAPPING,
    "location_of_retention": LOCATION_MAPPING,
    "exonic_status": EXON_STATUS_MAPPING,
    "intronic_status": INTRON_STATUS_MAPPING,
    "retains_rt_intron": RETAINS_RT_MAPPING,
    "intron_retention": RETAINS_ANY_MAPPING,
    "is_intrapriming": IS_INTRAPRIMMING_MAPPING,
    "is_poly_a_supported": IS_POLYA_SUPPORTED_MAPPING,
    "is_read_truncated": IS_TRUNCATED_MAPPING,
    "is_truncation_supported": IS_TRUNCATION_SUPPORTED_MAPPING,
}

# =============================================================================
# DATA PROCESSING CONFIGURATIONS
# =============================================================================

DESCRIPTOR_MAPPING_MISSING_VALUES: Dict[str, Union[bool, str, int, float]] = {
    # Categorical/string columns
    "intronic_status": "PASS",
    "exonic_status": "PASS",
    "retains_rt_intron": "NULL",
    "intron_retention": "NO_RETENTION",
    # Numeric columns
    "number_of_retentions": 0,
    "component_retention_ratio": 0.0,
    "query_intron_component_size": 0,
    "ref_introns_component_size": 0,
    # Boolean column
    "is_dirty_intron_component": False,
}

FLOAT_TO_INT_COLS: List[str] = [
    "number_of_retentions",
    "query_intron_component_size",
    "ref_introns_component_size",
]

# =============================================================================
# BED12 COLS, FINAL OUTPUT, BIGBED SCHEMA AND BINARIES
# =============================================================================

BED_TO_BIG_BED: str = "../../c/bedToBigBed"
SCHEMA_SQL: str = "../../schema/schema.as"
READ_DECIDE_COLUMNS: List[str] = [
    "id",
    "R_read_status",
    "P_read_status",
    "T_read_status",
]

MAP_IGNORE_CATEGORY: Dict[str, str] = {
    "rtintron": "R_has_rt",
    "retention": "R_read_status",
    "intraprimming": "P_read_status",
    "truncation": "T_read_status",
}

MAPPING_KEYS_TO_READ_NAME: Dict[str, str] = {
    "R_has_rt": "rt_reads.bed",
    "pass": "pass.bed",
    "trash": "trash.bed",
    "R_read_status": "retention.bed",
    "P_read_status": "intraprimming.bed",
    "T_read_status": "truncation.bed",
}

BED12_COLS: List[str] = [
    "chr",
    "start",
    "end",
    "id",
    "score",
    "strand",
    "cds_start",
    "cds_end",
    "rgb",
    "exon_count",
    "exon_sizes",
    "exon_starts",
]

BB_SCHEMA_COLUMNS: List[str] = [
    "id",
    "R_read_status",
    "R_read_intron_status",
    "R_metadata_html",
    "R_retentions_html",
    "R_retains_rt_intron_html",
    "R_rt_html",
    "P_read_status",
    "P_metadata_html",
    "T_read_status",
    "T_metadata_html",
    "O_read_orf_score",
    "O_metadata_html",
    "TAG_five_clip_len",
    "TAG_three_clip_len",
    "TAG_polya_len",
    "TAG_polya_read_len",
    "TAG_mapping_identity",
    "TAG_singleton",
    "TAG_orf_number",
    "TAG_nested_orf_number",
    "TAG_fake_fusion",
    "TAG_review_fusion",
    "TAG_strong_nmd",
    "TAG_weak_nmd",
    "TAG_unique_tai",
    "TAG_duplicated_non_best_orf",
]
