#!/usr/bin/env python3

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.1"

import re

# =============================================================================
# REFERENCE INTRON COLUMNS AND CONFIGURATIONS
# =============================================================================

REF_INTRONS_COLUMNS = [
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

REFERENCE_INTRON_ALL_COLUMNS = [
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

REFERENCE_INTRON_ADD_COLUMNS_WITH_SPLICEAI = [
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

REFERENCE_INTRON_ADD_COLUMNS = [
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

DESCRIPTOR_COLUMNS_TO_VECTOR = [
    "coords_of_retention",
    "is_intron_retained_in_frame",
    "location_of_retention",
    "retention_donor_score",
    "retention_acceptor_score",
    "retention_support_type",
    "retains_rt_intron_map",
    "has_rt_intron_map",
]

INTRON_METADATA_COLUMNS_FROM_DESCRIPTOR = [
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

RETENTION_COLUMNS_FROM_DESCRIPTOR = [
    "id",
    "coords_of_retention",
    "is_intron_retained_in_frame",
    "location_of_retention",
    "retention_donor_score",
    "retention_acceptor_score",
    "retention_support_type",
]

RETAINS_RT_COLUMNS_FROM_DESCRIPTOR = ["id", "retains_rt_intron_map"]
HAS_RT_COLUMNS_FROM_DESCRIPTOR = ["id", "has_rt_intron_map"]

RETENTION_HTML_COL = "R_retentions_html"
RETAINS_RT_HTML_COL = "R_retains_rt_intron_html"
HAS_RT_HTML_COL = "R_rt_html"
INTRON_METADATA_HTML_COL = "R_metadata_html"
INTRON_READ_STATUS_COL = "R_read_status"
RETENTION_NESTED_HINT_COLUMN = "R_genomic_coords"
RETENTION_COMMON_KEY = "id"
INTRON_STATUS_COL = "R_read_intron_status"
HAS_RT_INTRON = "HAS_RT_INTRON"
HAS_RT_BUCKET = "R_has_rt"

RETENTION_METADATA_EXCLUDE = ["id", "R_read_status"]
RETENTION_NESTED_EXCLUDE = ["id"]

INTRON_MOD_METADATA_RENAME = {
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

INTRON_MOD_RETENTION_SUBMOD_RENAME = {
    "coords_of_retention": "R_genomic_coords",
    "is_intron_retained_in_frame": "R_intron_in_frame",
    "location_of_retention": "R_gene_location",
    "retention_acceptor_score": "R_acceptor_score",
    "retention_donor_score": "R_donor_score",
    "retention_support_type": "R_support_type",
}

INTRON_MOD_RETAINS_RT_SUBMOD_RENAME = {"retains_rt_intron_map": "R_genomic_coords"}
INTRON_MOD_HAS_RT_SUBMOD_RENAME = {"has_rt_intron_map": "R_genomic_coords"}

# =============================================================================
# POLY-A TAIL ANALYSIS
# =============================================================================

POLYA_METADATA_COLUMNS_FROM_DESCRIPTOR = [
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

POLYA_MOD_METADATA_RENAME = {
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

POLYA_METADATA_HTML_COL = "P_metadata_html"
POLYA_METADATA_EXCLUDE = ["id", "P_read_status"]
POLYA_RESULT_COLS = ["id", "P_read_status", "P_metadata_html"]

# =============================================================================
# TRUNCATION ANALYSIS
# =============================================================================

TRUNCATION_METADATA_COLUMNS_FROM_DESCRIPTOR = [
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

TRUNCATION_MOD_METADATA_RENAME = {
    "component_truncation_ratio": "T_component_truncation_ratio",
    "is_dirty_utr_component": "T_was_component_reviewed",
    "is_novel_start": "T_is_novel_start",
    "is_read_truncated": "T_read_status",
    "is_truncation_supported": "T_is_truncation_supported",
    "query_utr_component_size": "T_query_utr_component_size",
    "ref_utr_component_size": "T_ref_utr_component_size",
    "truncation_support_ratio": "T_truncation_support_ratio",
}

TRUNCATION_METADATA_HTML_COL = "T_metadata_html"
TRUNCATION_METADATA_EXCLUDE = ["id", "T_read_status"]
TRUNCATION_RESULT_COLS = ["id", "T_read_status", "T_metadata_html"]

# =============================================================================
# ORF (OPEN READING FRAME) PREDICTION
# =============================================================================

ORF_PREDICTOR_COLS = [
    "O_blast_id",
    "O_read_orf_score",
    "O_class_0_prob",  # model
    "O_blast_pid",
    "O_blast_e_value",
    "O_blast_offset",
    "O_blast_percentage_aligned",  # Blast
    "O_tai_start_score",
    "O_tai_stop_score",  # Tai
    "O_toga_masked",
    "O_toga_label",
    "O_toga_pid",
    "O_toga_blosum",
    "O_toga_overlap_bp",  # TOGA
]

ORF_MOD_METADATA_RENAME = {
    "blast_id": "O_blast_id",
    "class1_prob": "O_read_orf_score",
    "class0_prob": "O_class_0_prob",
    "blast_pid": "O_blast_pid",
    "blast_e-value": "O_blast_e_value",
    "blast_offset": "O_blast_offset",
    "blast_percentage_aligned": "O_blast_percentage_aligned",
    "translationAI_orf_start": "O_tai_start_score",
    "translationAI_orf_stop": "O_tai_stop_score",
    "masked": "O_toga_masked",
    "toga_label": "O_toga_label",
    "toga_pid": "O_toga_pid",
    "toga_blosum": "O_toga_blosum",
    "toga_overlap_bp": "O_toga_overlap_bp",
}

ORF_METADATA_HTML_COL = "O_metadata_html"
ORF_METADATA_EXCLUDE = ["id", "O_read_orf_score"]
ORF_RESULT_COLS = ["O_blast_id", "O_read_orf_score", "O_metadata_html"]

# =============================================================================
# TAG PARSING AND READ METADATA
# =============================================================================

TAG_MAP = {
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
}

TAG_REGEX = re.compile(r"([A-Z]{2})(\d+)?")

# =============================================================================
# DATA TYPE MAPPINGS AND TRANSFORMATIONS
# =============================================================================

FRAME_MAPPING = {
    True: "IN_FRAME",
    False: "OUT_OF_FRAME",
    "true": "IN_FRAME",
    "false": "OUT_OF_FRAME",
}

LOCATION_MAPPING = {
    "Mixed": "UTR-CDS",
    "CDS": "CDS",
    "UTR": "UTR",
}

EXON_STATUS_MAPPING = {
    "KEEP": "PASS",
    "DISCARD": "RETAINS_SPLICED_INTRON",
    "UNCLEAR": "UNCLEAR",
}

INTRON_STATUS_MAPPING = {"KEEP": "PASS", "DISCARD": "HAS_RT_INTRON"}

RETAINS_RT_MAPPING = {
    True: "RETAINS_RT_INTRON",
    False: "NULL",
    "True": "RETAINS_RT_INTRON",
    "False": "NULL",
}

RETAINS_ANY_MAPPING = {
    True: "RETAINS_INTRON",
    False: "NO_RETENTION",
    "True": "RETAINS_INTRON",
    "False": "NO_RETENTION",
}

IS_INTRAPRIMMING_MAPPING = {
    True: "INTRAPRIMMING",
    False: "PASS",
    "true": "INTRAPRIMMING",
    "false": "PASS",
}

IS_POLYA_SUPPORTED_MAPPING = {
    True: "SUPPORTED",
    False: "NOT_SUPPORTED",
    "true": "SUPPORTED",
    "false": "NOT_SUPPORTED",
}

IS_TRUNCATED_MAPPING = {
    True: "TRUNCATED",
    False: "PASS",
    "true": "TRUNCATED",
    "false": "PASS",
}

IS_TRUNCATION_SUPPORTED_MAPPING = {
    True: "SUPPORTED",
    False: "NOT_SUPPORTED",
    "true": "SUPPORTED",
    "false": "NOT_SUPPORTED",
}

DESCRIPTOR_MAPPING_TYPES = {
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

DESCRIPTOR_MAPPING_MISSING_VALUES = {
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

FLOAT_TO_INT_COLS = [
    "number_of_retentions",
    "query_intron_component_size",
    "ref_introns_component_size",
]

# =============================================================================
# BED12 COLS, FINAL OUTPUT, BIGBED SCHEMA AND BINARIES
# =============================================================================

BED_TO_BIG_BED = "../../c/bedToBigBed"
SCHEMA_SQL = "../../schema/schema.as"
READ_DECIDE_COLUMNS = ["id", "R_read_status", "P_read_status", "T_read_status"]

MAP_IGNORE_CATEGORY = {
    "rtintron": "R_has_rt",
    "retention": "R_read_status",
    "intraprimming": "P_read_status",
    "truncation": "T_read_status",
}

MAPPING_KEYS_TO_READ_NAME = {
    "R_has_rt": "rt_reads.bed",
    "pass": "pass.bed",
    "trash": "trash.bed",
    "R_read_status": "retention.bed",
    "P_read_status": "intraprimming.bed",
    "T_read_status": "truncation.bed",
}

BED12_COLS = [
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

BB_SCHEMA_COLUMNS = [
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
]
