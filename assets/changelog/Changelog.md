# Changelog

## [v2.0.25] - 2026-07-29

This release replaces the previous two-entrypoint model (`isoseq` / `map`) with three granular modes — `subreads`, `ccs`, and `flnc` — that map directly to the real PacBio processing stages users are working with. The new model gives you control over where your data enters the pipeline: raw subreads, already-called CCS reads, or full-length non-chimeric reads that have already gone through primer removal and clustering. The internal ISOSEQ subworkflow was refactored to branch on these entrypoints, running PBCCS and chunk merging only when needed (i.e., for `subreads` and `ccs`), while `flnc` reads skip straight to LIMA and refinement. The default entrypoint was changed from `isoseq` to `subreads`, and the workflow banner now prints the active entrypoint at launch to make it immediately visible which mode is in use. Additional structural comments were added throughout the main workflow and the ISOSEQ subworkflow to clarify the stage boundaries in the code.

### New entrypoint model

- Replaced the previous `isoseq` and `map` entrypoints with three new options: `subreads`, `ccs`, and `flnc`. Validation in `main.nf` was updated to reject unknown values with a clear error message listing the valid options.
- Changed the default `entrypoint` parameter from `"isoseq"` to `"subreads"` in `nextflow.config` and `params.json`.
- The workflow banner now displays `Entrypoint: ${params.entrypoint}` alongside the input, output, and genome paths, making it easier to confirm which mode is active at a glance.

### ISOSEQ subworkflow refactoring

- The `ISOSEQ` subworkflow now accepts an `entrypoint` input parameter (expected values: `ccs` or `flnc`; `subreads` is handled upstream). A `switch` statement dispatches execution into two branches:
  - **`ccs`**: runs `PBCCS` on chunked BAMs to generate circular consensus sequences, groups chunks by parent sample via `groupTuple`, merges them with `PBMERGE`, and passes the merged BAMs to `LIMA` for primer removal.
  - **`flnc`**: passes the raw input BAMs directly to `LIMA`, skipping CCS generation and merging entirely.
- An explicit error is thrown if an unrecognized entrypoint reaches the subworkflow.

### Preprocessing dispatch update

- The preprocessing subworkflow now checks for `subreads` or `ccs` to route into the `ISOSEQ` branch (the `ISOSEQ` call now also passes the `entrypoint` value), and checks for `flnc` to enter the direct FASTQ-reading branch (formerly the `map` branch).

### Codebase organization

- Added ASCII-delimited section headers throughout `workflows/ark.nf`, `subworkflows/isoseq/main.nf`, and `subworkflows/preprocessing/main.nf` to visually separate stages such as Autosql, Preprocessing, Alignment, ORF calling, NMD calling, Polishing, and Tracking. These are purely cosmetic and have no effect on pipeline logic.

## [v2.0.24] - 2026-07-24

This release consolidates the pipeline's alignment architecture under a single unified subworkflow, renames the project from isopipe to ARK, and introduces FLAIR as a sixth aligner option. The previously separate desalt, ultra, and pbmm2 alignment subworkflows have been merged into a single `split_align` subworkflow that dispatches to the correct backend via a switch statement, significantly reducing code duplication. Additional work includes a new ORF renaming engine for xORF, programmatic AutoSQL schema generation, and several channel-wiring corrections across downstream subworkflows.

### Pipeline rename to ARK

- Renamed the primary workflow from `ISOPIPE` to `ARK` and the corresponding entry file from `workflows/isopipe.nf` to `workflows/ark.nf`. The pipeline description was updated to "A Reference pipeline to annotate euKaryotes at high resolution" and all manifest metadata (name, homePage) now point to the `alejandrogzi/ark` repository.
- Updated `main.nf` to include `ARK` from the new workflow path and to log the Ark banner (version, description, authors, lab) instead of the former isopipe banner.

### Unified alignment architecture

- Removed three dedicated alignment subworkflows — `desalt_align`, `ultra_align`, and `pbmm2_align` — and consolidated all six supported aligners into the single `split_align` subworkflow. The subworkflow now accepts an `aligner` string parameter and uses a switch statement to route chunked reads through the correct alignment module, BAM conversion, and multi-sample merge path.
- Added dedicated index modules for each aligner: `ARK_INDEX` and `FLAIR_INDEX` (both minimap2-based with appropriate `-x` presets), alongside the existing `MINIMAP2_INDEX`, `ULTRA_INDEX`, `DESALT_INDEX`, and `PBMM2_INDEX`. Each aligner now has its own publish directory under `06_<ALIGNER>_INDEX`.
- Each aligner's BAM output now routes through a uniquely named samtools conversion process (`SAMTOOLS_BAM_ARK_ALIGN`, `SAMTOOLS_BAM_MINIMAP2_ALIGN`, `SAMTOOLS_BAM_DESALT_ALIGN`, `SAMTOOLS_BAM_PBMM2_ALIGN`, `SAMTOOLS_BAM_FLAIR_ALIGN`) to prevent channel collisions in multi-sample merging.
- Multi-sample BAM merging was similarly split into per-aligner `SAMTOOLS_MERGE_BAM_MULTI_SAMPLE_*` processes, ensuring correct version tracking when multiple aligners are exercised across pipeline invocations.

### New aligner: FLAIR

- Added **FLAIR** as a selectable aligner (`aligner = "flair"`). The integration includes a `FLAIR_INDEX` module that builds the minimap2 index with the `-x splice` preset and a `FLAIR_ALIGN` process aliased from `MINIMAP2_ALIGN` that runs with the standard splice preset. The `split_align` subworkflow routes FLAIR-aligned output through a dedicated `SAMTOOLS_BAM_FLAIR_ALIGN` conversion step.

### Ark as default aligner

- Changed the default `aligner` parameter from `"minimap2"` to `"ark"`. The `"mm2"` aligner runs minimap2 with the PacBio Kinnex standard preset (`-uf -ax splice:hq`), replacing the previous verbose flag set (`-a -c --eqx -uf -C5 -G ... -ax splice:hq --secondary=... -cs ... --junc-bonus ... --junc-pen ...`) that now runs through `"ark"`.
- The second-pass re-alignment (cigar extension + fragment detection) is now gated exclusively behind `params.aligner == 'ark'`, since other aligners handle fragment resolution differently or do not benefit from the additional pass.

### ORF renaming engine

- Added `rename_predictions.py` (v0.0.4), a standalone Python script that rewrites ORF prediction identifiers in BED12/TSV files with human-readable, information-rich names. Supported features include hash-based rebasing of root IDs (`--rebase`), ORF score appending (`--append-orf-score`), protein name appending (`--append-protein-name`), custom prefix injection (`--custom-prefix`), and full deactivation (`--deactivate`).
- Added five new xORF parameters: `xorf_rename_deactivate` (default: `false`), `xorf_rename_rebase` (default: `false`), `xorf_rename_append_orf_score` (default: `true`), `xorf_rename_append_protein_name` (default: `true`), and `xorf_rename_custom_prefix` (default: `null`). These are propagated through both `XORF_PREDICT_ORFS` and `XORF_PREDICT_FUSION_ORFS` process invocations.

### AutoSQL schema generation

- Added `AUTOSQL_BASE` and `AUTOSQL_SCHEMA` Nextflow processes that programmatically generate AutoSQL `.as` schema files at runtime, replacing the previous reliance on static files from `assets/as/`. `AUTOSQL_BASE` emits a minimal BED12-compatible schema, while `AUTOSQL_SCHEMA` emits the extended schema including read-status, metadata, collapsed-reads, and ORF metadata columns. Both are imported in the main ARK workflow and their outputs are wired to the `autosql` and `schema` channels used by downstream bed-to-bigbed conversion.

### Parameter and validation changes

- Renamed the `minimap2` aligner option to `mm2` throughout the parameter schema, documentation, and validation logic. The valid aligner set is now: `mm2`, `ultra`, `desalt`, `pbmm2`, `flair`, `ark`.
- Removed the `ultra_do_second_pass` parameter entirely. The Ultra two-pass logic from v2.0.22 is superseded by the unified second-pass gating in `split_align`, which only activates for the `ark` aligner.
- Improved validation error messages with explicit `ERROR:` prefixes and structured multi-line formatting via `stripIndent()`.

### Channel wiring fixes

- Corrected the annotation channel format in `PREPOLISH` and `POLISH` subworkflows: the annotation is now expected as `[ val(meta), [ annotation ] ]` from preprocessing rather than being re-mapped with `[id:annotation.baseName]` at the subworkflow boundary. This prevents inconsistent channel shapes when the annotation flows through multi-sample paths.
- The `POLISH` subworkflow now receives the `autosql` channel as an explicit input parameter rather than constructing it from a static asset path, ensuring the AutoSQL schema stays in sync with the dynamically generated output.

### Chores

- Added `.gitattributes` to the repository root.
- Updated `params.json` schema documentation to reflect the new aligner list and the `xorf_rename_deactivate` parameter.
- Removed dedicated `ext.args` blocks for `MINIMAP2_ALIGN` fragment-detection processes, replacing them with the simplified Ark preset configuration.

## [v2.0.23] - 2026-07-18

This release expands the pipeline's alignment options with two new backends — deSALT and pbmm2 — giving users four aligners to choose from depending on their sequencing platform and sensitivity requirements. The minimap2 branch also receives several correctness fixes, particularly around second-pass gating and index preset configuration, along with minor resource tuning and a new xORF parameter.

### New aligner backends

- Added **deSALT** as a selectable aligner (`aligner = "desalt"`), a splice-aware long-read aligner well-suited for isoform discovery. The integration includes a new `DESALT_INDEX` module for building the deSALT genome index, a `DESALT_ALIGN` module for alignment with configurable annotation guidance via `desalt_use_annotation`, and a dedicated `desalt_align` subworkflow that handles chunking, alignment, BAM conversion, multi-sample merging, adapter removal, cigar extension, polyA segmentation, twin collapse, and fusion detection.
- Added **pbmm2** as a selectable aligner (`aligner = "pbmm2"`), PacBio's official minimap2 wrapper optimized for CCS reads. The integration includes a `PBMM2_INDEX` module that builds the index using the `splice` preset, a `PBMM2_ALIGN` module using the `ISOSEQ` preset, and a dedicated `pbmm2_align` subworkflow covering the full post-alignment processing pipeline.
- Updated the `aligner` parameter documentation and schema to reflect all four supported options: `minimap2`, `ultra`, `desalt`, and `pbmm2`.
- Extended the preprocessing subworkflow with conditional index-building branches for deSALT and pbmm2, accepting new `desalt_index` and `pbmm2_index` parameters to support pre-built index paths.

### Minimap2 indexing and second-pass fixes

- Corrected the `MINIMAP2_INDEX` process configuration by adding `ext.args` with the appropriate `-x` splice alignment preset, ensuring the minimap2 index is built with the correct parameters for transcript-aware alignment.
- Fixed `ext.use_junc_bed` assignment for both `MINIMAP2_ALIGN` and `MINIMAP2_ALIGN_FRAGMENTS` to properly reflect `!params.minimap2_align_use_splice_scores`, preventing the junction BED from being applied in scenarios where splice scores should take precedence.
- Gated the entire fragment detection and second-pass re-alignment logic in the `split_align` subworkflow behind the new `minimap2_align_do_second_pass` parameter (default: `true`). When disabled, the pipeline skips fragment detection and minimap2 re-alignment entirely, proceeding directly to segmentation with the primary BAM. This mirrors the behavior already present in the Ultra backend and avoids unnecessary computation when the second pass is not needed.

### xORF improvements

- Added the `xorf_skip_netstart` parameter (default: `true`) to allow users to bypass NetStart-based start codon predictions during ORF calling, streamlining results when only canonical start codons are of interest.

### Resource tuning

- Reduced the `process_medium_high_memory` resource label from 6 CPUs and 96 GB to 4 CPUs and 36 GB per task attempt, providing a more cost-effective baseline for medium-memory workloads without compromising stability.

### Chores

- Bumped `isotools` and `xorf` submodules to their latest versions.

## [v2.0.22] - 2026-06-19

This release introduces a two-pass alignment strategy for the Ultra backend that significantly improves sensitivity for fragmented reads, along with finer-grained control over minimap2 junction scoring during splice site annotation. Several channel wiring bugs introduced in v2.0.21 have also been corrected.

### Ultra two-pass alignment

- Implemented a second pass in the Ultra alignment branch: after the initial Ultra alignment, reads undergo optional cigar extension and fragment detection via isotools, and the resulting fragments are re-aligned using minimap2 for higher-resolution junction placement. This is controlled by the new `ultra_do_second_pass` parameter (default: `true`).
- Added the `ultra_max_intron_size` parameter (default: 300 kb) to configure the maximum intron length passed to the Ultra aligner.
- Added new process modules for the ultra second pass: `MINIMAP2_ALIGN_FRAGMENTS_ULTRA`, `SAMTOOLS_BAM_FRAGMENTS_ULTRA`, and `ISOTOOLS_FIND_FRAGMENTS_ULTRA`, each with dedicated process configuration in `nextflow.config`.
- Updated the preprocessing subworkflow to conditionally build a minimap2 index when Ultra is selected and `ultra_do_second_pass` is enabled, ensuring the second pass has the required index available.
- Changed the Ultra merged BAM output directory from `06_ULTRA_ALIGN/MERGED` to `06_ULTRA_ALIGN/ULTRA_MERGED` to avoid ambiguity with other merge paths, and added `FRAGMENTS` publish directories for fragment-level SAM and BAM outputs.
- The Ultra second pass fragments now flow through `ISOTOOLS_SEGMENT_POLYA_FRAGMENTS` and are mixed with the primary segmented output for downstream processing.

### Minimap2 scoring granularity

- The `--junc-bed` flag in `MINIMAP2_ALIGN` is now conditionally applied based on the new `ext.use_junc_bed` configuration value, allowing the pipeline to omit the junction bed when splice scores are available and should take precedence.
- The `SPLICEAI_DERIVE` process arguments (`--include-ss-from-regions`, `--position-for-ss-regions`, `--bonus`) are now only applied when both `minimap2_align_use_junc_bed` and `minimap2_align_use_splice_scores` are enabled, preventing redundant or conflicting scoring modes.

### Bug fixes

- Corrected the annotation channel format in the `SPLICING` subworkflow to prevent the annotation path from being incorrectly wrapped in a meta tuple when passed to `SPLICEAI_DERIVE`.
- Fixed `BED2GTF` invocation in the preprocessing subworkflow to supply a proper meta map instead of relying on the bed file's base name directly.
- Ensured `global_output_dir` is properly referenced in the Slurm scheduler script (`assets/sh/do_isopipe.sh`).

## [v2.0.21] - 2026-06-18

This release introduces a major architectural shift with the integration of the Ultra aligner alongside the existing minimap2 pipeline, giving users the flexibility to choose between alignment backends. It also brings significant improvements to splicing score annotation, Veredict classification, and HPC job orchestration.

### Ultra aligner integration
- Added a dedicated ultra aligner branch with new subworkflow `ultra_align/main.nf` and nf-core modules for `ultra/index` and `ultra/align`. The pipeline now branches at alignment-time: if `aligner = "minimap2"` the existing minimap2 path is used; if `aligner = "ultra"` the new Ultra path is taken, including Ultra-specific index preparation, alignment, BAM merging, and adapter removal.
- Updated `main.nf`, `nextflow.config`, and preprocessing subworkflows to conditionally dispatch channels based on the selected aligner, and added `ultra_index` and `ultra_use_annotation` parameters.
- Fixed merging step to be Ultra-specific and corrected input Ultra index channel wiring.
- Updated collateral modules to align with the new aligner branching logic.

### Splicing score annotation
- Bumped `splicing-rs` to v0.0.6, adding `--bonus` and `--bonus-score` flags to reward annotation-derived spliced sites during spliceai derive, improving the quality of junction-level scoring for annotated transcripts.
- Propagated the `--bonus` flag through the spliceai derive module.

### Veredict classification
- Veredict.py updated to v0.0.8, implementing a new `ARTIFACT` output category that allows the pipeline to explicitly flag and separate technical artifacts from genuine biological signals in the polishing step.

### xloci improvements
- Added `--unmask` option to xloci intron extraction, providing finer control over repeat masking during intron boundary delineation.

### Infrastructure and HPC
- Added a new HPC matrix scheduler script (`assets/sh/do_isopipe.sh`) to streamline job array submission and resource allocation across Slurm-based clusters.
- Fixed annotation channel tuple propagation across downstream processes to prevent channel misalignment in multi-sample runs.

## [v2.0.20] - 2026-05-28

### Breaking Changes
- v2.0.20 -> conditional cigar extension step!
- v2.0.20 -> splicing v0.0.5, include SS from --regions in final derived score tsv
- fix collapse impl
- static -> add additional collapsing step to remove extremely deep duplicates after segmenting
- V2.0.19 static -> iso-align impl + pooled reads to find fragments
- v2.0.19 -> re-align implementation of fragments

### Features
- SS from regions + SS CDS on default in spliceAi derive

### Fixes
- remove --collapse-mode from collapse mod
- ensure collapses fills up new channel
- avoid grep code 1 on non-duplicates vs duplicates
- update schema + include RT in output

### Chores
- bump minimap2 image to latest version
- pre-release iso-align + fragment logic
- bump submodules

## [v2.0.19] - 2026-04-24

### Breaking Changes
- v2.0.19 -> impl collpase for passes
- v2.0.19 -> fix veredict RT masking
- v2.0.19 -> include collapse in schemas

### Features
- v2.0.18 static -> replace gawk join begraph with bigwig merge and single bg to bw
- static -> add global_prefix to avoid colliding files in scratch

### Chores
- drop ucsc binaries
- bg baseName for bw + bigwigmerge publish

## [v2.0.18] - 2026-04-22

### Breaking Changes
- v2.0.18 -> add bigWig to track + increase classify core count + add process_high_low_memory_fast
- v2.0.18 -> drop validate_id return error, keep as warn

### Features
- add genepred lint in preprocessing

### Chores
- modify process for adapter

## [v2.0.17] - 2026-04-22

### Breaking Changes
- v2.0.17 -> relax validation id allowing asymmetry + warnings
- v2.0.17 -> add isotools adapter

## [v2.0.16] - 2026-04-20

### Breaking Changes
- v2.0.16 -> fix APARENT padding lower limit to 10bp

### Chores
- bump isotools intron process
- rm scrap

## [v2.0.15] - 2026-04-19

### Breaking Changes
- v2.0.15 -> staging options in config + split align merging on file name stripped of cigar extension
- v2.0.15 -> aparent interval logic re-impl
- isotools cigar impl -> duplicate isotools segment [add bai to input channel]
- static v2.0.14 -> add splicing scores for orphan + port logic into preprocess + splicing
- static v2.0.14 -> strip WGET duplicated processes + allow local path avoiding download

### Fixes
- clone meta for extended in isotools-cigar
- update WGET calls for output file

### Chores
- bump submodules
- bump xorf
- copyright + drop extract
- drop minimap correction -> replaced by iso-cigar
- bump isotools
- add FILL marker to xorf_selenocysteine
- Change condition to clean up temporary BAM files

## [v2.0.14] - 2026-04-10

### Breaking Changes
- v2.0.14 -> splicing v0.0.4; automate bigwig naming tokens [plus/minus + donor/acceptor]
- v2.0.14 -> isotools orphan impl + bump updates + loading track statements
- v2.0.14 -> trackDb impl + nf mod
- BREAKING CHANGE: v2.0.14 -> spliceAi impl [chunk, prediction, publish] + side tools bigwigmerge, wigtobigwig + rust/py impl + docker img + sync subworkflow inputs
- v2.0.14 -> spliceai pipeline impl
- v2.0.14 -> spliceai impl

### Features
- add gxf2bed mod
- feat/fix: add detach + join fusions/nmd -> bb + publish

### Fixes
- add procps to spliceai img

### Chores
- bump xorf
- add spliceai img
- fmt + correct output channel names
- fmt + img version

## [v2.0.13] - 2026-04-09

### Notes
- WARN: revert versioning, typo -> v2.0.13

## [v2.0.12] - 2026-04-06

### Breaking Changes
- v2.0.12 -> track upload impl

### Features
- add ssh impl + include rsync/ssh in rs img
- --autoseql + bed12 base schema
- publish joined bigbeds

### Fixes
- remove empty outputs + add sorting step on bed extension

## [v2.0.11] - 2026-04-05

### Breaking Changes
- v2.0.11 veredict nf mod
- v2.0.11 -> isotools + generic join + bedtobigbed impl
- v2.0.11 -> veredict impl v0.0.3 + new bb schema

## [v2.0.10] - 2026-04-04

### Breaking Changes
- v2.0.10 -> v0.0.6 aparent, fix threshold type to float
- v2.0.10 -> merge bed + introns as single input for polishing step

## [v2.0.9] - 2026-04-04

### Breaking Changes
- v2.0.9 -> isotools polish + re-factoring xloci, bigtools, aparent (--threshold)

## [v2.0.8] - 2026-04-02

### Breaking Changes
- v2.0.8 -> aparent v0.0.4, fix ghost outputs

### Fixes
- force single channel for aparent joins
- cover lima multi bam output

### Chores
- upgrade isoseq cluster2 reqs

## [v2.0.7] - 2026-03-26

### Breaking Changes
- v2.0.7 -> bump aparent, fix printing statement

## [v2.0.6] - 2026-03-26

### Breaking Changes
- v2.0.6 -> aparent reverse sorted output + bg join + bigwig conversion impl
- v2.0.6 -> bg/bw + ghost outputs fix + conversion impl

### Chores
- emit chromsizes

## [v2.0.5] - 2026-03-26

### Breaking Changes
- v2.0.5 -> aparent predict impl!
- v2.0.5 -> prepolish placeholders from other subworkflows

### Features
- aparent mod inclusion [chunk + predict]
- wget mod

### Fixes
- include pandas

### Chores
- fix CLI tool name
- intronIC fixed args
- update args in xloci intron
- update entry command
- bump submodule
- bump xorf

## [v2.0.4] - 2026-03-25

### Breaking Changes
- v2.0.4 -> prepolish subworkflow impl
- update img; include aparent, publish py img

### Features
- aparent ini + docker py
- add intronIC mod
- add isotools classify intron mod
- add xloci intron mod
- multi-sample to BAM before merging -> change SAMTOOLS_MERGE input

### Chores
- modify align process specs
- fmt
- samtools merge as process_low_long
- watermark

## [v2.0.3] - 2026-03-16

### Breaking Changes
- v2.0.3 -> splicing v0.0.3; adjusting coords to match minisplice coords, avoids spurious alignment shift

### Fixes
- derive junc-bed as value Channel instead of direct config path + config publishes minimap2 sam/bam

## [v2.0.2] - 2026-03-15

### Breaking Changes
- v2.0.2 -> samtools merge in mapping branch

### Features
- open port for minisplice preloaded scores

### Fixes
- match minimap2 index emitted channel preloaded/raw; explicit emptiness with null splice_scores

### Chores
- fix entrypoint docs + update latest img in spliceai derive

## [v2.0.1] - 2026-03-14

### Breaking Changes
- bump v2.0.1 -> isox-rs/splicing v0.0.2

### Fixes
- correct building path
- update bigwig arg

### Chores
- add docs
- change version to first release of isox-rs

## [v2.0.0] - 2026-03-13

### Breaking Changes
- v2.0.0 -> naming, re-factoring + CI publishing workflow
- v2.0.0 -> splicing scores prediction through spliceai/minisplice + compatibility
- v2.0.0 -> split align/bam process; bumps to latest minimap2 version + --spsc
- v2.0.0 codebase upload
- v2.0.0 init!

### Features
- updated assets
- isolate all params in main.nf

### Fixes
- drop tmps

### Chores
- bump time, bytes, slab dependencies

## [v1.0.0] - 2025-12-20

### Fixes
- missing comma leading to reading column names out of frame
