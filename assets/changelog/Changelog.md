# Changelog

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
