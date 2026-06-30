# `longshot` Codebase Audit Plan

**Date:** 2026-06-30
**Scope:** `main.nf`, all files under `subworkflows/`, `modules/`, `conf/`, and `nextflow.config`. `CLAUDE.md` was read first for terminology/architecture conventions.
**Methodology:** Manual cross-reference of every `process` definition in `modules/*.nf` against every `include { ... }` and pipe invocation in `subworkflows/*.nf` and `main.nf`, against the canonical process list from `bash scripts/test.sh`, and against `withName:` selectors in `conf/process.config`. This is a research/documentation pass only — no `.nf` or `.config` files were modified.

---

## A. Missing Code

- **A1**: `params.tenx_5kit_primers` is required by `LIMA_ISOSEQ`/`ISOSEQ_REFINE`/`CONSTRUCT_MULTIPLEX_PRIMERS` 5′-chemistry branches (`subworkflows/preprocess.nf:79,89`, `modules/preprocess.nf:196`) but `main.nf`'s pre-flight `check_file()` block (`main.nf:65–86`) never validates `params.tenx_5kit_primers` exists, only `tenx_3kit_primers`. A 5′-chemistry run will fail at runtime deep inside `PREPROCESS` instead of failing fast at pre-flight. → Add `check_file(params.tenx_5kit_primers, "10x 5kit primers FASTA")` alongside the existing `tenx_3kit_primers` check in `main.nf:79`.
- **A2**: `modules/quantify.nf` defines `SQANTI3_QC` (line 43) emitting a `junctions` output (`modules/quantify.nf:56`), but the actual sharded/chunked path used in production (`SQANTI3_QC_CHUNK`, `modules/quantify.nf:312-325`) does not emit a `junctions` file at all. If per-experiment junction QC is a wanted output, it is silently missing from the real execution path. → Either add a `junctions` emit to `SQANTI3_QC_CHUNK`/`SQANTI3_MERGE_CHUNKS`, or drop the unused `junctions` emit when removing the dead `SQANTI3_QC` process (see B1).
- **A3**: `docs/plan-ase-pipeline.md` (referenced from `CLAUDE.md:93,107`) describes `subworkflows/ase.nf` and `modules/ase.nf` as part of the "Future Plans" ASE sub-pipeline, but neither file exists under `subworkflows/` or `modules/`. This is expected/by design (explicitly scoped as future work in `CLAUDE.md`), not an oversight — flagged only for completeness. → No action needed.

---

## B. Dead Code

- **B1**: `modules/quantify.nf:1-40` (`ISOQUANT_DISCOVERY`) and `modules/quantify.nf:43-81` (`SQANTI3_QC`) are unsharded/unchunked process definitions that are never invoked by any subworkflow — `subworkflows/classify.nf` only calls `ISOQUANT_DISCOVERY_SHARD` and `SQANTI3_QC_CHUNK`/`SQANTI3_MERGE_CHUNKS`. These are leftover pre-sharding-refactor implementations (~80 lines combined). → Delete both processes from `modules/quantify.nf`.
- **B2**: `conf/process.config:66-67` still defines `withName: 'ISOQUANT_DISCOVERY' { ... }` and `withName: 'SQANTI3_QC' { ... }` resource blocks for the dead processes in B1 — including a 600 GB / 48h allocation that matches nothing in the real DAG. These are orphaned selectors left over from the resource-tier-to-explicit-allocation refactor. → Delete both `withName` blocks once B1 is applied.
- **B3**: `conf/resources.config:40` declares `container_seqkit = "quay.io/biocontainers/seqkit:2.13.0--he881be0_0"`, but no process anywhere in `modules/` references `params.container_seqkit` — `seqkit` is never invoked by any process. → Remove the unused `container_seqkit` param.
- **B4**: `nextflow.config:16-17`'s `params.polya_peak` / `params.intropolis` and their optional-arg handling in `modules/quantify.nf:60-61` are wired into the dead `SQANTI3_QC` (B1) as well as the live `SQANTI3_QC_CHUNK` (`modules/quantify.nf:328-329`) — once B1 is deleted, the duplicate wiring in `SQANTI3_QC` goes away automatically; only the chunk-path wiring remains (and is legitimate, in active use). → No separate action beyond completing B1.

---

## C. Overcomplicated Code

- **C1**: `subworkflows/quantify.nf` (21 lines total) is a single-process passthrough: it maps/combines a channel and calls `ISOQUANT_QUANTIFY` with no additional branching, QC, or multi-step logic — effectively a thin wrapper subworkflow around one process. Per `CLAUDE.md:82` ("group multi-step logic in subworkflows"), a single-process subworkflow with no real orchestration is an unnecessary abstraction layer. → Fold `ISOQUANT_QUANTIFY` directly into `main.nf` as a top-level process call, removing `subworkflows/quantify.nf` entirely (~21 lines + 1 file) — unless it's deliberately kept separate as a future branch point for the ASE sub-pipeline, in which case document that rationale inline.
- **C2**: `modules/quantify.nf:254-309` (`SQANTI3_SPLIT_GTF`) and `modules/quantify.nf:352-389` (`SQANTI3_MERGE_CHUNKS`) implement ~135 lines of bespoke Python GTF-block-parsing logic purely to parallelize SQANTI3 QC into `params.sqanti_chunks` (default 8) pieces. This is a second, independent sharding mechanism layered on top of the chromosome-based IsoQuant sharding already computed earlier in the same subworkflow (`subworkflows/classify.nf:11-31`, `isoquant_shards()`). → Reuse the existing chromosome shard boundaries for SQANTI3 chunking instead of a second transcript-count-based custom splitter/merger pair, eliminating one bespoke Python implementation.
- **C3**: `subworkflows/classify.nf:11-31` (`isoquant_shards()`) hardcodes a 17-shard chromosome partition as a Groovy function embedded directly in the subworkflow file, mixing pipeline orchestration with static configuration data. This makes shard topology invisible to anyone scanning `nextflow.config`/`conf/` for tunables, unlike `params.sqanti_chunks` which is a proper config-level parameter. → Move `isoquant_shards()` into config as a `params.isoquant_shards` list, consistent with how `sqanti_chunks` is already exposed.
- **C4**: `main.nf:90-140` embeds ~50 lines of inline Groovy (chemistry auto-detection, adapter-kit shorthand resolution, null/NA/none string normalization) directly inside the channel `.map{}` closure in the main workflow block. The "is this value null/NA/none" check is duplicated nearly verbatim three times (`main.nf:97-98` for `shortread_barcodes`, `main.nf:100-101`/`106` for `tenx_index`/`chemistry`, `main.nf:117-118` for `adapter_kit`). → Extract a single `is_null_value(String v)` helper next to `check_file`/`preflight_samplesheet` (`main.nf:17-58`) and call it from all sites, removing ~10 duplicated lines.
- **C5**: `main.nf:83-86` unconditionally validates that `mas8_primers.fasta`, `mas12_primers.fasta`, and `mas16_primers.fasta` all exist under `assets/adapters/` on every pipeline run, regardless of whether any samplesheet row actually uses those shorthand kits — while custom `adapter_kit` paths supplied directly in the samplesheet are never validated at all (the resolved `adapter_file` in `main.nf:121-123` is passed downstream with no existence check). This is inconsistent defensive coding: it over-validates fixed shorthand files unconditionally but under-validates the actual per-row values that vary at runtime. → Only validate shorthand kit files for kits actually referenced by samplesheet rows, and add an existence check for resolved custom `adapter_kit` paths to make validation coverage consistent.

---

## D. Naming

- **D1**: The `CRAMINO` family (`CRAMINO_RAW`, `CRAMINO_SEGMENTED`, `CRAMINO_FILTERED`, `CRAMINO_FL`, `CRAMINO_FLTNC` aliased in `subworkflows/preprocess.nf:2`, plus a bare `CRAMINO` in `subworkflows/align.nf:8`) are all `CRAMINO as X` aliases of the single underlying `process CRAMINO` in `modules/qc.nf:56`, which already takes a `stage` string argument (`modules/qc.nf:62`) to label its output. The alias-per-stage pattern is justified — it lets `conf/process.config:35-39` give each pipeline stage its own resource tier, and it matches the documented "read-loss funnel" QC design (`CLAUDE.md:73`) — so this is **not unjustified duplication**. However, the naming is inconsistent: 5 of 6 invocations use a `CRAMINO_<STAGE>` alias while the alignment-stage one is just bare `CRAMINO` (`subworkflows/align.nf:8,77`), despite being invoked 3 times with different `stage` values (`merged`/`dedup`/`aligned`, `subworkflows/align.nf:71-77`) in the same subworkflow. → For consistency, alias it as `CRAMINO as CRAMINO_ALIGN` in `subworkflows/align.nf:8` (matching the `PREPROCESS:CRAMINO_*` convention), or standardize on bare process names plus the `stage` tuple argument everywhere — pick one convention repo-wide.
- **D2**: No redundant subworkflow-name-in-process-name violations were found: processes inside `PREPROCESS`, `ALIGN`, `CLASSIFY`, `QUANTIFY`, and `EXPORT` are all named for what they do (e.g. `SKERA_SPLIT`, `ISOSEQ_CORRECT`) without repeating their parent subworkflow's name as a prefix. No action needed.
- **D3**: `modules/exporter.nf` is named after the *subworkflow* it backs (`EXPORT`) using an agent-noun form (`-er`), while the sibling files `modules/align.nf`, `modules/preprocess.nf`, `modules/quantify.nf` all use the bare stage noun. → Rename `modules/exporter.nf` → `modules/export.nf` for consistency with the other three module files (cosmetic, low priority).
- **D4**: `CLAUDE.md:3` still parenthetically notes the project was "formerly `longreadr`", and the directory-layout diagram at `CLAUDE.md:58` still shows the tree rooted at `longreadr/` instead of `longshot/`. This is a documentation-only naming leftover from the project rename, not a code issue, but could mislead a new contributor about the actual repo directory name. → Update the tree root in `CLAUDE.md:58` to `longshot/` in a future docs cleanup (out of scope for this code-only audit).

---

## Estimated Impact

Applying all of Section B (dead code) and Section C (overcomplicated code) fixes:

- **B1+B2**: Remove `ISOQUANT_DISCOVERY` + `SQANTI3_QC` processes (~80 lines, `modules/quantify.nf`) and their two `withName` blocks (~4 lines, `conf/process.config`) — **~84 lines removed**, 0 files removed.
- **B3**: Remove unused `container_seqkit` param (`conf/resources.config:40`) — **1 line removed**.
- **C1**: Fold `subworkflows/quantify.nf` into `main.nf` — **~21 lines removed, 1 file removed**.
- **C2**: Reuse chromosome shards for SQANTI3 chunking instead of a second custom GTF splitter/merger — **~100 lines removed** if implemented (largest single simplification, but higher design/testing effort than the other items, so treat as a candidate rather than a quick win).
- **C3**: Move `isoquant_shards()` into config — roughly line-count-neutral (moved, not deleted); improves discoverability rather than reducing code, not counted in the total below.
- **C4**: Deduplicate the null/NA/none check in `main.nf` — **~10 lines removed**.

**Total estimated reduction if all of B + C (including the higher-effort C2) are applied: ~215 lines removed across 4 files (`modules/quantify.nf`, `conf/process.config`, `conf/resources.config`, `main.nf`), plus 1 file removed (`subworkflows/quantify.nf`).**

A conservative "quick wins only" pass (B1–B3 + C1 + C4, deferring the riskier C2 GTF-splitter redesign) yields **~115 lines removed across 4 files, 1 file removed**.
