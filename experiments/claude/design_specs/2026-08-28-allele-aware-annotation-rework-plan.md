# Allele-Aware Annotation Rework — Implementation Plan

> **For agentic workers:** execute task-by-task with superpowers:subagent-driven-development. Steps use `- [ ]` checkboxes.

**Goal:** Replace the allele-less annotation source (`.annot.gz` + arbitrary `unique(subset=SNP)`) with the allele-bearing `.annot.parquet` (has `A1/A2`), so every variant — including the ~0.2%/chr multiallelic sites whose annotations genuinely differ — carries its own correct annotations. Make the ridge fit and the explainability join allele-aware. Add a dedicated hg19 genetic-map asset for the plot's recombination track (the allele-bearing parquet has no `CM`).

**Why:** Measured on real chr15 data — `(CHR,BP)` does not uniquely identify a variant (~1,242 multiallelic rsIDs/chr, annotations differ up to 49/189 columns); `snpvar_meta` has 18,859 multi-allele `[CHR,BP,SNP]` groups. The current SNP-dedup silently discards one allele's annotations, and the previously-added `_assert_annotation_positions_unique` guard would break real runs. Evidence: `experiments/claude/polyfun_explain_probe/` + scratchpad `check_*.py`/`investigate.py`.

**Design decisions (settled with the user):**
- Annotation source → allele-bearing `.annot.parquet` from the 30GB `baselineLF_v2.2.UKB.polyfun.tar.gz`; no dedup (already unique on `[CHR,BP,A1,A2]`).
- Ridge annotation↔snpvar_meta join → allele-aware (unordered allele key); re-fit.
- Explainability join (contrast + plot) → allele-aware on `[CHR,POS,unordered-allele-key]`, replacing the guard.
- CM → new dedicated **hg19** genetic-map reference asset (polyfun panel is hg19/build 37).

## Global Constraints

- Everything via pixi (`pixi r ...`, `pixi r invoke green` after each task; capture to a logfile and check the exit code — testmon may report "no tests ran").
- Prefer polars; Path/PurePath; Literal for enumerable strings; column names from constants. Docstrings: no backticks/RST.
- Tests Task-level, dependency-injected (no monkeypatch/skipif); no assertions on error-message text; share creation across tests.
- Reuse the existing `unordered_allele_key(ea, nea)` helper (in `susie_r_finemap_task.py`, used by `align_data` for the prior join) for every allele-set match, so annotation↔snpvar↔gwas all match the same way regardless of which allele is labelled "effect".
- Column-name constants for alleles: use `GWASLAB_EFFECT_ALLELE_COL="EA"`, `GWASLAB_NON_EFFECT_ALLELE_COL="NEA"` on the gwas side; the annotation/snpvar side use `A1`/`A2` — introduce `ANNOT_A1_COL="A1"`, `ANNOT_A2_COL="A2"` constants (in `build_baseline_lf_annotation_parquet_task.py`) rather than literals.

## Ordering & dependency

A (genetic map) is independent. B (annotation source) → C (ridge refit) → D (contrast join) → E (plot join + CM). Then the previously-written Task 4 generators / Task 5 demonstrator get the new gene-map task wired in. A can run in parallel with B/C.

---

## Task A: hg19 genetic-map reference asset + recomb helper

**Files:**
- Create: `mecfs_bio/assets/reference_data/genetic_map/genetic_map_hg19.py` (Download + parse task instances).
- Create: `mecfs_bio/build_system/task/genetic_map/parse_genetic_map_task.py` — `ParseHg19GeneticMapTask`.
- Create: tests under `test_mecfs_bio/unit/build_system/task/genetic_map/`.

**Interfaces:**
- Produces `GENETIC_MAP_HG19` (a Task whose output parquet has columns `CHR` (Int), `POS` (Int, hg19/build37 bp), `RECOMB_RATE_CM_PER_MB` (Float), `GENETIC_MAP_CM` (Float)), and module constants for those column names (`GMAP_CHR_COL`, `GMAP_POS_COL`, `GMAP_RATE_COL`, `GMAP_CM_COL`).
- Source: the Eagle hg19 genetic map `genetic_map_hg19_withX.txt.gz` (whitespace-separated, columns `chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)`), a few MB. Pin `md5_hash=None` first, then pin after first download (mirror `BASELINE_LF_ANNOTATION_TARBALL`). If that exact URL is unavailable at build time, any standard hg19/GRCh37 genetic map with per-position cM works — keep the parse task's output schema fixed and adapt only the parser.

**Steps (TDD):**
- [ ] Write a failing task test: feed a tiny synthetic map file (a handful of rows across 2 chromosomes) via `FakeTask`/injected `fetch`, run `ParseHg19GeneticMapTask`, assert output columns/dtypes and that `RECOMB_RATE_CM_PER_MB` is carried through (or derived via `np.gradient(cM, pos)*1e6` if the source lacks a rate column — pick based on the real columns).
- [ ] Implement the download instance + parse task (stream the gz, `scan_csv` whitespace, rename to the fixed schema, sort by `[CHR,POS]`, write parquet). Use `execute_command`/repo download task patterns; do not read from `~`.
- [ ] `invoke green` + commit.

Note: the plot currently derives `dCM/dBP` from the annotation `CM` column; this asset replaces that source. Prefer the source's own `COMBINED_rate(cM/Mb)` column if present (already a rate) over re-differentiating.

---

## Task B: annotation source → allele-bearing `.annot.parquet`

**Files:**
- Modify: `mecfs_bio/assets/reference_data/polyfun/annotations/baseline_lf_annotations.py` (URL + docstring).
- Modify: `mecfs_bio/build_system/task/annotation_weights/build_baseline_lf_annotation_parquet_task.py`.
- Modify: `test_mecfs_bio/unit/build_system/task/annotation_weights/test_build_baseline_lf_annotation_parquet_task.py` (or create if absent).

**Changes:**
- `BASELINE_LF_ANNOTATION_TARBALL.url` → `https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.polyfun.tar.gz` (the 30GB tarball; verified via the probe to contain per-chromosome `baselineLF2.2.UKB.<chr>.annot.parquet` members with `A1/A2`). Reset `md5_hash=None` (re-pin after first real download — the hash changes with the new tarball).
- `BuildBaselineLFAnnotationParquetTask`:
  - Member regex → `baselineLF2\.2\.UKB\.(\d+)\.annot\.parquet$`; extract those members.
  - `_scan_one_chromosome` → `pl.scan_parquet(member_path)` (not `scan_csv`). The parquet already carries `A1/A2` (String) and 187 Float annotation columns; there is **no `CM`** column.
  - `ANNOT_KEY_COLUMNS` → `["CHR", "BP", "SNP", "A1", "A2"]` (drop `CM`, add alleles). Cast only the true annotation columns (everything not in `ANNOT_KEY_COLUMNS`) to Float32; keep `A1/A2` as String.
  - **Remove** `unique(subset=_SNP_COL)`. Instead, add an unordered allele key `ak = unordered_allele_key(A1, A2)` and `unique(subset=[CHR, BP, ak], keep="first")` to collapse only the lossless ordering-duplicates (~8/chr; their annotations are identical — verified). Do **not** dedup by SNP.
  - Sort by `[CHR, BP]` as today; write one parquet.
  - Keep the `assert members` completeness check; additionally assert the extracted set covers chromosomes 1..22 (fail fast if a member is missing from the 30GB stream).
- Update the module docstrings (the "no A1/A2" and "191" notes are now false).

**Steps (TDD):**
- [ ] Failing task test: synthetic 2-chromosome input mimicking `.annot.parquet` (CHR/BP/SNP/A1/A2 + a couple real annotation columns), including a multiallelic site (same rsID, two allele pairs, different annotation values) and an ordering-duplicate (A1/A2 swapped, identical annotations). Assert: output retains `A1/A2`; the multiallelic site keeps **both** rows; the ordering-duplicate collapses to one; output unique on `[CHR,BP,ak]`.
- [ ] Implement; run the test.
- [ ] `invoke green` + commit. (The real 30GB build is a separate manual/CI step; unit test uses synthetic input.)

---

## Task C: ridge weights — allele-aware join + re-fit

**Files:**
- Modify: `mecfs_bio/build_system/task/annotation_weights/ridge_annotation_weights_task.py`.
- Modify: `test_mecfs_bio/unit/build_system/task/annotation_weights/test_ridge_annotation_weights_task.py`.

**Changes:**
- `snpvar_meta` load: keep `A1/A2`; **remove** `unique(subset=_SNP_COL)`.
- Replace `_JOIN_KEYS = [CHR, BP, SNP]` with an allele-aware join: add `ak = unordered_allele_key(A1,A2)` on both the annotation slice and the snpvar_meta slice, join on `[CHR, BP, ak]` (`how="inner"`). This pairs each annotation allele with its own `snpvar_bin`. Update the stale comment (lines 59-64) to describe the allele-aware join.
- `_annotation_columns` must now also exclude `A1/A2` (they are in `ANNOT_KEY_COLUMNS` after Task B, so `[c for c in schema if c not in ANNOT_KEY_COLUMNS]` already drops them — verify).
- Assert the annotation↔snpvar join does not multiply rows beyond the annotation slice (each annotation `[CHR,BP,ak]` matches ≤1 snpvar row) — fail fast, mirroring Task B's uniqueness posture.

**Steps (TDD):**
- [ ] Update the existing ridge test: its synthetic annotation + `snpvar_meta` fixtures currently have no alleles. Add `A1/A2` to both, include one multiallelic position with two alleles carrying different `snpvar_bin`, and assert the fit pairs each allele's annotations with its own target (e.g. the recovered weights still match a known linear truth built allele-aware). Keep the closed-form/`mean_heldout_r2` assertions meaningful.
- [ ] Implement; run the test.
- [ ] `invoke green` + commit.

Note: this re-fit changes `BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS` values; downstream contrast gammas come from here, so B+C must land before D is verified on real data.

---

## Task D: rework the contrast join — allele-aware, remove the guard

**Files:**
- Modify: `mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_contrast_task.py`.
- Modify: `test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_contrast_task.py`.

**Changes:**
- `_load_annotations` now returns an allele-bearing slice (`CHR, POS(=BP), A1, A2, <annots>`). Remove `_assert_annotation_positions_unique` (the `(CHR,POS)`-uniqueness invariant is false). Instead:
  - add `ak = unordered_allele_key(A1, A2)` to the annotation slice; `unique(subset=[CHR,POS,ak])` (lossless ordering-dedup) then assert `[CHR,POS,ak]` unique.
- The run-variant↔annotation joins (`uni_annot`, `pf_annot`) change from the CHR/POS-only `_ANNOT_KEY` to `[CHR, POS, ak]`, with `ak = unordered_allele_key(EA, NEA)` on the run side. Assert each run variant matches ≤1 annotation row.
- Keep everything downstream (contrast, family, display, selection) unchanged — the join now attaches each run variant's own annotations.

**Steps (TDD):**
- [ ] Update the fixture: annotation gains `A1/A2`; replace `test_contrast_raises_on_duplicate_annotation_position` with `test_contrast_resolves_multiallelic_by_allele` — two annotation rows at one position with different alleles and different annotation values; two run variants at that position with the matching alleles; assert each run variant receives its own (distinct) annotation-derived contrast, and the task does NOT raise.
- [ ] Keep the closed-form test (add trivial `A1/A2` to its single-allele fixture so the allele join is exercised).
- [ ] Implement; run the tests.
- [ ] `invoke green` + commit.

---

## Task E: rework the plot — allele-aware join + CM from the genetic-map asset

**Files:**
- Modify: `mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_plot_task.py`.
- Modify: its test.

**Changes:**
- Annotation join for the family-scaled panels → allele-aware (reuse the contrast module's now-allele-aware `_load_annotations` + `[CHR,POS,ak]` join), so panels still match the tables.
- Recombination track: add a `genetic_map_task` dependency (`GENETIC_MAP_HG19` from Task A). Read its parquet, slice to the locus window, plot `RECOMB_RATE_CM_PER_MB` vs `POS` on the Manhattan secondary axis (replaces the annotation-`CM` `np.gradient` path — the annotation no longer has `CM`). Update `create()`/`deps` to include `genetic_map_task`.
- Remove the old `_load_cm`/annotation-CM recomb path.

**Steps (TDD):**
- [ ] Update the smoke test: add a synthetic genetic-map `FakeTask` + fetch entry; pass `genetic_map_task`; assert both PNG and SVG still written.
- [ ] Implement; run the test.
- [ ] `invoke green` + commit.

---

## Task F: rewire generators + demonstrator for the gene-map dep

**Files:**
- Modify: `mecfs_bio/asset_generator/polyfun_explain_fine_mapping_asset_generator.py` (Task 4 output) — add `GENETIC_MAP_HG19` into `generate_polyfun_explain_group`'s `PolyfunExplainPlotTask.create(...)` call; thread through `SharedFineMapInputs` if per-locus, else reference the module-level asset directly.
- Modify: its test (the 8-run assertions are unchanged; just construction must still succeed with the new dep).
- Modify: the Task 5 demonstrator module if needed (no change if the generator references `GENETIC_MAP_HG19` directly).

**Steps:**
- [ ] Update generator + test; run; `invoke green` + commit.

---

## Manual verification (post-merge, documented — not a unit test)

Requires building: the 30GB annotation tarball → allele-bearing parquet (Task B, assert 22 chromosomes), the re-fit ridge weights (Task C), the hg19 genetic map (Task A). Then build one demonstrator plot target and confirm: no allele-join assertion fires on real data (or, if it does, that is the signal to bring to the user); family panels match the display table; recomb track renders from the genetic map; PNG+SVG written.

## Out of scope / risks

- 30GB download and ridge re-fit are one-time build costs (annotation build streams the tarball once, keeps `.annot.parquet` members ~1.5GB). `.annot.parquet` presence for all 22 chromosomes is asserted at build time (cannot be cheaply pre-listed — non-seekable gzip).
- The old `.annot.gz`-derived `CM` is abandoned; the recomb track now comes from the dedicated hg19 map.
- Superseded: the retracted `snpvar_meta`-bridge idea and Task 2's `_assert_annotation_positions_unique` guard (removed in Task D).
