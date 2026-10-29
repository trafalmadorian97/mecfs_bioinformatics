# Large-scale Trait × Proteomics Colocalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** For each of four trait GWAS, colocalize (coloc.abf) its genome-wide significant loci with
the UKB-PPP proteins that have a pQTL in the same locus, within a selection budget, and emit a raw
results table, a docs-ready display table, and stacked locus plots.

**Architecture:** A reusable per-source *protein hit table* (HapMap3 slim files, configurable p
threshold) drives candidate selection. One Task per trait selects (locus, protein) pairs by score
within a pair budget, fetches dense protein data on demand (download, slice windows, discard),
orients each window with the existing untrusted genome-reference harmonizer, aligns exactly to the
trait, and runs coloc.abf behind a method interface. Two small downstream Tasks produce the display
table and the plots from the per-trait result directory.

**Tech Stack:** Python 3.13, polars, numpy/scipy, attrs (frozen, slots), matplotlib OO API, the
repo's Task build system; pixi for everything (`pixi r ...`).

**Spec:** `experiments/claude/design_specs/2026-10-29-large-scale-protein-coloc-design.md`

## Global Constraints

- Analysis build is GRCh38 throughout; every join is in GRCh38; assert build from metadata.
- Trait input must be a `FilteredGWASDataMeta` with `harmonization_info.build == "38"`; assert in `.create()`.
- Protein windows are oriented with `resolve_chromosome(..., ChromosomeContext(trusted=False, ...), ParquetPanelLoader(...))`; panel `GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES`, ancestry `"nfe"`.
- Trait ↔ protein matching is an exact join on (CHR, POS, EA, NEA) after both sides are "NEA = hg38 REF". Never use `unordered_allele_key`.
- PPP v1 cohort: EUR discovery (`EUR_DISCOVERY_PPP_MANIFEST_PATH`, `HAPMAP_3_PPP_DATABASE`).
- Defaults (all parameters): hit-table p < 5e-8; trait lead p < 5e-8; trait window ±500 kb; cis = gene within ±1 Mb of window; pair budget 30; per-locus cap 5; −log10 p clipped at 50; MHC (`extended_mhc_interval("38")`) excluded; coloc priors p1 = 1e-4, p2 = 1e-4, p12 = 1e-5; display PP.H4 > 0.8.
- Memory bounded by one window / one chromosome file (16 GB runner).
- No edits under `docs/`. Docs-facing parquet: Int32 integers, zstd, byte-stream-split floats (`write_parquet_table`).
- Repo conventions: polars; `@frozen(slots=True)`; `.create()` derives meta from deps; Literal aliases for enumerable strings; no `Any`/`Callable[...]`; no mocking (inject fakes); imports at top; column-name constants, not literals; docstrings without backticks or RST; `pixi r invoke green` after each task, output captured to a log file.

## Review Focus

1. **Mirrored indel collision** — a protein deletion TCA→T stored flipped as (NEA=T, EA=TCA) next to a trait insertion T→TCA must not join as the same variant. Test in Task 3.
2. **A protein that is a candidate at two loci on different chromosomes** must be downloaded once and both pairs tested. Test in Task 6.
3. **A selected pair with too few aligned variants** (or none) must get status "failed" with a reason, not crash the run. Test in Task 6.
4. **Duplicate (CHR, POS, EA, NEA) keys** on either side of a window must be dropped before coloc (coloc assumes one row per variant). Test in Task 3.
5. **No pair passes the display threshold** → an empty display parquet and an empty plot index, not an error. Test in Task 7.

---

## File Structure

```
mecfs_bio/constants/protein_coloc_constants.py          # output column names, status Literal + values
mecfs_bio/build_system/task/protein_coloc/
    __init__.py
    coloc_abf.py               # Task 1: trait models, priors, ColocMethod protocol, ColocAbf
    ppp_window_fetcher.py      # Task 2: LocusWindow, DenseWindowFetcher protocol, PppWindowFetcher
    orient_and_align.py        # Task 3: orient_window, align_windows
    protein_hit_table_task.py  # Task 4: PppProteinHitTableTask
    trait_loci.py              # Task 5: find_trait_loci
    candidate_selection.py     # Task 5: build_candidates, select_candidates, SelectionConfig
    trait_protein_coloc_task.py# Task 6: TraitProteinColocTask, ColocRunConfig
    coloc_display_table_task.py# Task 7
    coloc_locus_plot_task.py   # Task 7
mecfs_bio/assets/reference_data/magma_gene_locations/processed/magma_entrez_gene_locations_data_build_38_unzipped.py  # Task 7
mecfs_bio/assets/gwas/ukbb_ppp/ppp_database/hapmap3/eur_discovery_hapmap3_ppp_hit_table.py                         # Task 4
mecfs_bio/assets/gwas/<trait>/.../coloc/   (four trait-prep + coloc modules)                                        # Tasks 8, 9
mecfs_bio/analysis/ppp_coloc_analysis.py                                                                           # Task 9
test_mecfs_bio/unit/build_system/task/protein_coloc/test_*.py
```

---

### Task 1: coloc.abf core and method interface

**Files:**
- Create: `mecfs_bio/build_system/task/protein_coloc/__init__.py`, `mecfs_bio/build_system/task/protein_coloc/coloc_abf.py`
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_coloc_abf.py`

**Interfaces:**
- Produces:
  - `@frozen class QuantitativeTraitModel: sd_y: float | None` (None → estimate sdY from data)
  - `@frozen class CaseControlTraitModel: case_fraction: float` (assert 0 < s < 1)
  - `TraitModel = QuantitativeTraitModel | CaseControlTraitModel`
  - `@frozen class ColocPriors: p1: float = 1e-4; p2: float = 1e-4; p12: float = 1e-5`
  - `@frozen class ColocInput: beta1, se1, beta2, se2: np.ndarray; maf1, n1: np.ndarray | None; maf2, n2: np.ndarray | None; model1: TraitModel; model2: TraitModel` — post_init asserts equal lengths, 1-D, float, se > 0, and that maf/n are present when a QuantitativeTraitModel has sd_y None.
  - `@frozen class ColocPosteriors: pp_h0..pp_h4: float; n_variants: int` — post_init asserts sum ≈ 1.
  - `class ColocMethod(Protocol): name: str; def run(self, data: ColocInput) -> ColocPosteriors`
  - `@frozen class ColocAbf(ColocMethod): priors: ColocPriors = ColocPriors(); name: str = "coloc_abf"`

- [ ] **Step 1: Write failing tests**

```python
def test_shared_causal_variant_gives_high_h4():
    # 200 variants, both traits have z≈10 at index 50, z≈0 elsewhere
    post = ColocAbf().run(_input(z1_at={50: 10.0}, z2_at={50: 10.0}))
    assert post.pp_h4 > 0.9

def test_distinct_causal_variants_give_high_h3():
    post = ColocAbf().run(_input(z1_at={50: 10.0}, z2_at={150: 10.0}))
    assert post.pp_h3 > 0.9

def test_no_signal_gives_high_h0():
    assert ColocAbf().run(_input(z1_at={}, z2_at={})).pp_h0 > 0.9

def test_posteriors_sum_to_one():
    p = ColocAbf().run(_input(z1_at={50: 6.0}, z2_at={52: 4.0}))
    assert abs(p.pp_h0 + p.pp_h1 + p.pp_h2 + p.pp_h3 + p.pp_h4 - 1) < 1e-12

def test_larger_p12_increases_h4():
    data = _input(z1_at={50: 5.0}, z2_at={50: 5.0})
    assert ColocAbf(ColocPriors(p12=1e-4)).run(data).pp_h4 > ColocAbf(ColocPriors(p12=1e-6)).run(data).pp_h4

def test_single_variant_matches_hand_computed_abf():
    # one variant, case-control trait 1 (W=0.2) and quantitative trait 2 with sd_y=1 (W=0.15):
    # assert pp_h4 equals the closed form from the formulas in Step 3, to 1e-12
    ...

def test_sd_y_estimated_when_absent():
    # quantitative model with sd_y=None and maf/n supplied; data simulated with sd_y=2
    # assert the estimate used is within 5% of 2 (expose estimate_sd_y as a module function)
```

`_input` is a test helper building `ColocInput` with se = 0.05, beta = z·se, quantitative sd_y = 1.

- [ ] **Step 2: Run and verify they fail** — `pixi r pytest test_mecfs_bio/unit/build_system/task/protein_coloc/test_coloc_abf.py -v` → import errors.

- [ ] **Step 3: Implement `coloc_abf.py`** following coloc (Giambartolomei 2014; coloc R package `approx.bf.estimates`, `combine.abf`, `sdY.est`). Pinned math (put it in the module docstring per repo convention):
  - prior SD W: 0.15·sdY (quantitative), 0.2 (case-control); V = se²; r = W²/(W²+V); z = β/se;
    lABF = ½(log(1−r) + r·z²).
  - L1 = logsumexp(lABF1), L2 = logsumexp(lABF2), L12 = logsumexp(lABF1 + lABF2).
  - log H0 = 0; H1 = log p1 + L1; H2 = log p2 + L2; H3 = log p1 + log p2 + logdiff(L1 + L2, L12); H4 = log p12 + L12; normalize with logsumexp.
  - `estimate_sd_y(se, maf, n) -> float`: least squares without intercept of 2·n·maf·(1−maf) on 1/se²; sdY = sqrt(coefficient).

- [ ] **Step 4: Run tests** → PASS.
- [ ] **Step 5: `pixi r invoke green > <scratch>/green.log 2>&1`, check exit 0 and "passed" in log; commit** `git commit -m "Add coloc.abf core behind a ColocMethod interface"`.

---

### Task 2: PPP dense window fetcher

**Files:**
- Create: `mecfs_bio/build_system/task/protein_coloc/ppp_window_fetcher.py`, `mecfs_bio/constants/protein_coloc_constants.py`
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_ppp_window_fetcher.py`

**Interfaces:**
- Consumes: `PppProteinFile`, `extract_protein_tar` (`build_system/task/ppp_database/build_slim_protein_parquet_task.py`); `WF.download_from_synapse(synid, dest_dir, filename) -> Path`.
- Produces:
  - `@frozen class LocusWindow: locus_id: str; chrom: int; start: int; end: int` (post_init: start ≤ end, 1 ≤ chrom ≤ 23)
  - `class DenseWindowFetcher(Protocol): def fetch(self, protein_id: str, windows: Sequence[LocusWindow], scratch_dir: Path, wf: WF) -> dict[str, pl.DataFrame]` — keyed by locus_id; frames in the common schema CHR (Int64), POS (Int64, GRCh38), EA, NEA (String), BETA, SE, EAF, N (Float64), *not yet oriented*.
  - `@frozen class PppWindowFetcher(DenseWindowFetcher): proteins: tuple[PppProteinFile, ...]` (tuple keeps the frozen Task hashable; index by OID inside `fetch`) — downloads the tar once per call, reads only the chromosome files the windows need (columns CHROM, GENPOS, ALLELE0, ALLELE1, A1FREQ, N, BETA, SE; EA=ALLELE1, NEA=ALLELE0, EAF=A1FREQ; chrom code X→23), filters each window, then deletes the tar and extracted directory before returning.
  - Column constants in `protein_coloc_constants.py` reuse gwaslab names (`GWASLAB_*`) for the common schema; add `LOCUS_ID_COL = "locus_id"`, `PROTEIN_ID_COL = "protein_id"`.

- [ ] **Step 1: Write failing tests** using a synthetic tar (two chromosome gz files in the real regenie layout, as `test_build_slim_protein_parquet_task.py` writes them) and a WF built with a fake `SynapseDownloader` that copies the tar and counts calls (no mocking; inject via WF's composition, see `build_system/wf/synapse_downloader.py`).

```python
def test_fetch_returns_only_window_rows_in_common_schema(tmp_path): ...
    # windows on chr1 and chr2; assert each frame's POS ⊂ [start, end], columns == COMMON_SCHEMA
    # and EA == ALLELE1, NEA == ALLELE0, EAF == A1FREQ for a known row
def test_fetch_downloads_once_for_windows_on_two_chromosomes(tmp_path):
    assert fake_downloader.n_calls == 1
def test_fetch_deletes_download_and_extraction(tmp_path):
    assert not any((tmp_path / "scratch").rglob("*.tar")) and not any((tmp_path / "scratch").rglob("*.gz"))
def test_window_with_no_rows_returns_empty_frame_with_schema(tmp_path): ...
```

- [ ] **Step 2: Run, verify fail.** **Step 3: Implement.** **Step 4: Run, verify pass.**
- [ ] **Step 5: green (logged); commit** `"Add PPP dense window fetcher"`.

---

### Task 3: Orientation and alignment of windows

**Files:**
- Create: `mecfs_bio/build_system/task/protein_coloc/orient_and_align.py`
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_orient_and_align.py`

**Interfaces:**
- Consumes: `resolve_chromosome`, `ChromosomeContext`, `DROP_REASON_COL` (`genome_reference_harmonization/resolve_chromosome.py`); `ParquetPanelLoader` (`genome_reference_harmonization_task.py`); `resolve_column_rules`; `GenomeReferenceHarmonizationOptions`; `IndexedFasta`.
- Produces:
  - `@frozen class OrientedWindow: frame: pl.DataFrame; n_input: int; drop_counts: Mapping[str, int]` (frame = kept rows only, common schema, NEA = hg38 REF)
  - `orient_window(rows: pl.DataFrame, chrom: int, fasta: IndexedFasta, panel_path: Path, panel_af_col: str, options: GenomeReferenceHarmonizationOptions) -> OrientedWindow` — always `trusted=False`.
  - `@frozen class AlignedWindow: frame: pl.DataFrame; n_trait: int; n_protein: int; n_duplicate_keys_dropped: int` — frame columns CHR, POS, EA, NEA, BETA_TRAIT, SE_TRAIT, EAF_TRAIT, N_TRAIT, BETA_PROTEIN, SE_PROTEIN, EAF_PROTEIN, N_PROTEIN (constants in `protein_coloc_constants.py`).
  - `align_windows(trait: pl.DataFrame, protein: pl.DataFrame) -> AlignedWindow` — drop every key duplicated within either side, then inner join on (CHR, POS, EA, NEA).
  - Strand-flip (reverse-complement) classes are left as resolve_chromosome handles them; the spec leaves measuring or dropping them to the builder. Do not add special handling unless Task 9's validation shows a problem.

- [ ] **Step 1: Write failing tests** with a tiny FASTA and panel built the way the existing `test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/` tests build theirs (reuse their helpers).

```python
def test_swapped_snv_is_flipped():            # EA==REF → after orient NEA==REF, BETA negated, EAF = 1-EAF
def test_unorientable_snv_is_dropped_and_counted():
def test_flipped_deletion_does_not_join_trait_insertion():
    # Genome reads "TCA" from pos p, so both T>TCA (insertion) and TCA>T (deletion) are on-reference.
    # Trait (already REF-oriented) carries the insertion: NEA=T, EA=TCA.
    # Panel has only the deletion: REF=TCA, ALT=T, AF=0.1.
    # Protein carries the deletion stored flipped: NEA=T, EA=TCA, EAF=0.9.
    # orient_window must resolve the protein row to NEA=TCA, EA=T (or drop it);
    # assert align_windows(trait, oriented.frame).frame has no row at pos p
def test_duplicate_keys_are_dropped_on_both_sides():
    assert aligned.n_duplicate_keys_dropped == 2 and aligned.frame.height == expected
def test_alignment_is_exact_on_ordered_alleles():
```

- [ ] **Step 2: fail. Step 3: implement. Step 4: pass.**
- [ ] **Step 5: green (logged); commit** `"Add on-the-fly untrusted orientation and exact alignment of coloc windows"`.

---

### Task 4: Protein hit table Task (configurable threshold)

**Files:**
- Create: `mecfs_bio/build_system/task/protein_coloc/protein_hit_table_task.py`, `mecfs_bio/assets/gwas/ukbb_ppp/ppp_database/hapmap3/eur_discovery_hapmap3_ppp_hit_table.py`
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_protein_hit_table_task.py`

**Interfaces:**
- Consumes: `BuildSlimProteinParquetTask` (`.protein.oid`, `.protein.gene`, slim columns BETA, SE in index row order); `HAPMAP_3_PPP_DATABASE_INDEX` (CHR, POS hg38).
- Produces: `PppProteinHitTableTask(GeneratingTask)` with
  `create(asset_id: str, protein_tasks: tuple[BuildSlimProteinParquetTask, ...], index_task: Task, p_threshold: float = 5e-8) -> PppProteinHitTableTask`
  (assert 0 < p_threshold < 1). Output `ResultTableMeta` parquet, one row per significant (protein, variant): PROTEIN_ID (OID), GENE, CHR, POS, MLOG10P (Float64). P computed from two-sided normal on BETA/SE in log space (avoid underflow: use `scipy.stats.norm.logsf`). Reads one protein at a time.
  Asset: `HAPMAP_3_PPP_HIT_TABLE_5E8 = PppProteinHitTableTask.create(asset_id="ppp_hit_table_hapmap_3_p5e-8", protein_tasks=HAPMAP_3_PPP_DATABASE.protein_tasks, index_task=HAPMAP_3_PPP_DATABASE_INDEX)`.

- [ ] **Step 1: failing Task-level tests** (fetch from file assets the way `test_ppp_protein_heritability_task.py` does):

```python
def test_hit_table_keeps_only_variants_below_threshold(tmp_path):   # 2 proteins, known z's
def test_threshold_is_configurable(tmp_path):                       # 1e-5 keeps more rows than 5e-8
def test_extreme_z_does_not_underflow(tmp_path):                     # z=40 → finite MLOG10P > 300
```

- [ ] **Steps 2–4: fail, implement, pass.**
- [ ] **Step 5: green (logged); commit** `"Add PPP protein hit table task with configurable threshold"`.

---

### Task 5: Trait loci and candidate selection (pure functions)

**Files:**
- Create: `mecfs_bio/build_system/task/protein_coloc/trait_loci.py`, `mecfs_bio/build_system/task/protein_coloc/candidate_selection.py`
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_trait_loci.py`, `test_candidate_selection.py`

**Interfaces:**
- Consumes: `GenomicInterval`, `extended_mhc_interval`; `read_gene_coords` output `dict[Oid, GenomicInterval]`; hit table from Task 4.
- Produces:
  - `find_trait_loci(trait: pl.DataFrame, p_threshold: float, half_window_bp: int, exclude: Sequence[GenomicInterval]) -> pl.DataFrame` — greedy: take min-p variant below threshold, window = lead ± half_window, remove variants inside, repeat; then merge overlapping windows (merged locus keeps the min-p lead). Columns: LOCUS_ID ("chr{c}:{start}-{end}"), CHR, START, END, LEAD_POS, LEAD_EA, LEAD_NEA, LEAD_MLOG10P, STATUS (`"candidate"` or `"excluded_mhc"` when the window overlaps any `exclude` interval).
  - `@frozen class SelectionConfig: max_pairs: int = 30; max_per_locus: int = 5; mlog10p_clip: float = 50.0; trans_penalty_weight: float = 2.0; cis_window_bp: int = 1_000_000`
  - `build_candidates(loci: pl.DataFrame, hit_table: pl.DataFrame, gene_coords: Mapping[Oid, GenomicInterval], cis_window_bp: int) -> pl.DataFrame` — one row per (non-excluded locus, protein with ≥1 hit in window): PROTEIN_ID, GENE, PROTEIN_MLOG10P (max in window), CIS_STATUS (`Literal["cis","trans","unknown"]`; unknown = no ST3 coordinates), N_OTHER_PROTEINS (distinct other proteins with hits in the window).
  - `select_candidates(candidates: pl.DataFrame, loci: pl.DataFrame, config: SelectionConfig) -> pl.DataFrame` — adds SCORE, SCORE_TRAIT, SCORE_PROTEIN, SCORE_PENALTY and STATUS (`"selected"` / `"not_selected_budget"`). Score = min(trait mlog10p, clip) + min(protein mlog10p, clip) − (w·log2(1 + N_OTHER) if CIS_STATUS ≠ "cis" else 0). Walk the ranking (score desc, then PROTEIN_ID, LOCUS_ID asc) taking a pair while selected < max_pairs and that locus has < max_per_locus.
  - Status values as a `ColocStatus` Literal + constants in `protein_coloc_constants.py`: candidate, excluded_mhc, selected, not_selected_budget, tested, failed.

- [ ] **Step 1: failing tests**

```python
def test_nearby_hits_merge_into_one_locus():
def test_mhc_window_marked_excluded():
def test_cis_pair_not_penalized_at_hotspot():      # cis protein with 50 other proteins outranks a trans one with equal p's
def test_unknown_gene_coords_treated_like_trans():
def test_per_locus_cap_respected():                  # 8 candidates at one locus, cap 5 → 5 selected
def test_pair_budget_respected_and_rest_marked():   # 40 candidates, budget 30 → 30 selected, 10 not_selected_budget
def test_clipping_limits_protein_term():            # mlog10p 300 and 60 score the same protein term
def test_selection_is_deterministic_under_ties():
```

- [ ] **Steps 2–4: fail, implement, pass.**
- [ ] **Step 5: green (logged); commit** `"Add trait locus finding and budgeted candidate selection"`.

---

### Task 6: The per-trait colocalization Task

**Files:**
- Create: `mecfs_bio/build_system/task/protein_coloc/trait_protein_coloc_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_trait_protein_coloc_task.py`

**Interfaces:**
- Consumes: Tasks 1–5; `load_fasta`, `load_panel_path`, `resolve_panel_af_col` (`genome_reference_harmonization_task.py`); `read_gene_coords`; `call_with_retries` (`mecfs_bio/util/retry`) around each protein fetch.
- Produces:
  - `@frozen class ColocRunConfig: lead_p_threshold: float = 5e-8; half_window_bp: int = 500_000; exclude_mhc: bool = True; selection: SelectionConfig = SelectionConfig(); min_aligned_variants: int = 50; harmonization: GenomeReferenceHarmonizationOptions = GenomeReferenceHarmonizationOptions()`
  - `TraitProteinColocTask(Task)` with `create(asset_id: str, trait_task: Task, trait_model: TraitModel, hit_table_task: PppProteinHitTableTask, gene_coords_task: Task, fasta_task: Task, panel_task: Task, panel_ancestry: PanelAncestry, fetcher: DenseWindowFetcher, method: ColocMethod = ColocAbf(), config: ColocRunConfig = ColocRunConfig()) -> TraitProteinColocTask`.
    `.create()` asserts the trait meta is `FilteredGWASDataMeta` with `harmonization_info.build == "38"`, that fasta/panel are build 38 (`resolve_harmonized_build`), and derives trait/project from the trait meta. The protein side is always `QuantitativeTraitModel(sd_y=1.0)` (PPP levels are inverse-rank normalized).
  - Output `ResultDirectoryMeta` directory containing `raw_results.parquet` (one row per locus × protein candidate plus one per excluded locus; columns per spec §Outputs (a): locus fields, protein fields, CIS_STATUS, score components, n_trait / n_protein / n_aligned / n_duplicate_keys_dropped, orientation drop counts as one column per drop reason (0 if absent), PP_H0..PP_H4, P1, P2, P12, METHOD, STATUS, FAILURE_REASON; integers Int32) and `windows/{locus_id}__{protein_id}.parquet` (the AlignedWindow frame of every tested pair, for plotting).
  - `execute()` reads as the algorithm: loci → candidates → selection → group selected pairs by protein → for each protein: fetch its windows, orient each, align with the trait window, coloc (or fail with reason if n_aligned < min_aligned_variants) → write outputs. Trait window rows come from the harmonized trait filtered to the window.

- [ ] **Step 1: failing Task-level tests** with a synthetic harmonized trait parquet (FilteredGWASDataMeta, harmonization_info build 38), a synthetic hit table, tiny FASTA + panel (as Task 3), and a fake `DenseWindowFetcher` implementation that returns prepared frames and counts calls per protein (inject; no mocking).

```python
def test_shared_signal_is_tested_with_high_h4(tmp_path):
    row = results.filter(STATUS == "tested", PROTEIN_ID == "OID1")
    assert row[PP_H4].item() > 0.8
def test_protein_at_two_loci_on_different_chromosomes_is_fetched_once(tmp_path):
    assert fake_fetcher.calls == {"OID1": 1} and results.filter(PROTEIN_ID == "OID1", STATUS == "tested").height == 2
def test_too_few_aligned_variants_marks_failed(tmp_path):
    assert row[STATUS].item() == "failed" and row[FAILURE_REASON].item() is not None
def test_unselected_and_mhc_rows_are_reported(tmp_path):
def test_create_rejects_unharmonized_trait():
    with pytest.raises(AssertionError): TraitProteinColocTask.create(trait_task=<meta without harmonization_info>, ...)
def test_window_parquet_written_for_each_tested_pair(tmp_path):
```

- [ ] **Steps 2–4: fail, implement, pass.**
- [ ] **Step 5: green (logged); commit** `"Add per-trait protein colocalization task"`.

---

### Task 7: Display table, GRCh38 gene track, and stacked locus plots

**Files:**
- Create: `mecfs_bio/build_system/task/protein_coloc/coloc_display_table_task.py`, `mecfs_bio/build_system/task/protein_coloc/coloc_locus_plot_task.py`, `mecfs_bio/assets/reference_data/magma_gene_locations/processed/magma_entrez_gene_locations_data_build_38_unzipped.py`
- Test: `test_coloc_display_table_task.py`, `test_coloc_locus_plot_task.py` (same test dir)

**Interfaces:**
- Consumes: Task 6 result directory; `write_parquet_table` (`dataframe_output.py`); `plot_gene_tracks`, `GENE_INFO_*` (`susie_stacked_plot_task.py`); `MAGMA_GENE_LOCATION_REFERENCE_DATA_BUILD_38_RAW` + `ExtractFromZipTask.create_from_zipped_reference_file` (mirror the build-37 unzipped module; the extracted `.loc` is whitespace-separated entrez_id, chrom, start, end, strand, symbol → map to `GENE_INFO_*` via its read spec).
- Produces:
  - `ColocDisplayTableTask.create(asset_id: str, coloc_task: TraitProteinColocTask, h4_threshold: float = 0.8)` → `ResultTableMeta` parquet: GENE, PROTEIN_ID, LOCUS (chr:start–end), LEAD_VARIANT, CIS_STATUS, PP_H4; sorted by PP_H4 desc; Int32 ints; zstd + byte-stream-split on PP_H4.
  - `ColocLocusPlotTask.create(asset_id: str, coloc_task: TraitProteinColocTask, display_task: ColocDisplayTableTask, gene_info_task: Task)` → `ResultDirectoryMeta` with `{locus_id}__{protein_id}.png` + `.svg` per display row, and `plot_index.parquet` (PROTEIN_ID, GENE, LOCUS_ID, PNG_FILENAME, SVG_FILENAME). Figure: OO matplotlib (`Figure`, no pyplot), three panels sharing x: trait −log10 p, protein −log10 p (each highlighting its lead), genes; title "GENE (cis|trans) PP.H4=0.93".

- [ ] **Step 1: failing tests**

```python
def test_display_keeps_only_rows_above_threshold_sorted(tmp_path):
def test_display_integer_columns_are_int32(tmp_path):
def test_empty_display_when_nothing_passes(tmp_path):        # height 0, schema intact
def test_plot_files_and_index_for_each_display_row(tmp_path):
def test_no_display_rows_gives_empty_index_and_no_images(tmp_path):
```

- [ ] **Steps 2–4: fail, implement, pass.**
- [ ] **Step 5: green (logged); commit** `"Add coloc display table and stacked locus plots"`.

---

### Task 8: GRCh38-harmonized trait inputs (four traits)

Wiring only (no structural tests, per repo convention); verified by import, typecheck, and the real
runs in Task 9. Each trait: GRCh38 gwaslab sumstats (no liftover, `genome_build="38"`) →
`GwasLabSumstatsToTableTask.create_from_source_task` → `GenomeReferenceHarmonizationTask.create(fasta_task=UCSC_HG38_INDEXED_FASTA, panel_task=GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES, panel_ancestry="nfe")`. Pattern: `liu_et_al_2023_eur_liftover_to_37_sumstats_harmonized.py`.

**Files (create, one module per trait under its `processed_gwas_data/` or `processed/` dir, named `*_build_38_harmonized.py`):**
- DecodeME: source `DECODE_ME_GWAS_1_SUMSTATS_MINIMAL_FILTERING` (native b38); trait model `CaseControlTraitModel(case_fraction=DECODE_ME_PREVALENCE_INFO.sample_prevalence)`.
- Seropositive RA: source `SEROPOS_RA_EAF_SCALED` (NOT `SEROPOS_RA_FILTERED_FOR_FREQ`); same pre-pipes as `ra_seropositive_standard_analysis.py` (ACGT filter, `EffectiveNFromCohortStringPipe`, then `ComputeBetaPipe`/`ComputeSEPipe`); pos `PosB38`; model `CaseControlTraitModel(case_fraction=0.5)` (effective N).
- Willer LDL: source `WILLER_LDL_EUR_DATA_RAW`; positions `hm_chrom`/`hm_pos`; ORIGINAL `effect_allele`/`other_allele` (uppercased by a pre-pipe) with `beta`, `standard_error`, `n`, and `eaf_ref` nulled for palindromic SNVs (placeholder values; spec) ; model `QuantitativeTraitModel(sd_y=None)`.
- IBD: source `LIU_ET_AL_2023_IBD_META_EUR_ONLY` with the existing column mapping, no liftover; model `CaseControlTraitModel(case_fraction=LIU_ET_AL_IBD_PREVALENCE_INFO.sample_prevalence)`.

- [ ] **Step 1:** Write the four modules. **Step 2:** `pixi r invoke green` (logged) passes, including import checks. **Step 3:** Build one (Willer, the smallest) via `DEFAULT_RUNNER.run([...])` and check the harmonizer's logged trust decision and that output build is 38. This also builds the gnomAD v4.1 panel the first time; if that build is impractically large, stop and report (fallback per spec: 1000G EUR hg38 30x panel, ancestry "eur").
- [ ] **Step 4: commit** `"Add GRCh38-harmonized inputs for the four coloc traits"`.

---

### Task 9: Coloc assets per trait, analysis runner, and real-run validation

**Files:**
- Create: `mecfs_bio/assets/gwas/ukbb_ppp/ppp_database/hapmap3/ppp_coloc_runs.py` (one `TraitProteinColocTask` + display + plot task per trait, `fetcher=PppWindowFetcher(proteins=<from EUR discovery manifest via the existing manifest-row helper in ukbb_ppp_slim_protein_asset_generator.py>)`, `gene_coords_task=SUN_ET_AL_2023_ST3_EXTRACTED`, `gene_info_task=` Task 7's build-38 gene locations), `mecfs_bio/analysis/ppp_coloc_analysis.py` (runs all terminal tasks).

- [ ] **Step 1:** Write modules; green (logged).
- [ ] **Step 2:** Run Willer LDL end to end (`pixi r python mecfs_bio/analysis/ppp_coloc_analysis.py` restricted to LDL, or a scratch runner). Tee the log under `experiments/claude/logs/` (uncommitted).
- [ ] **Step 3: Positive-control check.** At least one tested pair with PP.H4 > 0.8, CIS_STATUS "cis", at a classic lipid locus (e.g. PCSK9). If none, investigate the pipeline (orientation drop counts, n_aligned, sign agreement of lead betas) before blaming biology; report findings.
- [ ] **Step 4:** Run the other three traits; record per-trait counts (loci, candidates, selected, tested, failed, displayed) in the PR description or `experiments/claude/` notes.
- [ ] **Step 5: commit** `"Wire PPP coloc runs for four traits and add analysis runner"`, then hand back the suggested docs snippets (`{{ data_table(...) }}` per trait and the plot paths) to the user instead of editing `docs/`.
