# Fast All-Pairs CT-LDSC Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an asset generator that computes all-pairs CT-LDSC genetic correlations plus per-trait heritabilities for about 23 traits in a few minutes, reusing the ppp_ldsc kernels, with fail-fast guards against upstream rsID/build bugs.

**Architecture:** A shared HapMap3 ∩ LD-score variant index Task. One aligned-trait Task per trait, which reads a gwaslab pickle or a dataframe, aligns z to the index and enforces the guards. An all-pairs rg Task and a heritability Task, which stack the aligned traits into an (S, T) matrix and call `estimate_trait_context` / `batched_rg`. The kernels gain small backward-compatible extensions (intercept jackknife, per-SNP N on the second side, MHC build parameter, alignment counts).

**Tech Stack:** Python, numpy, polars, narwhals, attrs, the repo build system (Task/Meta/Fetch), pytest, pixi.

**Spec:** experiments/claude/design_specs/2026-10-25-fast-all-pairs-ct-ldsc-design.md

## Global Constraints

- All commands via pixi: `pixi r pytest <explicit unit paths>`; never pass the test_mecfs_bio root to pytest.
- After each task: `pixi r invoke green > <scratchpad>/green.log 2>&1; echo $?`, then read the log (exit 0 alone is not proof under testmon).
- PPP behaviour must not change: every existing test under test_mecfs_bio/unit/build_system/task/ppp_ldsc/ passes unmodified. The constant-N `batched_rg` path stays bit-identical.
- Estimator: GenomicSEM-style (free intercepts, chi² filter max(0.001·N_max, 80), contiguous delete-block jackknife), `n_blocks` default 200.
- Regression SNPs: HAPMAP3_SNPLIST_FOR_GENOMIC_SEM ∩ THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE, minus strand-ambiguous SNPs and the extended MHC in build "19"; M = total M_5_50 before the join.
- Repo conventions: imports at the top of the file; no monkeypatching or mocks (use FakeTask, inject dependencies); docstrings without backticks or RST; column-name constants rather than literals; @frozen(slots=True) attrs; invariants in __attrs_post_init__; no bare tuple returns; Literal types; no `match=` on pytest.raises that duplicates the error text; Task-level tests preferred; no structural tests for asset generators.
- Don't edit anything under docs/. Never commit logs or PDFs; stage explicit paths.
- The existing gwaslab generator, CT_LDSC_INITIAL_ASSET_GENERATOR, its plot and LAVA wiring keep their asset ids and behaviour.

## Review Focus

1. A trait whose N column is null on rows with a finite z: NaN would silently poison h² and n_bar. Expected: the aligned-trait task fails with the alias in the message (test in Task 6).
2. gwaslab pickles store EA/NEA (and sometimes rsID) as pandas Categoricals. Expected: the pickle path yields exactly the same aligned z as the parquet path for the same data (test in Task 6).
3. A trait whose h² estimate is ≤ 0. Expected: its pairs get NaN rg/se/z/p (no crash), and the heritability row is still written (test in Task 7).
4. A trait with correct rsIDs but positions in hg38 (or no declared build). Expected: the position guard is skipped, not tripped (test in Task 6).
5. Duplicate aliases, or aligned tasks built against a different index. Expected: an assert at create time (aliases) or at execute time (the rsID column differs from the context) (both tested in Task 7).

## Interface note: where position-guard "declared build" comes from

Metas generally do not carry a build, so the spec's "source meta declares build 19" resolves as follows (the builder's call, recorded here). For a gwaslab pickle, use `sumstats.build` when it is "19" or "38", else None. For a dataframe input, use `meta.harmonization_info.build` when the meta is a `FilteredGWASDataMeta` with harmonization_info, else None. The position guard applies only when the declared build is "19".

---

### Task 1: Intercept jackknife and heritability diagnostics in the rg kernel

**Files:**
- Modify: `mecfs_bio/build_system/task/ppp_ldsc/batched_ldsc_rg.py`
- Test: `test_mecfs_bio/unit/build_system/task/ppp_ldsc/test_batched_ldsc_rg.py`

**Interfaces:**
- Produces:
  - `BlockedWlsResult` (frozen attrs): `slope (K,)`, `intercept (K,)`, `slope_delete (n_blocks, K)`, `intercept_delete (n_blocks, K)`; `_blocked_wls(w, y, ld, n_blocks) -> BlockedWlsResult` (replaces the 3-tuple).
  - `_delete_jackknife_se(delete: np.ndarray) -> np.ndarray`: `sqrt((B-1)/B * sum((d - d.mean(0))**2, axis=0))` over axis 0. Use it for rg_se with the exact same operation order, so rg_se is unchanged.
  - `TraitLdscContext` gains `h2_se: float`, `intercept: float`, `intercept_se: float`, `mean_chi2: float` (over kept SNPs), `lambda_gc: float` (median kept chi² / 0.4549364231195724; reuse `_CHI2_1DF_MEDIAN` from batched_ldsc_h2), `n_snps: int` (kept count).
  - `BatchedRgResult` gains `gcov_se (K,)` and `gcov_intercept_se (K,)`, float-kind, validated in post_init.

- [ ] **Step 1: Write failing tests**
```python
def test_trait_context_reports_h2_se_and_intercept_matching_exact() -> None:
    # constant-N synthetic trait (reuse _synthetic), no NaNs
    ctx = estimate_trait_context(z, np.full(s, _N), ld, _M, n_blocks=_N_BLOCKS)
    exact = exact_h2_single(z**2, ld, _N, _M, n_blocks=_N_BLOCKS)
    assert ctx.intercept == pytest.approx(exact.intercept, rel=1e-8)
    assert ctx.h2_se == pytest.approx(exact.h2_se, rel=0.05)  # block membership differs slightly
    assert ctx.intercept_se > 0 and ctx.n_snps == int(ctx.keep.sum())

def test_rg_reports_finite_gcov_and_intercept_ses() -> None:
    res = batched_rg(ctx, z_protein, ld, n_protein, _M, n_blocks=_N_BLOCKS)
    assert np.all(np.isfinite(res.gcov_se)) and np.all(res.gcov_se > 0)
    assert np.all(np.isfinite(res.gcov_intercept_se)) and np.all(res.gcov_intercept_se > 0)
```
- [ ] **Step 2:** `pixi r pytest test_mecfs_bio/unit/build_system/task/ppp_ldsc/test_batched_ldsc_rg.py -v`. Expected: FAIL (AttributeError on the new fields).
- [ ] **Step 3:** Implement as above. Update the module docstring's "rg standard error" paragraph to mention the intercept jackknife. Update the one external `estimate_trait_context` caller only if needed (ppp_protein_rg_task.py doesn't construct the context directly).
- [ ] **Step 4:** Run all of `test_mecfs_bio/unit/build_system/task/ppp_ldsc/`. Expected: PASS, including the unmodified existing tests.
- [ ] **Step 5:** Commit: `git add` the two files; message "Add intercept jackknife and h2 diagnostics to batched LDSC rg kernel".

### Task 2: Per-SNP N for the second side of batched_rg

**Files:**
- Modify: `mecfs_bio/build_system/task/ppp_ldsc/batched_ldsc_rg.py`
- Test: `test_mecfs_bio/unit/build_system/task/ppp_ldsc/test_batched_ldsc_rg.py`

**Interfaces:**
- Consumes: Task 1's `TraitLdscContext`.
- Produces: `batched_rg(trait_ctx, z_protein, ld, n_protein, m, *, n_blocks, exclude)` where `n_protein` is `(K,)` (unchanged path) or `(S, K)` float with NaN allowed where z is NaN. For the 2-D case: assert `n_protein.shape == z_protein.shape`; `n_p = np.where(finite_p, n_protein, 0.0)`; threshold per column `_chisq_threshold(np.nanmax(n_protein, axis=0))`. This mirrors `estimate_trait_context`'s handling exactly. Update the docstring.

- [ ] **Step 1: Write failing tests**
```python
def test_per_snp_n_second_side_matches_its_own_trait_context() -> None:
    # Z (S, 3), N (S, 3) with N varying per SNP (e.g. 8_000..12_000) and ~5% NaN in each column
    ctx0 = estimate_trait_context(Z[:, 0], N[:, 0], ld, _M, n_blocks=_N_BLOCKS)
    res = batched_rg(ctx0, Z[:, 1:], ld, N[:, 1:], _M, n_blocks=_N_BLOCKS)
    for j in (1, 2):
        own = estimate_trait_context(Z[:, j], N[:, j], ld, _M, n_blocks=_N_BLOCKS)
        assert res.h2_protein[j - 1] == pytest.approx(own.h2, rel=1e-10)

def test_constant_matrix_n_equals_vector_n() -> None:
    vec = batched_rg(ctx, Zp, ld, np.array([_N, _N]), _M, n_blocks=_N_BLOCKS)
    mat = batched_rg(ctx, Zp, ld, np.where(np.isfinite(Zp), _N, np.nan), _M, n_blocks=_N_BLOCKS)
    np.testing.assert_allclose(mat.rg, vec.rg, rtol=1e-12)
    np.testing.assert_allclose(mat.rg_se, vec.rg_se, rtol=1e-12)
```
- [ ] **Step 2:** Run test_batched_ldsc_rg.py. Expected: FAIL (shape/broadcast error).
- [ ] **Step 3:** Implement the 2-D branch, leaving the 1-D branch byte-for-byte unchanged.
- [ ] **Step 4:** Run `test_mecfs_bio/unit/build_system/task/ppp_ldsc/`. Expected: PASS.
- [ ] **Step 5:** Commit: "Accept per-SNP sample size on the second side of batched_rg".

### Task 3: MHC build parameter and alignment counts

**Files:**
- Modify: `mecfs_bio/build_system/task/ppp_ldsc/ppp_ldsc_context.py`, `mecfs_bio/build_system/task/ppp_ldsc/trait_alignment.py`
- Test: `test_mecfs_bio/unit/build_system/task/ppp_ldsc/test_ppp_ldsc_context.py`, `test_mecfs_bio/unit/build_system/task/ppp_ldsc/test_trait_alignment.py`

**Interfaces:**
- Produces:
  - `build_ppp_ldsc_context(index_df, ld_df, *, drop_strand_ambiguous=True, exclude_mhc=True, build: GenomeBuild = "38")`, using `extended_mhc_interval(build)`. Fix the stale hg38-only comments.
  - `TraitAligned` gains `n_found: int` (context rsIDs present in the trait after dedup), `n_matched: int` (found and alleles match directly or swapped), `n_position_compared: int | None`, `n_position_concordant: int | None`. Post_init: `0 <= n_matched <= n_found <= len(z)`; the position fields are both None or both int, with `concordant <= compared`.
  - `align_trait_to_context` computes the position fields only when BOTH `context_variants` and the trait carry GWASLAB_CHROM_COL and GWASLAB_POS_COL. compared = matched rows with non-null trait CHR/POS after `cast(pl.Int64, strict=False)`; concordant = those whose (CHR, POS) equal the context's. Otherwise both are None. PPP callers pass no CHR/POS, so they get None, and `min_trait_snps` behaves as before.

- [ ] **Step 1: Write failing tests**
```python
def test_context_excludes_hg19_mhc_when_build_19() -> None:
    # index row at chr6 pos inside EXTENDED_MHC_BUILD_37 but outside the build-38 interval
    ctx = build_ppp_ldsc_context(index, ld, build="19")
    assert mhc_rsid not in kept_rsids(ctx)

def test_alignment_counts_found_and_matched() -> None:
    # reuse the orientation test's trait: one match, one flip, one allele mismatch, one absent
    aligned = align_trait_to_context(...)
    assert (aligned.n_found, aligned.n_matched) == (3, 2)
    assert aligned.n_position_compared is None

def test_alignment_position_concordance() -> None:
    # context_variants with CHR/POS; trait with CHR/POS where 1 of 2 matched rows has a shifted POS
    aligned = align_trait_to_context(...)
    assert (aligned.n_position_compared, aligned.n_position_concordant) == (2, 1)
```
- [ ] **Step 2:** Run both test files. Expected: FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run `test_mecfs_bio/unit/build_system/task/ppp_ldsc/`. Expected: PASS.
- [ ] **Step 5:** Commit: "Add MHC build parameter and alignment counts to ppp_ldsc helpers".

### Task 4: Constants, alignment guards, SumstatsSource override

**Files:**
- Create: `mecfs_bio/constants/fast_ct_ldsc_constants.py`, `mecfs_bio/build_system/task/fast_ct_ldsc/__init__.py`, `mecfs_bio/build_system/task/fast_ct_ldsc/alignment_guards.py`
- Modify: `mecfs_bio/build_system/task/gwaslab/gwaslab_genetic_corr_by_ct_ldsc_task.py` (SumstatsSource)
- Test: `test_mecfs_bio/unit/build_system/task/fast_ct_ldsc/__init__.py`, `test_mecfs_bio/unit/build_system/task/fast_ct_ldsc/test_alignment_guards.py`

**Interfaces:**
- Produces, in constants (values are what downstream readers require; the rest are the builder's call but must be constants):
  - Pairs: `FAST_CT_LDSC_P1_COL="p1"`, `FAST_CT_LDSC_P2_COL="p2"`, `FAST_CT_LDSC_RG_COL="rg"`, `FAST_CT_LDSC_RG_SE_COL="se"`, `FAST_CT_LDSC_RG_Z_COL="z"`, `FAST_CT_LDSC_RG_P_COL="p"`, `FAST_CT_LDSC_GCOV_COL="gcov"`, `FAST_CT_LDSC_GCOV_SE_COL="gcov_se"`, `FAST_CT_LDSC_GCOV_INT_COL="gcov_int"`, `FAST_CT_LDSC_GCOV_INT_SE_COL="gcov_int_se"`, `FAST_CT_LDSC_H2_OBS_P1_COL="h2_obs_p1"`, `FAST_CT_LDSC_H2_OBS_P2_COL="h2_obs_p2"`, `FAST_CT_LDSC_H2_INT_COL="h2_int"`, `FAST_CT_LDSC_N_SNPS_COL="n_snps"`.
  - Heritability: `FAST_CT_LDSC_TRAIT_COL="p"`, `FAST_CT_LDSC_H2_OBS_COL="h2_obs"`, `..._H2_OBS_SE_COL="h2_obs_se"`, `..._H2_LIAB_COL="h2_liab"`, `..._H2_LIAB_SE_COL="h2_liab_se"`, `..._INTERCEPT_COL="Intercept"`, `..._INTERCEPT_SE_COL="Intercept_se"`, `..._LAMBDA_GC_COL="lambda_gc"`, `..._MEAN_CHI2_COL="mean_chi2"`, `..._N_BAR_COL="n_bar"`, `..._COVERAGE_COL="coverage"`, `..._ALLELE_MISMATCH_FRACTION_COL="allele_mismatch_fraction"`, `..._POSITION_CONCORDANCE_COL="position_concordance"`. Also reuse `FAST_CT_LDSC_N_SNPS_COL`.
  - Aligned-trait files: `FAST_CT_LDSC_Z_COL="z"`, `FAST_CT_LDSC_N_COL="n"`, `FAST_CT_LDSC_ALIGNED_FILENAME="aligned.parquet"`, `FAST_CT_LDSC_DIAGNOSTICS_FILENAME="diagnostics.json"`.
  - HapMap3 snplist columns: `HAPMAP3_SNPLIST_RSID_COL="SNP"`, `HAPMAP3_SNPLIST_A1_COL="A1"`, `HAPMAP3_SNPLIST_A2_COL="A2"`.
- Produces, in alignment_guards.py:
  - `AlignmentGuards` (frozen): `min_coverage: float = 0.5`, `max_allele_mismatch_fraction: float = 0.02`, `min_position_concordance: float = 0.95`; post_init asserts each is in [0, 1]. These defaults are provisional; Task 9 replaces them.
  - `AlignmentDiagnostics` (frozen): `n_context: int`, `n_found: int`, `n_matched: int`, `n_position_compared: int | None`, `n_position_concordant: int | None`, `declared_build: GenomeBuild | None`; properties `coverage -> float` (n_matched/n_context), `allele_mismatch_fraction -> float` (0.0 when n_found == 0), `position_concordance -> float | None` (None when not compared or compared == 0). `to_json_dict() -> dict[str, int | str | None]` and `@classmethod from_json_dict(d) -> AlignmentDiagnostics`.
  - `alignment_diagnostics(aligned: TraitAligned, declared_build: GenomeBuild | None) -> AlignmentDiagnostics`.
  - `check_alignment_guards(diagnostics: AlignmentDiagnostics, guards: AlignmentGuards, alias: str) -> None`: asserts in order coverage, mismatch, then position (only when declared_build == "19" and position_concordance is not None). Each message names the alias, the diagnostic, its value and the threshold.
- Produces, in SumstatsSource: a new last field `alignment_guards: AlignmentGuards | None = None`. The gwaslab task ignores it.

- [ ] **Step 1: Write failing tests** (pure functions, small):
```python
def test_guards_pass_on_clean_diagnostics() -> None: check_alignment_guards(clean, AlignmentGuards(), "T")
@pytest.mark.parametrize("diag", [low_coverage, high_mismatch, low_concordance_build19])
def test_guards_trip(diag) -> None:
    with pytest.raises(AssertionError): check_alignment_guards(diag, AlignmentGuards(), "T")
def test_position_guard_skipped_without_declared_build_19() -> None:  # same counts, declared_build None and "38"
def test_diagnostics_json_roundtrip() -> None: assert AlignmentDiagnostics.from_json_dict(d.to_json_dict()) == d
def test_guard_thresholds_validated() -> None:
    with pytest.raises(AssertionError): AlignmentGuards(min_coverage=1.5)
```
- [ ] **Step 2:** Run the test file. Expected: FAIL (ImportError).
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run the file plus `test_mecfs_bio/unit/build_system/task/gwaslab/test_gwaslab_genetic_corr_by_ct_ldsc_task.py`. Expected: PASS.
- [ ] **Step 5:** Commit: "Add fast CT-LDSC constants and alignment guards".

### Task 5: Shared variant index Task and asset

**Files:**
- Create: `mecfs_bio/build_system/task/fast_ct_ldsc/variant_index_task.py`, `mecfs_bio/assets/reference_data/linkage_disequilibrium_score_reference_data/extracted/hapmap3_eur_ldsc_variant_index.py`
- Test: `test_mecfs_bio/unit/build_system/task/fast_ct_ldsc/test_variant_index_task.py`

**Interfaces:**
- Consumes: Task 3's `build_ppp_ldsc_context(..., build="19")`; Task 4's HapMap3 constants.
- Produces:
  - `FastCtLdscVariantIndexTask(Task)` (frozen): `meta: Meta`, `hapmap3_task: Task`, `ld_scores_task: Task`; `deps = [hapmap3_task, ld_scores_task]`. execute: inner-join the snplist (SNP, A1, A2) to the LD scores (CHR, SNP, BP) on rsID. Write parquet with GWASLAB_CHROM_COL, GWASLAB_POS_COL (from BP), GWASLAB_RSID_COL, GWASLAB_EFFECT_ALLELE_COL (= A1), GWASLAB_NON_EFFECT_ALLELE_COL (= A2), PPP_INDEX_IS_STRAND_AMBIGUOUS_COL (`is_palindromic_expr`), sorted by (CHR, POS, rsID). `create(asset_id, hapmap3_task, ld_scores_task)` uses `ReferenceFileMeta(group="linkage_disequilibrium_score_reference_data", sub_group="fast_ct_ldsc", sub_folder=PurePath("processed"), id=..., extension=".parquet", read_spec=DataFrameReadSpec(DataFrameParquetFormat()))`. Assert both deps `produces_dataframe`.
  - `IndexedLdscContext` (frozen): `context: PppLdscContext`, `variants: pl.DataFrame` (the index rows gathered by `context.row_pos`: rsID, EA, NEA, CHR, POS; post_init asserts the height equals `context.n_snps`).
  - `load_indexed_ldsc_context(index_task: FastCtLdscVariantIndexTask, fetch: Fetch) -> IndexedLdscContext`: reads the index and `index_task.ld_scores_task` and calls `build_ppp_ldsc_context(index_df, ld_df, build="19")`.
  - Asset: `HAPMAP3_EUR_LDSC_VARIANT_INDEX = FastCtLdscVariantIndexTask.create("hapmap3_eur_ldsc_variant_index", HAPMAP3_SNPLIST_FOR_GENOMIC_SEM, THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE)`.

- [ ] **Step 1: Write failing tests** (FakeTask + parquet/tsv inputs, `make_wf()`):
```python
def test_index_joins_snplist_to_ld_scores_with_a1_as_effect_allele() -> None:
    # snplist has one rsID absent from LD scores; LD has one rsID absent from the snplist
    assert out[GWASLAB_RSID_COL].to_list() == expected_sorted_rsids
    assert out[GWASLAB_EFFECT_ALLELE_COL].to_list() == expected_a1
def test_loaded_context_drops_hg19_mhc_and_strand_ambiguous() -> None:
    ictx = load_indexed_ldsc_context(task, fetch)
    assert mhc_rsid not in ictx.variants[GWASLAB_RSID_COL].to_list()
    assert palindromic_rsid not in ictx.variants[GWASLAB_RSID_COL].to_list()
```
- [ ] **Step 2:** Run. Expected: FAIL.
- [ ] **Step 3:** Implement the task, the helper and the asset.
- [ ] **Step 4:** Run. Expected: PASS.
- [ ] **Step 5:** Commit: "Add shared HapMap3 LDSC variant index task".

### Task 6: Aligned-trait Task (both input formats, guards)

**Files:**
- Create: `mecfs_bio/build_system/task/fast_ct_ldsc/aligned_trait_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/fast_ct_ldsc/test_aligned_trait_task.py`

**Interfaces:**
- Consumes: Tasks 3–5.
- Produces:
  - `FastCtLdscAlignedTraitTask(Task)` (frozen): `meta: Meta`, `source: SumstatsSource`, `index_task: FastCtLdscVariantIndexTask`, `guards: AlignmentGuards`; `deps = [source.task, index_task, index_task.ld_scores_task]`; property `alias -> str` (source.alias).
  - execute (as a high-level spec):
    load context → `_load_trait(source, fetch) -> LoadedTrait(frame: narwhals.LazyFrame, declared_build: GenomeBuild | None)` → `source.pipe.process` → select rsID/EA/NEA/BETA/SE plus N/CHR/POS when present → `align_trait_to_context(..., trait_total_sample_size=source.sample_info.total_sample_size, min_trait_snps=0)` → `alignment_diagnostics` → `check_alignment_guards(diag, source.alignment_guards or self.guards, alias)` → assert `n` is finite wherever `z` is finite → write `<dir>/aligned.parquet` (rsID, z, n in context row order) and `<dir>/diagnostics.json` → return `DirectoryAsset`.
  - `_load_trait`: for GWASLabSumStatsMeta, `read_sumstats` → `pl.from_pandas(sumstats.data)` with Categorical/Enum columns cast to `pl.String` → lazy narwhals; declared build per the interface note. For a dataframe meta, `scan_dataframe_asset(..., parquet_backend="polars")`; declared build per the interface note.
  - `create(asset_id, source, index_task, guards)`: asserts the source meta is `GWASLabSumStatsMeta`, or that `produces_dataframe(source.task)` holds and the meta is `FilteredGWASDataMeta | GWASSummaryDataFileMeta`; meta = `ResultDirectoryMeta(id=asset_id, trait=src.trait, project=src.project, sub_dir=PurePath("analysis/fast_ct_ldsc_aligned"))`.
  - `read_aligned_trait(task: FastCtLdscAlignedTraitTask, fetch) -> AlignedTraitArrays(rsid: pl.Series, z: np.ndarray, n: np.ndarray, diagnostics: AlignmentDiagnostics)`, used by Task 7.

- [ ] **Step 1: Write failing tests** (synthetic index of about 400 SNPs via a FakeTask index; traits as FakeTask parquet/pickle):
```python
def test_parquet_trait_is_aligned_with_flips_and_writes_diagnostics() -> None:
    # half the rows have EA/NEA swapped and BETA negated -> aligned z identical to unswapped truth
    np.testing.assert_allclose(arrays.z, expected_z)
    assert arrays.diagnostics.n_matched == n_present
def test_gwaslab_pickle_matches_parquet_input() -> None:      # Review Focus 2
    np.testing.assert_array_equal(from_pickle.z, from_parquet.z)
def test_shuffled_rsids_trip_mismatch_guard() -> None:         # pytest.raises(AssertionError)
def test_sparse_trait_trips_coverage_guard() -> None:
def test_shifted_positions_trip_concordance_guard_when_declared_build_19() -> None:
    # dataframe meta = FilteredGWASDataMeta(harmonization_info=HarmonizationInfo(build="19", ...))
def test_shifted_positions_pass_without_declared_build() -> None:   # Review Focus 4
def test_per_source_override_lets_a_tripping_trait_through() -> None:
def test_null_n_with_finite_z_fails() -> None:                 # Review Focus 1
def test_aligned_rsids_written_in_context_order() -> None:
    assert arrays.rsid.to_list() == ictx.variants[GWASLAB_RSID_COL].to_list()
```
- [ ] **Step 2:** Run. Expected: FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run. Expected: PASS.
- [ ] **Step 5:** Commit: "Add fast CT-LDSC aligned-trait task".

### Task 7: Trait matrix, all-pairs rg Task, heritability Task

**Files:**
- Create: `mecfs_bio/build_system/task/fast_ct_ldsc/trait_matrix.py`, `mecfs_bio/build_system/task/fast_ct_ldsc/all_pairs_rg_task.py`, `mecfs_bio/build_system/task/fast_ct_ldsc/heritability_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/fast_ct_ldsc/test_all_pairs_rg_task.py`

**Interfaces:**
- Consumes: Tasks 1, 2, 5, 6; `liability_h2` from gwaslab/ldsc_diagnostic_plot_task.py; `MULTI_TRAIT` from gwaslab_genetic_corr_by_ct_ldsc_task.py.
- Produces:
  - `TraitMatrix` (frozen): `context: PppLdscContext`, `aliases: tuple[str, ...]`, `z: (S, T)`, `n: (S, T)`, `diagnostics: tuple[AlignmentDiagnostics, ...]`, `sample_infos: tuple[PhenotypeInfo, ...]`; post_init checks shapes and lengths.
  - `load_trait_matrix(index_task, aligned_tasks: Sequence[FastCtLdscAlignedTraitTask], fetch) -> TraitMatrix`: asserts each aligned rsID column equals the context variants' rsIDs.
  - `estimate_trait_contexts(matrix: TraitMatrix, *, n_blocks: int) -> list[TraitLdscContext]`.
  - `FastCtLdscAllPairsRgTask(Task)` (frozen): `meta`, `index_task`, `aligned_tasks: tuple[FastCtLdscAlignedTraitTask, ...]`, `n_blocks: int = 200`, `batch_size: int = 8`. Post_init: at least 2 traits, unique aliases, aliases sorted by `str.lower`. `deps = [index_task, index_task.ld_scores_task, *aligned_tasks]`. execute: load → estimate → `all_pairs_table(matrix, contexts, n_blocks=..., batch_size=...) -> pl.DataFrame` → parquet. Row (i<j): p1 = aliases[i], p2 = aliases[j]; rg/se/z/p, gcov/gcov_se/gcov_int/gcov_int_se from `batched_rg(ctx_i, z[:, J], ld, n[:, J], m)` over column batches J of traits > i; h2_obs_p1 = ctx_i.h2; h2_obs_p2 = ctx_j.h2; h2_int = ctx_j.intercept; n_snps. `create(asset_id, index_task, aligned_tasks, n_blocks, batch_size)` uses `ResultTableMeta(id, trait=MULTI_TRAIT, project="genetic_correlation", sub_dir=PurePath("fast_ct_ldsc"), extension=".parquet", read_spec=DataFrameReadSpec(DataFrameParquetFormat()))`.
  - `FastCtLdscHeritabilityTask(Task)` (frozen): `meta`, `index_task`, `aligned_tasks`, `n_blocks: int = 200`; same deps and meta pattern. `heritability_table(matrix, contexts) -> pl.DataFrame` has one row per trait with the Task 4 heritability columns. Liability factor = `liability_h2(1.0, info)` (the conversion is linear): h2_liab = factor·h2_obs, h2_liab_se = factor·h2_obs_se, null for quantitative traits. The diagnostics columns come from `AlignmentDiagnostics`.

- [ ] **Step 1: Write failing tests** (synthetic: S = 20_000 SNPs on 2 chromosomes, LD scores drawn from a fixed distribution, 4 traits with z ~ multivariate normal so that true rg(A,B) = 0.6 and C is independent; D is pure noise (h² ≤ 0); trait B has per-SNP N; trait A is binary; fixed seed; build aligned-trait tasks via FakeTask fetch of pre-written aligned dirs, or run Task 6 for real):
```python
def test_all_pairs_outputs_each_unordered_pair_once() -> None:
    assert set(zip(df["p1"], df["p2"])) == {("A","B"),("A","C"),("A","D"),("B","C"),("B","D"),("C","D")}
def test_correlated_pair_recovers_rg() -> None:
    assert row("A","B")["rg"] == pytest.approx(0.6, abs=0.15); assert abs(row("A","C")["rg"]) < 0.2
def test_nonpositive_h2_trait_gives_nan_rg_without_crash() -> None:   # Review Focus 3
    assert all(np.isnan(r["rg"]) for r in rows_with("D"))
def test_pairs_and_heritability_tables_agree() -> None:
    # h2_obs_p2 == heritability h2_obs[p2]; h2_int == Intercept[p2]
def test_liability_columns_only_for_binary() -> None:
    # A: h2_liab_se / h2_obs_se == liability_h2(1.0, info_A); C: null
def test_outputs_feed_clustermap_and_lava_overlap_readers() -> None:
    ds = load_xr_corr_dataset(GeneticCorrSource(task=pairs_fake_task), fetch)
    assert ds[XR_GENETIC_CORR_ARRAY].shape == (4, 4)
    path = _get_sample_overlap_path(pairs_fake_task, herit_fake_task, fetch, ["A","B","C"], tmp_path)
    assert path is not None
def test_aligned_trait_from_a_different_index_rejected() -> None:   # Review Focus 5
    # aligned dir whose rsID column is a permutation of the context's
    with pytest.raises(AssertionError): load_trait_matrix(index_task, (bad_aligned, ok_aligned), fetch)
def test_duplicate_aliases_rejected() -> None:                 # Review Focus 5
    with pytest.raises(AssertionError): FastCtLdscAllPairsRgTask.create(..., aligned_tasks=(a, a_dup))
```
- [ ] **Step 2:** Run. Expected: FAIL.
- [ ] **Step 3:** Implement. Keep the per-pair loops in helpers so execute reads as a spec.
- [ ] **Step 4:** Run the file plus test_genetic_correlation_clustermap_plot.py. Expected: PASS. If the synthetic tolerance can't be hit with S = 20_000, raise S rather than loosen the tolerance.
- [ ] **Step 5:** Commit: "Add fast CT-LDSC all-pairs rg and heritability tasks".

### Task 8: Generator and assets

**Files:**
- Create: `mecfs_bio/asset_generator/fast_genetic_correlation_asset_generator.py`, `mecfs_bio/assets/gwas/multi_trait/genetic_correlation/ct_ldsc/ct_ldsc_initial_sources.py`, `mecfs_bio/assets/gwas/multi_trait/genetic_correlation/ct_ldsc/ct_ldsc_fast_initial_asset_generator.py`, `mecfs_bio/assets/gwas/multi_trait/genetic_correlation/ct_ldsc/ct_ldsc_fast_plot.py`, `mecfs_bio/analysis/fast_ct_ldsc_analysis.py`
- Modify: `mecfs_bio/assets/gwas/multi_trait/genetic_correlation/ct_ldsc/ct_ldsc_initial_asset_generator.py` (import the sources list; asset ids unchanged)

**Interfaces:**
- Produces:
  - `FastGeneticCorrTasks` (frozen): `aligned_trait_tasks: Mapping[str, FastCtLdscAlignedTraitTask]`, `pairs_task`, `heritability_task`, `pairs_markdown_task`; `terminal_tasks() -> list[Task]` returns [pairs, heritability, markdown].
  - `fast_genetic_corr_by_ct_ldsc_asset_generator(base_name: str, sources: Sequence[SumstatsSource], variant_index_task: FastCtLdscVariantIndexTask = HAPMAP3_EUR_LDSC_VARIANT_INDEX, guards: AlignmentGuards = AlignmentGuards(), n_blocks: int = 200, batch_size: int = 8) -> FastGeneticCorrTasks`. It sorts sources by `alias.lower()`. Asset ids: `f"{base_name}_aligned_{alias.lower()}"`, `f"{base_name}_fast_ct_ldsc_pairs"`, `f"{base_name}_fast_ct_ldsc_heritability"`, `f"{base_name}_fast_ct_ldsc_pairs_markdown"`. The markdown uses `ConvertDataFrameToMarkdownTask.create_from_result_table_task` with `DropColPipe([gcov, gcov_se, h2_obs_p1, h2_obs_p2, h2_int])` + `FormatFloatNumbersPipe(.4g)` on rg, se, z, p, gcov_int, gcov_int_se.
  - `CT_LDSC_INITIAL_SOURCES: list[SumstatsSource]`, moved verbatim with its comments.
  - `CT_LDSC_FAST_INITIAL_ASSET_GENERATOR = fast_genetic_corr_by_ct_ldsc_asset_generator("initial_rg_fast", sources=CT_LDSC_INITIAL_SOURCES)`.
  - `CT_LDSC_FAST_INITIAL_PLOT`: same as ct_ldsc_plot.py with `task=...pairs_task` and asset id "initial_genetic_correlation_by_fast_ct_ldsc_plot".
  - analysis script: `DEFAULT_RUNNER.run(gen.terminal_tasks() + [CT_LDSC_FAST_INITIAL_PLOT], incremental_save=True)` under `if __name__ == "__main__"`.

- [ ] **Step 1:** Implement (no structural test per convention; import-time construction is the check).
- [ ] **Step 2:** `pixi r python -c "import mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_fast_plot, mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_initial_asset_generator"`. Expected: no error. Confirm with a one-off check that the old generator's asset ids are unchanged: `sorted(t.asset_id for t in CT_LDSC_INITIAL_ASSET_GENERATOR.corr_tasks.values())` is the same before and after (record the "before" output prior to editing).
- [ ] **Step 3:** `pixi r invoke green > <scratchpad>/green.log 2>&1; echo $?`. Expected: 0, and the log shows the import, asset-import and typecheck steps passing.
- [ ] **Step 4:** Commit: "Add fast CT-LDSC asset generator and initial assets".

### Task 9: Real-data guard calibration and comparison experiment

**Files:**
- Create: `experiments/claude/fast_ct_ldsc_vs_gwaslab/compare.py` (tees its log to compare.log in the same dir, uncommitted)
- Modify: `mecfs_bio/build_system/task/fast_ct_ldsc/alignment_guards.py` (defaults), and any SumstatsSource in ct_ldsc_initial_sources.py that genuinely needs an override (with a comment saying why)

**Setup:** The asset-store config `default_runner_config.yaml` is gitignored and lives at the repo root. Symlink it into the worktree before running. The user has a run going in the repo root, so **ask the user before starting real-data builds**: they share the asset store, and two runners writing the same assets could collide.

- [ ] **Step 1:** Build the aligned-trait tasks only, with the guards temporarily permissive (`AlignmentGuards(0, 1, 0)` passed via the generator in the script). Tabulate coverage, allele_mismatch_fraction, position_concordance and declared_build per trait.
- [ ] **Step 2:** Set the defaults in AlignmentGuards with comfortable margin around the observed clean values. Add per-source overrides only for documented outliers. Update the Task 4 tests if any literal default appears there.
- [ ] **Step 3:** Run the fast generator end to end, recording wall time. Load the existing gwaslab aggregation and heritability outputs (build them if they're absent, which is slow). Tabulate per pair rg_fast, rg_gwaslab, Δrg, Δrg/se_gwaslab; per trait h2 and Δh2; and Bonferroni calls under each.
- [ ] **Step 4:** `pixi r invoke green > <scratchpad>/green.log 2>&1; echo $?`. Expected: 0.
- [ ] **Step 5:** Commit compare.py and the defaults change (not the log): "Calibrate fast CT-LDSC guards on real traits; add gwaslab comparison".
- [ ] **Step 6:** Report to the user: wall time, max |Δrg/se|, any Bonferroni call changes, and the chosen guard defaults with the observed values behind them.
