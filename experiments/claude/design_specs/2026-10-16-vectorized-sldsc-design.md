# Vectorized cell-type-specific S-LDSC — Design

Date: 2026-10-16
Branch: vectorized-sldsc

## 1. Goal and success criteria

The cell-type-specific S-LDSC Task (mecfs_bio/build_system/task/gwaslab/gwaslab_cell_analysis_by_sldsc.py,
`CellAnalysisByLDSCTask`), driven by `standard_sldsc_task_generator`
(mecfs_bio/asset_generator/concrete_sldsc_generator.py), takes many hours per trait. It runs one
regression per row of six LDSC-SEG `.ldcts` datasets (1,015 cell/tissue types in total: Multi_tissue_chromatin
489, ImmGen 292, Multi_tissue_gene_expr 205, GTEx_brain 13, Corces_ATAC 13, Cahoy 3).

We will add a new, vectorized S-LDSC Task that:

1. In its default `"exact"` mode, reproduces gwaslab's `cell_type_specific` output (Name, Coefficient,
   Coefficient_std_error, Coefficient_P_value) to floating-point tolerance.
2. Optionally runs a `"shared_controls"` fast-inexact mode whose speed and accuracy cost we measure.
3. Is agnostic to the source of the S-LDSC annotations: it consumes a format-neutral "prepared annotation
   set" asset, produced by per-format adapter Tasks. LDSC-SEG `.ldcts` directories are the first adapter.

Success criteria (measured on DecodeME, all six standard sets):

- exact mode: max |Δz| < 1e-4 vs the existing gwaslab results and no change in FDR calls at α = 0.01;
- end-to-end wall time per trait measured in minutes, not hours;
- peak RSS comfortably under 16 GB (the CI runner budget, which proxies laptop runnability).

shared_controls is measured, not gated; its numbers decide whether it becomes the default in a later change.

## 2. Why the existing Task is slow (hypothesis, to be confirmed by profiling)

Reading gwaslab/extension/ldsc/ldsc_sumstats.py `cell_type_specific`: for every `.ldcts` row it

- re-reads ~44 gzipped per-chromosome CSVs (target + control, 22 chromosomes) via pandas
  (`_read_chr_split_files` → `ldscore_fromlist`), sorts by CHR/BP, and string-merges on SNP against the
  sumstats;
- re-reads the control ("all genes"/ENTEX.control etc.) LD scores, identical for every row of a set;
- runs `reg.Hsq(..., twostep=None, old_weights=True)`: a single weighted least squares + 200-block
  jackknife, ~3e9 flops — well under a second.

So the hours are expected to be I/O and parsing, not linear algebra. The validation experiment (§7.2)
profiles the old Task to confirm this.

## 3. Exact computation being reproduced

From gwaslab's `cell_type_specific` + `_read_ld_sumstats` + `LD_Score_Regression.__init__` (old_weights
branch) + `LstsqJackknifeFast`:

Sumstats preprocessing (gwaslab `_estimate_h2_cts_by_ldsc` wrapper):

- HapMap3 restriction (gwaslab `Sumstats.estimate_h2_cts_by_ldsc` → `_get_hapmap3`, added during
  implementation): keep rows whose rsID is in gwaslab's bundled HapMap3 list for the sumstats' build and
  whose EA/NEA equal the HapMap3 A1/A2 in either orientation. The new Task reproduces this with
  gwaslab's own `_get_hapmap_full_polars(build, include_alleles=True)`; sumstats without an rsID column
  are rejected (gwaslab's CHR/POS fallback is not reproduced).
- gwaslab fills missing N with the median N and then casts every variant's N to int64
  (`.apply("int64")`, truncation toward zero, so a non-integer effective N such as 41,236.8 becomes
  41,236). The new Task reproduces the cast exactly. It deliberately does NOT reproduce the median fill:
  a missing N on any regression SNP (after alignment to the baseline SNP set) raises an error, because a
  HapMap3 SNP without a sample size indicates an upstream data-processing problem that should not be
  glossed over. This is the one intentional divergence from gwaslab; it never changes numbers on inputs
  where both succeed.
- Z = BETA/SE if Z absent; rsID renamed to SNP.
- dropna on required columns.

Merging and filtering:

- `smart_merge(ref_ld, sumstats)` then merge with w_ld: rows follow the baseline LD-score order
  (sorted by CHR, BP). The exact row order is pinned down empirically in validation and reproduced.
- Zero-variance baseline LD columns are dropped (`novar_cols`); reproduce this.
- Minimum SNP count (no gwaslab counterpart, fail-fast): if fewer than `min_hapmap3_snps` sumstats SNPs
  fall in the baseline's HapMap3 regression-SNP set (counted after the baseline/weights join, before the
  χ² filter), the Task raises. Default 100,000; tests pass a small value. Too few HapMap3 SNPs usually
  indicates an upstream processing problem (for example, lost rsIDs).
- χ² filter: keep Z² < max(0.001·max(N), 80) (strict inequality).
- Every regression SNP must have cts LD scores (ldsc raises ValueError otherwise).
- n_blocks = min(n_snp, 200); separators floor(linspace(0, n_snp, n_blocks + 1)).

Per cell type c, on the filtered rows:

- Unscaled LD columns L_c = [t_c, controls…, baseline…]; M_c = [M(t_c), M(controls)…, M(baseline)…]
  where each M is the genome-wide M_5_50 (sum over chromosomes, `M_fromlist`).
- x_tot,c = row-sum of L_c; M_tot,c = sum of M_c.
- h²_agg,c = M_tot,c · (mean χ² − 1) / mean(x_tot,c · N) (intercept = 1), clipped to [0, 1].
- w_c = 1 / (2 (1 + h²_agg,c · N · max(x_tot,c, 1) / M_tot,c)²) · 1 / max(w_ld, 1).
- Design X_c = L_c · N / N̄ with an appended intercept column of ones; y = χ².
- ldsc multiplies X and y by √w_c normalized to sum 1; the scalar cancels, so this is WLS with
  effective weight w_c (the same product-weight point as batched_ldsc_h2).
- Jackknife: per-block XᵀWX and XᵀWy; est = solve(totals); delete_j = solve(totals − block_j);
  pseudo = n_blocks·est − (n_blocks − 1)·delete; cov = cov(pseudo, ddof=1) / n_blocks.
- Report coef = est[0] / N̄, se = √cov[0,0] / N̄, p = norm.sf(coef / se) (one-sided). Rows sorted by p.

## 4. Data contract (prepared reference assets)

Defined once in `mecfs_bio/build_system/task/sldsc/prepared_sldsc_reference.py`: filename constants and
`@frozen` loader/writer types used by every adapter and by the consumer.

### 4.1 Prepared annotation set (directory asset, one per source dataset)

- `ld_scores.parquet`: one row per SNP with `SNP`, `CHR`, `BP`, sorted by CHR then BP; one float32 column
  per unique LD-score source (a control shared by many rows is stored once). Column names are internal ids.
  float32 is lossless relative to the 3-decimal precision of the source files and halves disk and memory.
- `m_5_50.parquet`: one row per LD column: column id, genome-wide M_5_50, and the source path it came from.
- `models.parquet`: one row per analysis label: `label`, `target_column`, `control_columns` (list[str]).
  Invariant (asserted on load/construction): every referenced column exists in the LD matrix; labels are
  unique; the target is not among its own controls.

### 4.2 Prepared baseline (directory asset, one per baseline/weights pairing)

- `ld_scores.parquet`: SNP, CHR, BP + baseline LD columns (~53 for baseline_v1.2), float32.
- `m_5_50.parquet`: per-column M_5_50.
- `w_ld.parquet`: SNP + weight LD score from weights.hm3_noMHC.@.
- The HapMap3 weights omit the MHC region (chr6: 75,359 baseline SNPs vs 72,387 weight SNPs). Because
  ldsc inner-joins the weights, the prepare Task restricts the baseline to SNPs with a weight LD score;
  `PreparedBaseline` then asserts that w_ld and the baseline LD scores share one SNP column, row for row.

### 4.3 Consumer-side alignment

- Inner-join sumstats onto the baseline SNP set by SNP (one join per trait, not per cell type), then join
  w_ld, then attach the annotation set and assert full coverage of the regression SNPs.
- Memory: the largest set (chromatin, ~490 unique columns × ~1.2M SNPs float32) is ~2.4 GB; fine for
  16 GB. The kernel streams cell types in chunks rather than materializing per-cell-type weighted copies.

## 5. Kernel: `mecfs_bio/build_system/task/sldsc/batched_sldsc_cts.py`

Pure numpy over plain arrays, no Task or file-format knowledge, in the style of
ppp_ldsc/batched_ldsc_h2.py. Entry point takes χ², N, w_ld, the shared (baseline) LD matrix and M, the
annotation LD matrix and M, the model table (as integer column indices), n_blocks, and
`weighting: SLDSCWeighting`. Returns a `@frozen` `BatchedSLDSCResult` (labels, coef, coef_se, p) with
shape/dtype invariants asserted in `__attrs_post_init__`.

`SLDSCWeighting = Literal["exact", "shared_controls"]` lives in `mecfs_bio/constants/`.

### 5.1 exact

- Once per trait: N̄, mean χ², baseline contributions to x_tot and M_tot, scaled shared columns.
- Per chunk of k cell types (k ≈ 32–64, chosen to bound memory): build the (k, S) weight matrix
  vectorized; per block, stack designs to (k, s_b, p) and form (k, p, p) Gram and (k, p) cross-products
  with batched `matmul`; one batched `np.linalg.solve` over (n_blocks, k, p, p) for delete values.
- Cost: ~1,015 × 1.2M × 56² × 2 ≈ 7.6e12 flops across all six sets, BLAS-bound — a few minutes.

### 5.2 shared_controls (fast, inexact)

- Group models by identical control set (all rows of an LDSC-SEG set form one group).
- Per group, one weight vector w_shared from the §3 formula with the target excluded: x_tot over
  baseline + controls, M_tot without M(t). That is, the weight of the model without the cell-type column.
- Shared design Z = [controls, baseline, 1] (scaled). Per block: G_b = Z_bᵀWZ_b and Z_bᵀWy once per group.
- Per cell type, per block: t_bᵀWZ_b (one (k × s_b)·(s_b × p) GEMM), t_bᵀWt_b, t_bᵀWy_b — O(S·K·p) total.
- Only the target coefficient is reported, and its SE uses only its own pseudovalues, so by Frisch–Waugh
  β_t = (tᵀWy − tᵀWZ G⁻¹ ZᵀWy) / (tᵀWt − tᵀWZ G⁻¹ ZᵀWt). Delete values factor G − G_j once per block
  (200 shared 55×55 factorizations) and apply scalar Schur complements for all K: ~6e8 flops.
- ≈ p-fold (~50×) fewer Gram flops than exact; runtime then dominated by parquet loading.
- Expected accuracy cost: target LD is a small fraction of x_tot, so weights move by a few percent;
  misspecified weights keep WLS unbiased and only cost efficiency. Guess ~1% of an SE; measured in §7.2.

## 6. Tasks, assets, generators

### 6.1 New package `mecfs_bio/build_system/task/sldsc/`

| Module | Contents |
|---|---|
| prepared_sldsc_reference.py | §4 contract: filename constants, loader/writer types. |
| prepare_ldcts_annotation_set_task.py | `PrepareLdctsAnnotationSetTask`: LDSC-SEG adapter. Fields: meta, source directory Task, `.ldcts` filename. Parses the index with the existing `LDCTSFileEntry` parser, reads each unique LD-score file once with polars, asserts identical SNP columns across files, writes the three parquets. `ReferenceDataDirectoryMeta`. |
| prepare_sldsc_baseline_task.py | `PrepareSLDSCBaselineTask`: baseline directory + `@` pattern, weights directory + `@` pattern → prepared baseline. |
| batched_sldsc_cts.py | §5 kernel. |
| vectorized_cell_analysis_by_sldsc_task.py | `VectorizedCellAnalysisBySLDSCTask`. Fields: meta, source_sumstats_task, baseline_task, annotation_set_task, weighting = "exact", pre_pipe = IdentityPipe(), min_hapmap3_snps = 100,000. execute() reads as a high-level spec: load sumstats + pipe → gwaslab-equivalent χ²/N prep → align to baseline → attach annotation set (assert coverage) → kernel → write CSV. `.create()` derives trait/project from the sumstats meta; output `ResultTableMeta` CSV with the same schema as `CellAnalysisByLDSCTask`. |

`LDCTSFileEntry` and its parser are reused in place from gwaslab_cell_analysis_by_sldsc.py (same
build_system layer, so no layering issue).

### 6.2 Asset instances

Under `mecfs_bio/assets/reference_data/linkage_disequilibrium_score_reference_data/prepared/`:

- six `PREPARED_*_ANNOTATION_SET` Tasks, one per existing extracted `_RD` Task (Cahoy, GTEx_brain, ImmGen,
  Corces_ATAC, Multi_tissue_chromatin, Multi_tissue_gene_expr);
- `PREPARED_BASELINE_V1_2_HM3_WEIGHTS` from `BASE_MODEL_PARTITIONED_LD_SCORES_EXTRACTED`
  (baseline_v1.2/baseline.@) and `PARTITIONED_MODEL_REGRESSION_WEIGHTS_EXTRACTED`
  (1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.@).

### 6.3 Generators

- `mecfs_bio/build_system/task_generator/sldsc_task_generator.py`: extract the downstream wiring
  (multiple testing → markdown → optional label join → scatter plot) into a free function accepting any
  cell-analysis Task that produces the result schema. Widen `CellAnalysisTaskGroup.cell_analysis_task` to
  `Task`. `SLDSCTaskGenerator` keeps its behaviour and asset ids via that helper.
- New `VectorizedSLDSCTaskGenerator.create(base_name, source_sumstats_task, baseline_task,
  entries: list[PreparedAnnotationSetRecord], weighting, multiple_testing_alpha, multiple_testing_method,
  pre_pipe)`. `PreparedAnnotationSetRecord` = entry_name + annotation_set_task + optional
  `CellOrTissueLabelRecord`. Source-agnostic.
- New `standard_vectorized_sldsc_task_generator(sumstats_task, base_name, pre_pipe, weighting="exact")` in
  `mecfs_bio/asset_generator/`: the fixed six-set list with the same entry names, label joins, and FDR
  settings as `standard_sldsc_task_generator`.
- Asset ids are new: every generated id carries a `_vec` marker (e.g.
  `{base_name}_{entry_name}_cell_analysis_by_s_ldsc_vec`), so old and new results coexist during comparison.

### 6.4 Migration

The old Task and generator remain: they are the validation oracle and back existing cached results.
Switching the 10 existing call sites is a follow-up after real-data validation.

## 7. Testing and validation

### 7.1 Unit tests (`test_mecfs_bio/unit/build_system/task/sldsc/`)

1. Synthetic LDSC-SEG fixture written by a helper into tmp_path: a few chromosomes, ~3,000 SNPs, ~3
   baseline annotations, a weights file, a 3-row `.ldcts` with one shared control, values printed with
   3 decimals, M_5_50 files, and sumstats with a planted enrichment in one cell type.
2. Exact-mode oracle (Task-level): prepare Tasks + `VectorizedCellAnalysisBySLDSCTask("exact")` vs the
   existing `CellAnalysisByLDSCTask` on the same fixture; coef, se, p agree to ~1e-6 relative (float32
   storage). Covers adapter, alignment, N handling, row order, weights and jackknife together.
3. shared_controls algebra (kernel-level; justified as the only independent oracle for the
   Frisch–Waugh/Schur path): vs a direct dense WLS + jackknife using the same shared weights and blocks,
   to floating-point tolerance.
4. Fail-fast: an annotation set not covering the regression SNPs raises; a model table referencing a
   missing column raises at construction; a regression SNP with missing N raises; fewer than
   min_hapmap3_snps HapMap3 SNPs raises.

### 7.2 Real-data experiment (`experiments/claude/vectorized_sldsc/`)

Reproducible script, output tee'd to an (uncommitted) log:

- Profile the old Task on one set (e.g. GTEx_brain) to confirm or refute the I/O-bound hypothesis.
- DecodeME, all six sets, three configurations: gwaslab (existing cached results), vectorized exact,
  vectorized shared_controls. Report max/median |Δz|, max |Δcoef|/se, FDR-call changes at α = 0.01,
  wall time, peak RSS.

## 8. Out of scope

- Switching existing call sites to the new generator.
- Adapters for other annotation sources (the contract leaves room for them).
- `print_all_cts` (reporting control-column coefficients).
- Making shared_controls the default.
