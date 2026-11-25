# Fast all-pairs CT-LDSC genetic correlation — Design

Date: 2026-10-25

## 1. Intent

`genetic_corr_by_ct_ldsc_asset_generator` (mecfs_bio/asset_generator/genetic_correlation_asset_generator.py)
drives the cross-trait genetic-correlation study in
mecfs_bio/assets/gwas/multi_trait/genetic_correlation/ct_ldsc/ct_ldsc_initial_asset_generator.py. It creates
one `GeneticCorrelationByCTLDSCTask` per trait pair plus one `SNPHeritabilityByLDSCTask` per trait. Each pair
task unpickles both gwaslab Sumstats, re-runs the same filtering (indels, palindromes, HLA, HapMap3,
duplicate rsIDs), re-reads the LD-score files, and has ldsc re-estimate both heritabilities. So every trait is
loaded and processed about T−1 times. With 13 traits (78 pairs) at about 20 s/pair, the study takes roughly
25 min. The user plans to add about 10 traits, giving 23 traits and 253 pairs, or about 85 min.

We will add a new asset generator that produces the same scientific outputs (an all-pairs rg table and a
per-trait heritability table) in **a few minutes total for all pairs**. It does each per-trait step once and
reuses the vectorized LDSC kernels we already validated for UKB-PPP
(mecfs_bio/build_system/task/ppp_ldsc/), so the repo does not gain another LDSC implementation.

## 2. Goals, non-goals, anti-goals

Goals
- All-pairs rg for about 23 traits in a few minutes, peak RSS well under 16 GB.
- Each trait's preparation is cached separately: adding one trait prepares only that trait, then reruns the
  cheap all-pairs step.
- Outputs readable unchanged by `GeneticCorrelationClustermapTask` (via `GeneticCorrSource` defaults
  p1/p2/rg/p) and by `LavaTask`'s sample-overlap reader (`lava_task.py` ~L404: p1, p2, gcov_int, h2_int in
  the pairs table; p, Intercept in the heritability table).
- Retain the heritability table, including liability-scale h² for binary traits.
- Accept two input formats: pickled gwaslab Sumstats (so the existing source list ports with no changes) and
  a plain dataframe asset (parquet etc.), which is the preferred format going forward.
- Reuse the ppp_ldsc helpers, extended where necessary, without changing PPP results.

Non-goals
- Reproducing gwaslab/ldsc output exactly. A fast, defensible CT-LDSC is enough.
- Replacing the existing generator. It stays, and `CT_LDSC_INITIAL_ASSET_GENERATOR`, the existing plot, and
  LAVA keep pointing at it. Any switch-over is a later decision, made after the comparison in §7.
- Moving or renaming the ppp_ldsc package to a generic location. That's a possible follow-up, not this
  change.

Anti-goals
- A second, divergent LDSC kernel. If something the kernels can't do comes up, extend the kernels; don't
  fork them.
- Silent behaviour changes for PPP. The PPP constant-N path must remain bit-identical.

## 3. Method (decided)

The estimator is the GenomicSEM-style estimator already implemented in batched_ldsc_rg.py / batched_ldsc_h2.py:
- free intercepts (not ldsc's two-step), chi-square filter max(0.001·N_max, 80) per trait, 200-block
  contiguous jackknife over the shared genome-sorted SNP set; missing or filtered SNPs are zero-weighted, not
  dropped;
- each trait's h² estimated on its own kept set, the covariance on the pair's kept intersection;
- rg SE = delete-block jackknife on the ratio (the classic ldsc rg SE).

Regression SNP set: HapMap3 (HAPMAP3_SNPLIST_FOR_GENOMIC_SEM, w_hm3.snplist: SNP, A1, A2) ∩ 1000G phase-3
EUR LD scores (THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE: CHR, SNP, BP, L2, M_5_50; the same
reference the gwaslab generator uses). Drop strand-ambiguous SNPs (w_hm3 is already palindrome-free, so
this is a defensive no-op) and the extended MHC in build 19 (positions come from the LD file's BP, which
is hg19). w_hm3 has no indels. M = total M_5_50 of the LD reference, taken before any join (as
build_ppp_ldsc_context already does).

Trait alignment: join on rsID to the index, keep z = BETA/SE where alleles match the index (A1 = effect
allele), negate when swapped, NaN otherwise; drop duplicate rsIDs deterministically. This is exactly
`align_trait_to_context`. Non-finite z (SE = 0) falls out through the kernels' isfinite masks, so no separate
degenerate-z step is needed.

Liability-scale h² (heritability table only; rg is scale-free): the existing
`liability_h2(h2_obs, phenotype_info)` in gwaslab/ldsc_diagnostic_plot_task.py (gwaslab's h2_obs_to_liab with
P = sample prevalence, K = population prevalence). The SE is scaled by the same factor (the conversion is
linear in h2_obs).

## 4. Task structure (decided)

```
LD scores (consolidated) ─┐
HapMap3 w_hm3.snplist ────┴─► CtLdscVariantIndexTask (one, shared)
                                   │
          SumstatsSource_t ────────┼─► CtLdscAlignedTraitTask_t  (one per trait)
                                   │          │
                                   ▼          ▼
                       CtLdscAllPairsRgTask      CtLdscHeritabilityTask
                                   │
                        ConvertDataFrameToMarkdownTask
```

### 4.1 CtLdscVariantIndexTask (shared)
Deps: HapMap3 snplist task, consolidated LD-score task. Output: a parquet table with columns CHR, POS
(hg19, from LD BP), rsID, EA (= w_hm3 A1), NEA (= A2), PPP_INDEX_IS_STRAND_AMBIGUOUS_COL (via
is_palindromic_expr from harmonize_gwas_with_reference_table_via_rsid.py, as the PPP index does). It's a reference-style asset, independent of any trait, built once per
(snplist, LD reference) pair.

Downstream tasks get the regression context by calling
`build_ppp_ldsc_context(index_df, ld_df, drop_strand_ambiguous=True, exclude_mhc=True, build="19")`.
Building the context takes seconds, and the result is a deterministic function of the index and LD
assets, so every consumer derives the same row set and order.

### 4.2 CtLdscAlignedTraitTask (one per trait)
Deps: the source's sumstats task, the index task, the LD-score task. Inputs: a `SumstatsSource` (task,
alias, sample_info, pipe) from gwaslab_genetic_corr_by_ct_ldsc_task.py, with one new optional field (default
None, so existing uses are unaffected): `alignment_guards: AlignmentGuards | None` (see 4.2.1).

Execute:
1. Read the trait as a narwhals LazyFrame, choosing by the source task's meta:
   - `GWASLabSumStatsMeta`: `read_sumstats(asset).data` → polars (cast categorical EA/NEA/rsID to String;
     gwaslab stores alleles as pandas Categoricals) → lazy narwhals;
   - a meta whose read_spec is a `DataFrameReadSpec`: `scan_dataframe_asset(..., parquet_backend="polars")`
     (lazy, so the rsID filter is pushed into the scan).
   - any other meta: assert with a clear message (fail at `.create()` time).
2. Apply `source.pipe.process(...)`.
3. `align_trait_to_context(trait_lazy, context_variants, trait_total_sample_size=
   source.sample_info.total_sample_size, ...)`, which now also returns the alignment counts (see §5).
4. Compute the alignment diagnostics and enforce the guards (4.2.1). Fail fast before writing anything.
5. Write a DirectoryAsset holding aligned.parquet (columns rsID, z, n in **context row order**, S rows)
   and diagnostics.json (the 4.2.1 diagnostics).

#### 4.2.1 Upstream-bug guards (decided: fail fast, per-source override)
These protect against upstream harmonization bugs, chiefly rsIDs assigned from positions in the wrong
genome build or otherwise wrongly assigned, and sumstats labelled with the wrong build. The join is by rsID,
so a correctly rsID'd trait in any build is fine; the guards target wrong rsIDs and mislabelled builds.

Diagnostics (all computed after the pipe, over the context SNP set):
- n_context: S. n_found: context rsIDs present in the trait. n_matched: found and alleles match the index
  directly or swapped. n_allele_mismatch: n_found − n_matched.
- coverage = n_matched / n_context.
- allele_mismatch_fraction = n_allele_mismatch / n_found. Correct rsIDs mismatch well under 1%; rsIDs that
  point at the wrong variant mismatch at a large rate, so this is the sharpest signal.
- position_concordance (only when the trait has CHR and POS columns): among matched SNPs, the fraction whose
  (CHR, POS) equals the index's hg19 (CHR, POS). Null when positions are absent.

Guards (asserts with a message naming the alias, the diagnostic value and the threshold):
- coverage >= min_coverage.
- allele_mismatch_fraction <= max_allele_mismatch_fraction.
- position_concordance >= min_position_concordance, applied only when positions exist AND the source
  meta declares build 19 (the old gwaslab path asserted infer_build() == "19"; this replaces that). Skip for
  other or unknown builds. Positions are only a diagnostic for those builds.

`AlignmentGuards` is a frozen attrs class holding min_coverage, max_allele_mismatch_fraction and
min_position_concordance, each validated to lie in [0, 1] in __attrs_post_init__. The generator takes a
default AlignmentGuards; a source's `alignment_guards`, when set, replaces it for that source (e.g. a
pre-restricted or unusual input). Provisional defaults: 0.5 / 0.02 / 0.95. The **first step of the §7
comparison experiment** measures the diagnostics on the current 13 traits, and the builder sets the defaults
from those numbers (with comfortable margin) before wiring the asset. This replaces PPP's absolute
min_trait_snps floor for this feature.

Meta: derived from the source task's meta (trait/project), per repo convention; the asset id includes the
generator base_name and the alias, because the same sumstats task can appear in two generators with
different pipes.

### 4.3 CtLdscAllPairsRgTask (one)
Deps: index task, LD-score task, all aligned-trait tasks. Fields: the ordered tuple of (alias, aligned task,
sample_info), `n_blocks` (200), `batch_size`.

Execute (should read as a high-level spec; loops live in helpers):
1. Build the context; read every aligned trait. Assert each one's rsID column equals the context's rsIDs
   (guards against a context/prep mismatch). Stack into Z (S, T) and N (S, T).
2. For each trait, `estimate_trait_context(Z[:, t], N[:, t], ld, m, n_blocks=...)`.
3. Traits are ordered by alias.lower(), as the existing generator does. For each i, compute the
   covariances against traits j > i in column batches of `batch_size` with `batched_rg(ctx_i, Z[:, J],
   ld, N[:, J], m, ...)`.
4. Emit one row per unordered pair (p1 = trait i, p2 = trait j).

Pairs-table columns: p1, p2, rg, se, z, p, gcov (observed scale), gcov_se, gcov_int, gcov_int_se, h2_obs_p1,
h2_obs_p2, h2_int (the intercept of p2's own h² regression, which LAVA's diagonal fallback reads),
n_snps (SNPs in the covariance). Format: parquet via ResultTableMeta (trait multi_trait, project
genetic_correlation, sub_dir e.g. ct_ldsc_fast).

### 4.4 CtLdscHeritabilityTask (one)
Same deps and steps 1–2 as 4.3. One row per trait. Columns: p (alias), h2_obs, h2_obs_se, h2_liab, h2_liab_se
(null for quantitative traits), Intercept, Intercept_se, lambda_gc, mean_chi2, n_snps, n_bar, plus the 4.2.1 diagnostics read from each aligned
trait's diagnostics.json (coverage, allele_mismatch_fraction, position_concordance). These let the user spot
a trait that passed the guards but looks off (mean_chi2 near 1, h2 <= 0). It's a separate
Task because a Task produces one asset. Recomputing `estimate_trait_context` costs seconds and is
deterministic, so it agrees with the h² values in the pairs table.

### 4.5 Generator
New file mecfs_bio/asset_generator/fast_genetic_correlation_asset_generator.py:
`fast_genetic_corr_by_ct_ldsc_asset_generator(base_name, sources, ld_scores_task=..._CONSOLIDATE,
hapmap3_task=HAPMAP3_SNPLIST_FOR_GENOMIC_SEM, ...)` returns a frozen attrs object holding: index_task,
aligned_trait_tasks (alias → task), pairs_task, heritability_task, pairs_markdown_task, and `terminal_tasks()`.
The index task may be built inside the generator or be a module-level asset passed as a default (builder's
call).

### 4.6 Assets
- Move the `sources=[...]` list out of ct_ldsc_initial_asset_generator.py into a shared constant (e.g.
  `CT_LDSC_INITIAL_SOURCES`) so that both generators consume the identical list.
- New ct_ldsc_fast_initial_asset_generator.py: `CT_LDSC_FAST_INITIAL_ASSET_GENERATOR =
  fast_genetic_corr_by_ct_ldsc_asset_generator("initial_rg_fast", sources=CT_LDSC_INITIAL_SOURCES)`.
- New clustermap asset for the fast pairs table (GeneticCorrSource with defaults), for side-by-side viewing.

## 5. Kernel reuse and the minimal extensions required

The following are reused as-is: `align_trait_to_context`, `TraitAligned`, `PppLdscContext`, `_blocked_wls`'s
block machinery, `estimate_trait_context`'s h² path, and `batched_rg`'s covariance/rg/jackknife path.

Required extensions (all backward compatible; the existing PPP tests must pass unchanged):

0. **Alignment counts from `align_trait_to_context`.** It currently drops allele-mismatched rows silently.
   Have it also return n_found and n_matched (as fields on the attrs result, alongside z and n). If the
   trait has CHR/POS, it should return the matched SNPs' trait positions too, or a position-concordance
   count against index columns passed in, so the task can compute position_concordance. The PPP caller and
   its min_trait_snps assert stay unchanged.
1. **`build_ppp_ldsc_context` build parameter.** MHC exclusion currently hardcodes build "38". Add
   `build: GenomeBuild = "38"` and use `extended_mhc_interval(build)`.
2. **Per-SNP N for the second side of `batched_rg`.** `n_protein` is currently (K,). Also accept (S, K)
   (NaN where absent), with the chi-square threshold from each column's nanmax N, exactly as
   `estimate_trait_context` does for the first side. With that change, the h2_p that `batched_rg`
   computes for trait j equals `estimate_trait_context(trait j).h2`, so diagonal h² stays consistent
   across pairs. Several of the current sources carry a real per-SNP N column (e.g. Fibromyalgia), so
   this is needed. The (K,) path must stay bit-identical.
3. **Intercept jackknife.** `_blocked_wls` currently returns only slope deletes. Also return
   intercept deletes (as an attrs result rather than a growing tuple), so gcov_int_se and Intercept_se
   can be computed.
4. **Heritability diagnostics on `TraitLdscContext`.** `estimate_trait_context` discards the intercept.
   Add h2_se (jackknife over h2_delete; algebraically the same as batched_h2's pseudo-value form),
   intercept, intercept_se, mean_chi2, lambda_gc, n_snps. Add gcov_se and gcov_intercept_se to
   `BatchedRgResult`.

`batched_h2` is not used for the heritability table because it assumes constant N per column. Item 4 makes
`estimate_trait_context` the single h² path for this feature.

## 6. Constraints
- Run everything via pixi; `pixi r invoke green` after changes (capture the log to a file).
- Peak memory: Z and N for 23 traits × ~1.1M SNPs ≈ 0.4 GB. `batched_rg` allocates on the order of 10–15
  (S, K) temporaries, so `batch_size` bounds peak memory (default 8, about 1 GB of temporaries). There is
  no need to bound how many traits are loaded.
- Repo conventions: imports at top of file, no monkeypatching, Literal types for enumerable strings,
  Task-level unit tests, docstrings without backticks or RST, column-name constants rather than repeated
  literals (put new ones in a constants module).
- Guard thresholds: see 4.2.1. Syncope's input is already HapMap3-restricted, which should not hurt
  coverage, since the index is itself HapMap3.

## 7. Validation and tests
- Unit (kernel): per-SNP-N `batched_rg` reproduces `estimate_trait_context` h² for the second side; the
  constant-N path is unchanged (the existing PPP tests); intercept SE matches `exact_rg_single` /
  genomic_sem_ldsc estimates on synthetic data to the same tolerance as the existing PPP probes.
- Unit (Task level): synthetic traits with known rg (see the CT-LDSC synthetic recipe memory): the
  aligned-trait task handles both input formats and allele flips; each guard trips on a synthetic trait
  built to violate it (rsIDs shuffled across variants → allele mismatch; most rows dropped → coverage;
  positions shifted → concordance) and a per-source override lets that trait through (no match= on the
  error text); the all-pairs task outputs every
  unordered pair once, with rg close to truth; feeding the pairs table to `GeneticCorrelationClustermapTask`
  and the overlap reader in lava_task runs.
- Experiment (not CI): experiments/claude/fast_ct_ldsc_vs_gwaslab/ first tabulates the 4.2.1
  diagnostics for every current trait (to set the guard defaults), then runs both generators on the current
  13 traits, records wall time, and tabulates per pair Δrg, Δrg/se and Δh². Expected: agreement within a
  fraction of an SE (the estimators differ in intercept handling and chi-square filtering), and no
  change in which pairs pass Bonferroni. Report the results to the user; don't gate on exact agreement.

## 8. Calls made without explicit user sign-off (check these)
- Two output Tasks (pairs, heritability) that each recompute the per-trait h², rather than one Task that
  writes a directory.
- Parquet for both tables (the gwaslab generator writes CSV). Consumers read through scan_dataframe_asset,
  so the format doesn't affect them.
- Duplicate rsIDs within a trait: keep one row deterministically (align_trait_to_context's existing
  behaviour). The gwaslab path drops all copies.
- Liability SE = observed SE × the conversion factor.
- The aligned-trait asset is a directory (aligned.parquet + diagnostics.json) so that diagnostics travel
  with the trait to the heritability table.
- Guard defaults are provisional until measured on the current traits.
- The source list moves to a shared constant, which touches ct_ldsc_initial_asset_generator.py (no change
  to its asset ids).

## 9. Left to the builder
- Exact module and class names, and the package location of the new Tasks (suggested:
  mecfs_bio/build_system/task/ct_ldsc_fast/).
- Exact heritability/pairs column constant names beyond the ones consumers require (p1, p2, rg, p,
  gcov_int, h2_int, p, Intercept).
- Whether the index task is a module-level asset or built by the generator.
- Default `batch_size`, chosen from a measured peak-RSS check.
