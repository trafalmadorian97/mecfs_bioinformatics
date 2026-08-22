# CSF protein heritability task — plan

Per-aptamer LDSC SNP-heritability over the HapMap3 CSF pQTL database (Western et al.
2024), modeled on `PppProteinHeritabilityTask`. All-variants h² only (cis-excluded
path deferred).

## Settled decisions

- **Sample size: median N per aptamer.** The `n_spread_probe` showed within-aptamer N
  is effectively constant over the HapMap3 regression SNPs (IQR/median = 0; median→per-SNP
  h² shift p95 ≈ 0.01%, max 0.03%). So one scalar N per aptamer, taken as the **median**
  of the per-variant N over the present context SNPs (not `constant_sample_size`, which
  asserts equality — CSF has a thin low-N tail that would trip it).
- **No cis-excluded row** (skip gene-coords dependency entirely for now).
- **Output identity: analyte + UniProt + entrez gene symbol.**
- **Reuse the same LD-score reference** (`ConsolidateLDScoresTask` →
  `THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE`).
- **Shared row order confirmed by assertion** in `create()` (holds by construction from
  `BuildSlimCsfAptamerParquetTask`, which aligns every aptamer onto the one index).

## Reused unchanged (no PPP edits)

- `PppLdscContext` + `build_ppp_ldsc_context` — the CSF index already exposes
  CHR/POS/rsID/`is_strand_ambiguous` (same column string), so the context builder works
  as-is. Import and reuse; do not fork.
- `batched_h2` — takes a scalar N per aptamer; reuse directly (all-variants, no `exclude`).

## New code

1. **`CsfAptamerFile` gains `uniprot`** (`mecfs_bio/build_system/task/csf_database/
   build_slim_aptamer_parquet_task.py`) — populated in `csf_slim_aptamer_asset_generator`
   from `UNIPROT_MANIFEST_COLUMN`. Add a `UniProtId` NewType to `csf_database_constants`
   for the signature. (CSF-only change; the slim parquet output is unaffected.)

2. **`mecfs_bio/constants/csf_ldsc_constants.py`** — output column-name constants:
   analyte, uniprot, gene_symbol, variant_set, h2, h2_se, intercept, mean_chi2,
   lambda_gc, n_snps, n_bar. The `variant_set` column is constant "all_variants" today
   (via `CsfVariantSet` literal) but kept for forward compatibility so a future
   cis-excluded set can be added as extra rows without a schema change.

3. **`mecfs_bio/build_system/task/csf_ldsc/csf_protein_heritability_task.py`** —
   `CsfProteinHeritabilityTask(GeneratingTask)`:
   - `deps` = aptamer_tasks + index_task + ld_scores_task (no sample-size / gene-coords task).
   - `execute`: build context once; batch aptamers (`batch_size` from config); per aptamer
     read BETA/SE/N at `context.row_pos`, chi² = (β/se)², `n = median_sample_size(N_at_context)`;
     `batched_h2(chi2, ld, n_arr, m)`; one output row per aptamer.
   - `median_sample_size(n_at_context, label)` free function — median of finite values;
     assert at least one finite (clear message naming the aptamer).
   - `create(...)`: assert `ld_scores_task` produces a dataframe; assert every aptamer task's
     `index_task.asset_id == index_task.asset_id` (shared-row-order guard); build a
     `ResultTableMeta` (trait `western_csf`, project `csf_heritability`).

4. **`mecfs_bio/assets/gwas/.../hapmap3/hapmap3_csf_heritability.py`** — wire
   `HAPMAP_3_CSF_DATABASE.aptamer_tasks` + `HAPMAP_3_CSF_DATABASE_INDEX` +
   the consolidated LD reference into `CsfProteinHeritabilityTask.create(...)`.

## Tests

- `median_sample_size`: median of a NaN-mixed array; low-N tail doesn't trip it; all-NaN
  raises. (Mirrors PPP's `constant_sample_size` tests.)
- Small Task-level smoke test with synthetic index + 2–3 tiny aligned aptamer parquets and
  a toy LD-score table, asserting one row per aptamer with the expected identity columns.
- `batched_h2` / `build_ppp_ldsc_context` already validated — not re-tested.

## Out of scope (later)

Cis-excluded h² (UniProt→hg38 coords), cross-trait rg, LCV, an analysis/build entry point.
