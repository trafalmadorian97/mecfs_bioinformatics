# PolyFun Approach 2 — trait-specific L2-regularized S-LDSC prior

Date: 2026-09-20. Status: design approved (brainstorming), ready for
implementation plan.

## Context

We fine-map GWAS loci with SUSIE, optionally supplying a PolyFun per-SNP
heritability prior. Today only PolyFun **Approach 1** is wired: a *precomputed*
`snpvar_bin` (`COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS`) learned from a
meta-analysis of 15 UK Biobank traits, floored via `create_prior_col_pipe(q)` and
handed to `SusieRFinemapTask` as a `PriorInfo`
(`mecfs_bio/asset_generator/polyfun_explain_fine_mapping_asset_generator.py`).

**Approach 2** instead estimates the per-SNP heritability (`snpvar`) *from the
GWAS being fine-mapped*, via an L2-regularized extension of S-LDSC over the
baseline-LF annotations. This spec adds Approach 2 as a drop-in alternative prior
so we can compare Approach-1 vs Approach-2 priors on the same loci (the immediate
DecodeME question: does DecodeME's own GWAS want to materially reorder the
generic functional prior?).

Two motivations (from the brainstorm):
1. **Prior-mismatch check.** We don't know whether the generic 15-trait prior
   fits our traits; running Approach 2 and diffing credible sets against
   Approach 1 is a direct way to find out.
2. **Path to custom annotations.** The longer-term goal is non-baseline-LF
   annotations (e.g. AlphaGenome). That needs Approach 3 (bin LD-scores from a
   reference panel), designed separately. Approach 2 is the standalone first step
   and establishes the estimator + wiring that Approach 3 will extend.

### The estimator, in one paragraph

Regress GWAS chi-square on the baseline-LF **annotation LD-scores** with LDSC
weights and an L2 (ridge) penalty; the fitted annotation coefficients `tau_c`
give each variant a predicted per-SNP heritability `snpvar_i = sum_c a_ic tau_c`.
To avoid winner's curse, `tau` is fit on even chromosomes to score odd
chromosomes and vice versa (`num_chr_sets = 2`) — and the ridge lambda (and the
weights) for each parity's fit are selected using only that parity's chromosomes,
so nothing that scores a chromosome is fit on data from its own parity (decision
8). The result is a per-variant
`snpvar` table that plugs into `PriorInfo` exactly where the Approach-1 prior does
(the floor/constrain is applied at the `PriorInfo` seam via the existing
`create_prior_col_pipe(q)`, not baked into `snpvar`).

### Decision: pragmatic weighted ridge, not a faithful PolyFun port (approved)

We deliberately do **not** chase bit-for-bit parity with PolyFun's
`ldsc_polyfun` `Hsq` internals. We keep the parts that matter for a *prior* —
even/odd LOCO for winner's-curse protection, ridge lambda by our own
leave-one-chromosome-out CV, per-column standardization, LDSC `omega = het * oc`
weighting — and we lean on the fact that **SUSIE normalizes the prior within each
locus** (`w_i = snpvar_i / sum_j snpvar_j`), so only *relative* `snpvar` within a
fine-mapping window matters, not absolute calibration. This mirrors the repo's
existing "intentionally NOT a faithful GenomicSEM port" precedent
(`genomic_sem_gwas_by_subtraction_full_python_task.py`). Validation is therefore
against *synthetic data with planted enrichment*, not against a PolyFun run.

### Decision: regress on the dense variant set, streamed by chromosome (approved)

baseline-LF is a *low-frequency* annotation model; restricting the regression to
HapMap3 would discard exactly the low-frequency regression SNPs it targets. The
dense LD-scores are ~29GB (~19.5M x 187), but — as in `RidgeAnnotationWeightsTask`
— the full design matrix is never materialized: the ridge normal equations need
only the weighted cross-products `X^T W X` (187x187) and `X^T W y` (187), which
accumulate additively across chromosomes. We stream one chromosome's
`.l2.ldscore` member at a time, fold it into per-chromosome sufficient-statistic
blocks, and discard the rows. A HapMap3 restriction remains available as an
*injected, optional* SNP-list filter (fast iteration), but dense is the default.

### Decision: factor a pure-numpy chromosome-blocked ridge submodule (approved)

The LOCO combine / standardize-solve-unstandardize / held-out-R^2 / lambda-by-LOCO
math is identical to what `RidgeAnnotationWeightsTask` already hand-rolls
(`_ChromStats`, `_standardized_system`, `_solve`, `_heldout_r2`,
`_select_alpha_loco`, `_combine`). Weighted is the general case; that task's
current code is the `w == 1` special case. We factor the linear algebra into one
**pure-numpy submodule with zero IO/Task/domain knowledge** (no parquet, no
annotations, no chi-square). Each task does its own scan + join + (new task only)
weight computation, hands the submodule per-chromosome `(X, y, w)` arrays, and
gets back blocks and fitted coefficients. Sequenced to de-risk (see Phasing): the
new estimator lands on the submodule first; migrating `RidgeAnnotationWeightsTask`
onto it is a separate change gated by a numeric-equivalence test.

## Definition of done

- A new task produces a durable per-variant `snpvar` table
  (`chrom, pos, nea, ea, snpvar`) from a munged GWAS + the baseline-LF annotation
  LD-scores, green under `pixi r invoke green`.
- On synthetic data with planted annotation enrichment, the recovered `snpvar`
  ranking recovers the planted ordering (Spearman near 1), and every input to
  scoring a chromosome's parity is fit only on the opposite parity (weights,
  lambda, and coefficients — see decision 8).
- The fine-mapping generator accepts an injectable prior source; passing the
  Approach-2 asset yields a matched SUSIE run whose credible sets can be diffed
  against the Approach-1 run on the same locus.
- The chromosome-blocked ridge submodule has its own focused unit tests.

## Design decisions (approved)

1. **Fidelity:** pragmatic weighted ridge; validate on synthetic, not vs PolyFun
   (see Context).
2. **Regression SNP set:** dense by default, streamed by chromosome; HapMap3 an
   optional injected filter.
3. **Shared submodule:** pure-numpy chromosome-blocked ridge, used by the new
   task now and `RidgeAnnotationWeightsTask` after a gated migration.
4. **Weights source — derive, don't download.** LDSC weights are
   `omega_i = 1 / (het_i * oc_i)` with `het_i = (1 + N * h2bar * l_i / M)^2` and
   `oc_i = max(l_i, 1)`, where `l_i` is the variant's total LD-score. The
   distributed baseline-LF LD-score members have **no `base` (all-ones)
   column** (schema spike, 2026-09-21), but the 20 `MAFbin_*` annotations
   partition every SNP (each SNP is in exactly one MAF bin, so their annotation
   values sum to 1 per SNP). Therefore `l_i = sum over the 20 MAFbin_* LD-score
   columns` is the exact total LD-score, and `M` = the reference SNP count
   (number of rows in the LD-score members / annotation matrix). This uses the
   LD-score members only (no `base` column, no separate `--w-ld-chr` weights
   download, no `.l2.M` file). Using `l_i` as the over-counting weight (rather
   than a dedicated regression-SNP weight file) is an acceptable pragmatic
   simplification consistent with decision 1 and with PolyFun's dense UKB setup,
   which uses essentially the same SNP set for LD-scores and regression. `h2bar`
   is a scalar from a first unweighted streaming pass — and is itself
   **per-parity** (`h2bar_odd`, `h2bar_even`), so the weights of the odd-fit model
   use only odd-chromosome data (decision 8). `M` and `l_i` carry no GWAS signal,
   so they stay global.
5. **Floor at the seam, not in `snpvar`.** The output is raw `snpvar`; the
   `max/q` floor is applied by reusing `create_prior_col_pipe(q)` inside the
   `PriorInfo`, identical to Approach 1. Keeps the two approaches swap-compatible.
6. **`M` and total-LD-score self-consistency.** `l_i` (MAFbin-sum) comes from
   the LD-score members; `M` (reference SNP count) is the row count of that same
   set. The LD-score column names are identical to the annotation-matrix column
   names (schema spike), so `tau_c` (fit on LD-scores) applies directly to `a_ic`
   (annotation matrix) by name. Assert the two column-name sets match at build.
7. **Exact `(chrom, pos, nea, ea)` join key — non-negotiable.** Every join
   between per-variant frames (munged sumstats, annotation LD-scores, annotation
   matrix, the snpvar output, the `PriorInfo` join in SUSIE) is on the exact
   four-part key `(chrom, pos, nea, ea)`. No rsid join, no unordered-allele key,
   no position-only join — any of those risk silent information loss. This aligns
   with the in-flight exact-join direction ([[project_unordered_allele_key_cleanup]],
   REF-oriented per [[project_baselinelf_and_broad_ld_are_ref_oriented]]). **If
   any required input lacks these four columns (e.g. the LD-score members turn out
   to carry only `SNP`/rsid), STOP and ask before proceeding** — do not fall back
   to another key.
8. **Ridge lambda AND weights nested within each parity fold (no leakage).** To
   score the even chromosomes we fit a model *entirely on odd chromosomes*: the
   weights use `h2bar_odd`, and lambda is chosen by leave-one-odd-chromosome-out
   CV over the **odd** chromosomes only, then the model is refit on all odd
   chromosomes at that lambda to produce `tau_odd` (which scores even
   chromosomes). The even-fit model is symmetric and scores odd chromosomes. So
   `select_alpha_loco` runs twice — once on the odd block dict, once on the even —
   never once globally. This is stricter than PolyFun (which selects lambda
   globally) and guarantees no aspect of the even-scoring model is fit on even
   data.

## Components

### A. Annotation LD-score members — new stream-extract (reference data)

- The `.l2.ldscore.parquet` members live in the **same** 30GB
  `baselineLF_v2.2.UKB.polyfun.tar.gz` bundle
  (`s3://broad-alkesgroup-ukbb-ld/UKBB_LD/...`) that
  `StreamExtractAnnotationParquetsTask` already streams to extract the
  `.annot.parquet` members — the LD-score members are currently *discarded*.
- Add a sibling stream-extract task that retains the
  `baselineLF2.2.UKB.<chr>.l2.ldscore.parquet` members instead. Prefer
  generalizing the existing task with an injected member-name predicate over
  copy-paste (both instantiations then differ only by the predicate + output
  meta). Output: `DirectoryAsset` of per-chromosome LD-score parquets;
  `ReferenceDataDirectoryMeta(group="polyfun", sub_group="annotations",
  sub_folder="raw", id="baseline_lf_2.2_ukb_ldscore_parquet_members")`. A
  `path_remap` candidate for `/mnt/d` (large, few-file, rarely read).
- Cost: re-streams the 30GB bundle once (the existing task streams-and-discards
  rather than storing the tarball at rest, so there is no cached copy to reuse).
  Acceptable one-time-per-machine cost; the alternative (store the 30GB tarball
  once and extract both member kinds) is noted but rejected to keep the existing
  annotation task untouched and avoid 30GB at rest.
- **Schema (RESOLVED by spike 2026-09-21,
  `experiments/claude/polyfun_approach2_ldscore_schema_spike/`):** each member has
  `CHR:int64, SNP:string, BP:int64, A1:string, A2:string` + the **187 annotation
  LD-score columns named identically to the annotation matrix** (192 columns
  total). So the four-part key is available (`CHR,BP` + `A1/A2`; baseline-LF is
  REF-oriented so `A1=REF=nea`, `A2=ALT=ea`). There is **no `base` column** —
  handled by the MAFbin-sum total LD-score (decision 4). No stop-and-ask needed.

### B. Chromosome-blocked ridge submodule (new, pure numpy)

- Location: `mecfs_bio/build_system/task/annotation_weights/`
  (e.g. `chromosome_blocked_ridge.py`) — colocated with both consumers; pure
  numeric, could move to `util/` only if a third consumer outside this area
  appears. No polars/parquet/Task imports.
- Public surface (weighted; `w=None` recovers the unweighted case):
  - `@frozen ChromRidgeBlock`: weighted cross-product sufficient stats over `p`
    features — `sw, swx (p), swxx (p,p), swxy (p), swy, swyy`. `attrs` post-init
    asserts shapes/dtype per repo convention.
  - `accumulate_block(x, y, w=None) -> ChromRidgeBlock`
  - `combine(blocks) -> ChromRidgeBlock` (sum subsets: LOCO / even-odd)
  - `standardized_system(block) -> StandardizedSystem` (center + per-column
    L2/std, weighted)
  - `solve(system, alpha) -> beta_std`
  - `heldout_r2(held, beta_std, train_moments) -> float`
  - `select_alpha_loco(blocks_by_chrom, alphas) -> AlphaSelection(alpha,
    mean_r2, r2_per_chrom)` (named attrs, not a bare tuple)
  - `fit(block, alpha) -> RidgeFit(beta_raw, beta_std, intercept)`
- Generalizes `RidgeAnnotationWeightsTask`'s existing private helpers; the derived
  weighted moment formulas extend the unweighted `_standardized_system` /
  `_heldout_r2` derivations already documented in that file.

### C. L2-regularized S-LDSC snpvar estimator (new Task)

- Location: `mecfs_bio/build_system/task/annotation_weights/`
  (e.g. `l2_sldsc_snpvar_task.py`). Working class name
  `L2RegularizedSldscSnpvarTask`.
- Injected deps (all Tasks; named kwargs):
  `munged_sumstats_task`, `annotation_ldscore_members_task` (Component A),
  `annotation_matrix_task` (default `BASELINE_LF_ANNOTATION_MATRIX`),
  optional `regression_snp_restriction_task` (default `None` = dense).
  The annotation matrix is a dep both to derive `M` / total LD-score and to score
  `snpvar` densely.
- `execute` (uniformly per-chromosome streaming; `p = 187`; all joins on
  `(chrom, pos, nea, ea)` per decision 7):
  1. **Load munged chi-square + N** joined to the LD-score members on
     `(chrom, pos, nea, ea)` (chi-square = Z^2; filter `chi2 < 80`, PolyFun's
     `MAX_CHI2`). Reuse the repo's degenerate-Z guard before squaring. Assert
     every frame carries the four-part key before any join.
  2. **Derive `M` and per-variant total LD-score `l_i`** from the LD-score
     members: `l_i = sum of the 20 MAFbin_* LD-score columns`, `M` = reference
     SNP count (decisions 4, 6) — global (no parity split; no GWAS signal).
  3. **Pass A (per-parity scalar h2bar):** stream chromosomes accumulating
     *unweighted* blocks (`w=None`) per chromosome; `h2bar_odd` from
     `combine(odd unweighted blocks)`, `h2bar_even` from the even ones (rough fit
     at a nominal alpha). Only the two scalars are retained.
  4. **Pass B (weighted per-chromosome blocks):** stream again; each chromosome is
     weighted with **its own parity's** h2bar — `omega_i = 1/(het_i * oc_i)`,
     `het_i = (1 + N * h2bar_parity(chrom_i) * l_i / M)^2`, `oc_i = max(l_i, 1)`
     (`omega = het * oc` per the repo's batched-LDSC note). Accumulate one weighted
     `ChromRidgeBlock` per chromosome; keep the odd and even block dicts separate.
  5. **Per-parity lambda + tau (nested, decision 8):**
     `sel_odd = select_alpha_loco(odd_blocks)`, refit
     `tau_odd = fit(combine(odd_blocks), sel_odd.alpha).beta_raw / Nbar`; symmetric
     `tau_even` from the even blocks. `select_alpha_loco` is thus called once per
     parity, its LOCO folds drawn only from that parity's chromosomes.
  6. **Score dense:** stream the annotation matrix per chromosome; an even
     chromosome is scored with `tau_odd`, an odd chromosome with `tau_even`
     (`snpvar_i = a_i . tau[opposite_parity(chrom_i)]`); write `snpvar.parquet`
     (`chrom, pos, nea, ea, snpvar`) via
     `write_df_according_to_format(..., ParquetOutFormat())`.
  7. **Diagnostics json:** chosen `alpha` per parity, per-chrom held-out R^2,
     total h2, `Nbar`, and the per-annotation `tau` table (odd + even) — feeds
     later comparison/explainability and is cheap insurance.
- Output: `DirectoryAsset` (snpvar parquet + diagnostics json),
  `ResultDirectoryMeta` derived in `create()` from the munged sumstats task's
  meta (trait/project pulled from the dep, per repo convention — never accepted as
  args). Filename constants are module-level for tests/consumers.

### D. Wire Approach 2 into the fine-mapping generator

- Today `generate_polyfun_explain_group` hardcodes the polyfun arm's `PriorInfo`
  to `COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS` +
  `create_prior_col_pipe(q)` + `POLYFUN_PRIOR_COL`.
- Introduce a small `@frozen` prior-source value object
  (`prior_task`, `prior_col`, `prior_pipe`) with a default constructed from the
  Approach-1 constants, and thread it through
  `generate_assets_polyfun_explain_fine_map` /
  `generate_polyfun_explain_group` (default preserves current behaviour exactly).
- An Approach-2 prior source points `prior_task` at the Component-C snpvar asset
  for that trait, `prior_col = "snpvar"`, `prior_pipe = create_prior_col_pipe(q)`.
  Because `snpvar` is dense over the annotation variants, the `PriorInfo` join
  cannot hard-error on missing coverage (asserted by a coverage test).
- The Approach-1-vs-Approach-2 comparison is then a matched pair of generator
  invocations for the same locus differing only in the prior source; credible-set
  agreement is read out with the existing UpSet-over-runs tooling. Explaining
  Approach 2 in annotation terms *via its own `tau`* (rather than the precomputed
  ridge weights) is a natural extension but **out of scope here**.

## Memory / compute notes

- Estimator peak ~ one chromosome of LD-scores (~1.5M x 187 f32 ~= 1.1GB) + the
  22 per-chromosome `p x p` f64 blocks (~6MB) + munged chi-square frame. Three
  streaming passes (unweighted, weighted, score) over per-chromosome parquet, one
  chromosome resident at a time. Comfortably < 16GB.
- Single-machine, no Docker, default pixi env. The ~29GB LD-score members
  directory is `path_remap` to `/mnt/d`, like the annotation matrix.

## Testing (dependency injection; no monkeypatch; no skipif; Task-level where possible)

- **Submodule (B):** rigged 3-chromosome system with a known `beta` and weights;
  assert `select_alpha_loco` picks the CV optimum, `combine` is additive,
  weighted `fit` recovers `beta` (tol), and `heldout_r2` matches a direct
  recomputation. Pure-numpy, fast, no fixtures.
- **Estimator (C):** synthetic annotations + LD-scores where a known annotation
  subset carries all heritability; simulate chi-square consistent with a chosen
  `tau`; assert recovered `snpvar` ranking recovers the planted ordering
  (Spearman near 1 — relative, not absolute, per decision 1). **No-leakage
  assertions (decision 8):** each chromosome is scored from the opposite parity;
  and the even-scoring `tau`/lambda are invariant to perturbing the chi-square of
  an even chromosome (only odd data may influence them), and symmetrically.
  Task-level, no network, no R (cf. the CT-LDSC synthetic recipe).
- **Coverage:** the snpvar output covers 100% of the annotation-matrix variants,
  so the SUSIE `PriorInfo` join cannot drop variants.
- **Migration gate (see Phasing):** numeric-equivalence test asserting the
  refactored `RidgeAnnotationWeightsTask` reproduces its pre-refactor
  coefficients on a rigged system.
- Follow repo conventions: share column-name constants; no assertions on
  log/error text; prefer isinstance; no structural tests for the generator wiring
  (import-time construction suffices).

## Phasing / sequencing

1. **Component A** (LD-score stream-extract) + schema spike.
2. **Component B** (submodule) + its unit tests.
3. **Component C** (estimator) on B; synthetic test green.
4. **Component D** (generator prior-source injection) + coverage test.
5. **Separately**, migrate `RidgeAnnotationWeightsTask` onto B, gated by the
   numeric-equivalence test. Independent of 1-4 so the new feature does not carry
   the refactor's risk.

## Risks / open questions

1. **LD-score member schema — RESOLVED** (spike 2026-09-21): members carry
   `CHR,SNP,BP,A1,A2` + 187 annotation-named LD-score columns; four-part key
   available; no `base` column (handled via MAFbin-sum, decision 4). See Component
   A.
2. **Munge path:** reuse the rpy2-free genomic_sem Python munge vs the gwaslab
   path — pick whichever already emits a clean HapMap3-independent chi-square + N
   keyed to the LD-score members. (The regression is dense, so munge need not
   restrict to HapMap3.)
3. **Weighting simplification** (decision 4): using the base LD-score as the
   over-counting weight rather than a dedicated regression-SNP weight file may
   shift absolute `snpvar` scale; irrelevant to within-locus-normalized priors,
   but noted for the diagnostics.
4. **Ancestry/LD match** is a *fine-mapping* concern (SUSIE's locus LD), not an
   Approach-2 prior concern; unchanged by this spec. For DecodeME the UKB-derived
   inputs are a defensible reference (~94% of the sample is UKB controls).

## Out of scope (future: Approach 3 + custom annotations)

Bin LD-score computation from a genotype reference panel (or precomputed UKB LD),
SNP binning by predicted `snpvar`, LOCO bin re-estimation (the non-parametric
recalibration), non-baseline-LF annotation ingestion (AlphaGenome etc.), and
annotation-term explanation of the Approach-2 prior via its own `tau`. Approach 2
is built so the estimator's annotation set and LD-scores are injected deps, so
Approach 3 extends rather than rewrites it.
