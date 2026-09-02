# Polyfun explainability — Spec 2: explainability layer

Date: 2026-08-27. Status: design approved (brainstorming), ready for
implementation plan. Consumes Spec 1
(2026-08-25-polyfun-annotation-weights-design.md).

## Context

We fine-map GWAS loci with SUSIE, optionally supplying the polyfun precomputed
per-SNP heritability prior (snpvar_bin, "Approach 1"). When that prior sharpens a
diffuse credible set, we want to explain, in annotation terms, why it favored the
variant it did.

Spec 1 built the surrogate weights: a ridge fit of the published prior on the 187
baseline-LF 2.2.UKB annotations produced a durable weights asset (gamma_raw and
gamma_standardized per annotation, plus a family label), and the annotation matrix
parquet (CHR, BP, SNP, CM + 187 annotations). This spec builds the layer that uses
those weights to explain a fine-mapping result.

The explainability method has three levels (from the design discussion in
chatgpt_conversation_polyfun_explainability.txt):

1. Prior lift: L_i = log(m * pi_i), how strongly the prior favored variant i over
   uniform, where m is the number of variants in the locus and pi_i is the
   normalized prior weight.
2. Local annotation contrast: C_c(i) = gamma_c * (a_ic - abar_c), how much
   annotation c made variant i look more functionally important than its
   plausible competitors.
3. Group ablation (counterfactual reruns).

This spec covers levels 1 and 2 plus a figure. Level 3 (ablation) is out of scope
and deferred to a later spec.

## Definition of done

Given a locus, the new machinery produces, for a matched pair of SUSIE runs (one
with the polyfun prior, one uniform):

- a prior-lift table and per-annotation / per-family contrast tables,
- a stacked explainability figure written as both PNG and SVG,

all green under repo `invoke green`, wired through a new inner + outer asset
generator. The outer generator, invoked at one locus, yields 8 SUSIE runs (4 run
configs x 2 priors). Verified correct on a demonstrator locus before rollout.

## Design decisions (approved)

1. Scope: levels 1 (prior lift) + 2 (local contrast) + the figure. Annotation
   ablation (level 3) is deferred.
2. Reference baseline abar_c: PIP-weighted mean of each annotation over ALL
   variants in the uniform-prior run:
   abar_c = sum_j PIP_j^uniform * a_jc / sum_j PIP_j^uniform.
   Chosen over the "credible-set competitors, j != i" variant because it is
   well-defined when the uniform run is diffuse (the interesting case), does not
   depend on the focal variant, and is stable.
3. Contrast coefficient scale: gamma_raw (raw annotation scale), so C_c(i) is a
   contribution to predicted SNPVAR in the prior's own units. (gamma_standardized
   is for global ranking only and is not used here.)
4. Attribution row set: the union of both runs' 95% credible-set variants.
5. Recombination-rate track: derived as the numerical slope dCM/dBP of the CM
   (centimorgan genetic-map position) column of the Spec 1 annotation parquet over
   the locus variants. No new data dependency. (The Recomb_Rate_10kb annotation is
   MAF-split and standardized, so it is NOT used for the rate line.)
6. Anchor variants: the Manhattan LD-coloring lead is the minimum-p (most
   significant) variant; the focal variant for selecting important families is the
   maximum-PIP-under-polyfun variant. These are deliberately allowed to differ.
7. New code only: a NEW plot task and NEW asset generators, inspired by the
   existing SusieStackPlotTask and fine_mapping_asset_generator. The existing plot
   task and generator are left untouched. The one modification to existing code is
   an additive output on SusieRFinemapTask (decision 8).
8. Per-variant prior provenance: SusieRFinemapTask already computes the aligned,
   filtered per-variant prior array actually passed to susie_rss, then discards it.
   It will additionally write prior.parquet (one row per retained variant, keyed by
   CHR/POS/A1/A2). This is the faithful single source of truth for pi_i and avoids
   re-deriving the allele-key alignment downstream. It is additive (a new file in
   the same DirectoryAsset) and touches no existing consumer.
9. Fail fast on prior coverage gaps: today SusieRFinemapTask's align_data
   inner-joins the polyfun prior, silently dropping any (gwas intersect ld) variant
   that lacks a prior. Since the LD reference and the polyfun prior come from the
   same authors, in theory every such variant is covered. So instead of accepting a
   silent drop, the prior alignment will assert full coverage and RAISE if any
   variant lacks a prior. This flags whether the mismatch ever actually happens. If
   a real locus trips it, we can later add an opt-in to suppress the error and
   accept a reduced variant set — but we do not start with that assumption, and no
   such flag is built now. Consequence: the polyfun run drops nothing, so with the
   prior-independent filtering the two runs retain an identical variant set.

## Structural facts relied on

- Diagnostic variant filtering (kriging logLR / z thresholds) depends only on
  z-scores and LD, not on the prior, so it removes the same variants in both runs.
  Combined with the decision-9 coverage guard (the polyfun prior must cover every
  (gwas intersect ld) variant, else raise), the polyfun run drops nothing, so the
  two runs retain an identical variant set. Downstream code still aligns the two
  runs by variant identity (CHR, POS, A1, A2) via join for ordering safety, but may
  rely on identical membership. Whether real loci actually satisfy full coverage is
  exactly what the decision-9 guard surfaces.
- The annotation parquet carries no alleles, so annotation values are joined to run
  variants by (CHR, BP) — the same allele-free join used by the Spec 1 ridge task.

## Components

### A. SusieRFinemapTask — prior coverage guard + additive prior.parquet

- File: mecfs_bio/build_system/task/r_tasks/susie_r_finemap_task.py.
- Prior coverage guard (decision 9): in the prior alignment (align_data), replace
  the silent inner-join drop with a fail-fast check — every (gwas intersect ld)
  variant must have a matching prior row; raise a clear error (count of missing plus
  a few example variant keys) otherwise. This is the one behavioral change to
  existing code, and it is what surfaces whether the polyfun run would ever drop
  variants. No opt-out flag is added now.
- Additive prior.parquet: after variant filtering, write the per-variant prior
  array alongside the existing outputs, in the same row order as
  filtered_gwas.parquet, keyed by
  GWASLAB_CHROM/POS/EFFECT_ALLELE/NON_EFFECT_ALLELE, value column a module-level
  PRIOR_FILENAME / prior-column constant.
- For the uniform run (prior_info is None) the values are constant (ones); the file
  is still written for a uniform interface.
- Aside from the guard, existing files, filenames, and consumers are unchanged.

### B. PolyfunExplainContrastTask (new)

- Package: mecfs_bio/build_system/task/polyfun_explain/ (new subpackage).
- Deps: [susie_uniform_task, susie_polyfun_task, ridge_weights_task,
  annotation_parquet_task].
- execute:
  1. Read both runs' pip.parquet, filtered_gwas.parquet, combined_cs.parquet, and
     the polyfun run's prior.parquet. Align uniform vs polyfun by variant identity
     (CHR, POS, A1, A2). With the decision-9 guard the two memberships are
     identical; the join guards ordering.
  2. Prior lift: normalize the prior over the locus so it sums to 1 --
     pi_i = w_i / sum_j w_j over the polyfun run's locus variants (so sum_i pi_i =
     1); m = number of locus variants; prior fold = m * pi_i; log fold =
     log(m * pi_i). The figure's prior-fold panel and the display table's lift
     column both use this same m * pi_i.
  3. Load gamma_raw and family from the ridge weights asset. Load annotation values
     from the annotation parquet, predicate-pushdown-filtered to the locus BP
     window, joined to the run variants on (CHR, BP).
  4. abar_c = PIP_uniform-weighted mean of annotation c over all uniform-run
     variants (decision 2).
  5. Per-annotation contrast C_c(i) = gamma_raw_c * (a_ic - abar_c) for each variant
     i in the union of the two runs' 95% credible sets. Per-family contrast =
     sum over annotations in the family.
  6. Focal variant = argmax PIP_polyfun. Top-x families (x default 3) = largest
     signed per-family contrast at the focal variant.
- Outputs (DirectoryAsset): prior-lift table (parquet), per-annotation contrast
  table (long parquet), per-family contrast table (parquet), the display table
  (parquet, see below), and a machine-readable record of the focal variant +
  selected families (so the plot task and tables agree on selection). Filename
  constants are module-level.
- Display table (parquet, docs-facing). This is the human-readable summary rendered
  in the docs via the data_table macro in main.py (Tabulator, client-side), so it
  is parquet, NOT markdown.
  - Rows: all variants in the union of the two runs' 95% credible sets (the
    attribution row set).
  - Short, abbreviated column names so columns fit on screen in the docs. Columns,
    in order (name -> meaning):
    - chr    -> chromosome (Int32)
    - pos    -> position (Int32)
    - ea     -> effect allele
    - nea    -> non-effect allele
    - cs_pf  -> polyfun-run credible-set number containing the variant (null if in
                none)
    - cs_u   -> uniform-run credible-set number containing the variant (null if in
                none)
    - pip_pf -> PIP, polyfun run
    - pip_u  -> PIP, uniform run
    - lift   -> prior fold m * pi_i (same quantity as the figure's prior-fold panel)
    - then x columns (x default 3), one per important family, each giving
      sum_{c in family} gamma_raw_c * a_ic for that variant (the raw scaled family
      value the family panels plot -- NOT the contrast). Named by a short family
      abbreviation. The x families are the same ones selected from the focal
      (max-PIP-polyfun) variant that drive the figure panels, so table and figure
      agree.
  - The credible-set number is the L-index (1-based) from that run's
    combined_cs.parquet cs column; if a variant somehow falls in more than one, use
    the lowest-numbered.
  - Sorted: descending by pip_pf.
  - chr and pos written as Int32 (the data_table macro requires Int32; Parquet
    INT64 decodes to JS BigInt, which Tabulator cannot format).
- Meta derived from the polyfun SUSIE run's ResultDirectoryMeta (reuse
  trait/project), raising on unknown meta.

### C. PolyfunExplainPlotTask (new, standalone)

- Package: mecfs_bio/build_system/task/polyfun_explain/. New file, inspired by
  susie_stacked_plot_task.py but independent of it; the old task is unchanged.
- Deps: [susie_uniform_task, susie_polyfun_task, contrast_task,
  annotation_parquet_task, gene_info_task]. (The contrast task supplies the focal
  variant and selected families so the figure and tables are consistent.)
- Shared-x (genomic position) stacked panels, top to bottom:
  1. Manhattan -log10 p, points colored by LD r^2 with the minimum-p lead variant
     (LD read from the run's filtered_ld.npy), with local recombination rate
     (dCM/dBP from the annotation CM column) overlaid on a secondary y-axis.
  2 .. (1+x). x family panels (default 3), one per selected important family:
     plot sum_{c in family} gamma_raw_c * a_ic across the locus (the raw scaled
     annotation value, NOT the contrast; the contrast is read off the profile,
     i.e. how the focal variant stands out from its neighbors).
  3. Prior fold m * pi_i (log scale).
  4. PIP, uniform-prior run.
  5. PIP, polyfun-prior run.
  6. Gene track (reuse the gene-track drawing approach from the old plot).
  With x = 3 this is 8 panels.
- The figure's variant spine is the polyfun run (where the focal variant lives);
  the uniform PIP panel is a left-join of uniform PIP onto that spine.
- Writes BOTH plot.png and plot.svg into its output directory (DirectoryAsset); the
  SVG is zoomable. Filename constants are module-level.
- x (number of family panels) is a task parameter defaulting to 3.

### D. Inner asset generator (new)

- File: mecfs_bio/asset_generator/polyfun_explain_fine_mapping_asset_generator.py,
  inspired by fine_mapping_asset_generator.py.
- Signature: (locus + shared harmonized-sumstats / LD inputs, one run config).
  A run config is the existing per-run knobs (max_credible_sets,
  z_score_filtering_threshold, and a label).
- Produces a matched pair of SusieRFinemapTask: one with prior_info = polyfun
  (PriorInfo over COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS), one with
  prior_info = None (uniform), identical otherwise. Then a PolyfunExplainContrastTask
  and a PolyfunExplainPlotTask over that pair.
- Returns a frozen group (the pair + contrast + plot) for the config.

### E. Outer asset generator (new)

- Same file. Reuses the existing shared setup (optimal UKBB LD interval, renamed LD
  labels, harmonized sumstats) once per locus, then calls the inner generator for
  each of the 4 existing run configs: L=1, L=2, L=10, L=10-strict
  (z_score_filtering_threshold = 1.0). Result: 8 SUSIE runs + 4 contrast/plot sets
  per locus.
- Returns a frozen task group with terminal_tasks() enumerating all terminal
  artifacts (the 8 runs, the 4 contrast tables, the 4 plots), mirroring
  BroadFineMapTaskGroup.

### F. Rollout (after the machinery is verified — separate step, not this build)

- Copy mecfs_bio/assets/gwas/me_cfs/decode_me/analysis/fine_mapping/with_palindromes
  to a new sibling folder and rewire it to call the new outer generator at every
  major DecodeME locus. Done only after the generator is verified correct on a
  demonstrator locus.

## Memory / compute notes

- Contrast task: reads two small run directories, a 187-row weights table, and a
  locus-windowed slice of the annotation parquet (predicate pushdown on CHR/BP, a
  few thousand rows x 187). Trivial memory.
- Plot task: matplotlib over a few thousand locus variants; PNG + SVG. Trivial.
- Both single-machine, default pixi env, no Docker. (SusieRFinemapTask itself runs
  susieR via rpy2 as today.)

## Testing (dependency injection; no monkeypatch; no skipif; Task-level)

- SusieRFinemapTask prior.parquet: extend the task's existing test to assert the
  new file is written with one row per retained variant and the expected key
  columns; for a polyfun-prior fixture, values match the aligned prior; for a
  uniform run, constant.
- SusieRFinemapTask prior coverage guard: a fixture whose prior table omits one
  (gwas intersect ld) variant makes execute raise. Assert that it raises, not the
  message text (repo convention).
- PolyfunExplainContrastTask: synthetic two-run fixtures (pip / filtered_gwas /
  combined_cs / prior parquet) + a synthetic ridge weights table + synthetic
  annotation parquet, with a known linear prior so abar_c and C_c(i) are checkable
  in closed form; assert the focal variant and top-family selection. Also assert the
  display table: rows = the union of the two credible sets, sorted descending by
  pip_pf, abbreviated column names, cs_pf/cs_u carrying each run's credible-set
  number (null when the variant is in the other run's CS but not this one's), the x
  family columns named by (and matching) the selected families, and chr/pos dtype
  Int32. Task-level, injected fetch.
- PolyfunExplainPlotTask: smoke test that plot.png and plot.svg both land in the
  output directory given synthetic inputs; no pixel assertions.
- Generators: construct the outer group and assert it wires 8 SUSIE runs and 4
  contrast/plot sets with distinct asset ids; no execution.
- Follow repo conventions: Task is the public API; no assertions on log/error text;
  share column-name / filename constants between producer and test; prefer
  isinstance over re-deriving.

## Out of scope

- Level 3 annotation ablation (counterfactual family reruns and delta-PIP /
  delta-CS-size / delta-entropy endpoints).
- Any change to the existing SusieStackPlotTask or fine_mapping_asset_generator.
  The only edits to existing code are on SusieRFinemapTask: the decision-9 prior
  coverage guard and the additive prior.parquet output.
- The per-locus rollout folder (Component F), performed after verification.
- An opt-in to suppress the decision-9 prior-coverage error and accept a reduced
  variant set. Added only if a real locus trips the guard.
