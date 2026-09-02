# Polyfun explainability — Spec 1: annotation weights layer

Date: 2026-08-25. Status: design approved (brainstorming), ready for implementation plan.

## Context

We fine-map GWAS loci with SUSIE, optionally supplying a polyfun precomputed
per-SNP heritability prior (`snpvar_bin`, "Approach 1"). We want explainability:
when the prior concentrates a diffuse credible set, on what annotation basis did
it do so?

The chosen mechanism (see the ChatGPT discussion and the brainstorm) is a linear
**surrogate model** of the published prior:

    SNPVAR_i  ~  alpha + sum_c gamma_c * a_ic

fit by ridge over the 187 baseline-LF 2.2.UKB annotations, genome-wide. The
per-annotation weights `gamma_c` then drive a local-contrast attribution in
Spec 2: `C_c(i) = gamma_c * (a_ic - baseline_c)`.

This document specifies **only Spec 1**: acquiring the annotations and producing
a durable `gamma_c` weights asset. Spec 2 (SUSIE prior output, two-run
generators, contrast/summary tables, 8-panel figure) is designed separately and
consumes this spec's weights asset.

### Feasibility (validated by throwaway spike, experiments/claude/polyfun_explain_probe)

- Raw annotations live in the **11GB** `baselineLF_v2.2.UKB.tar.gz`
  (`s3://broad-alkesgroup-ukbb-ld/UKBB_LD/...`), NOT the 30GB `.polyfun.tar.gz`
  (which holds only `.l2.ldscore.parquet`). Members
  `baselineLF2.2.UKB.<chr>.annot.gz` are LDSC text: `CHR, BP, SNP, CM` + 187
  annotation columns, one file per chromosome, BP-sorted.
- Join key to `snpvar_meta` is **rsid (SNP)** (annot.gz carries no A1/A2).
- The annotation set **covers all ~19M polyfun prior variants** (a slight
  superset; 100% of meta variants have annotations; a few thousand extra
  annotation variants and a handful of duplicate rsids per chromosome).
- `y = snpvar_bin` is a **binned** step function (~13 distinct values, ~27%
  nonzero, tiny scale ~1e-8). We are explaining the prior actually used, not the
  original S-LDSC tau_c; docs must say so.
- Cross-chromosome test R^2 ~= **0.87**; R^2 is essentially flat in alpha
  (0.1 .. 1e5); top |gamma_c| are biologically sane (conservation, coding,
  promoter/enhancer, eQTL/H3K27ac MaxCPP, TSS).
- Streaming primal-Gram fit keeps peak RSS ~9GB (dominated by the meta table,
  not the regression); genome-wide stays the same -> fits the 16GB budget.

## Definition of done

`RidgeAnnotationWeightsTask` is green (repo `invoke green`), and on the real data
reproduces the spike: held-out R^2 ~0.87 and biologically sane top annotations.
The weights asset is a small parquet table one row per annotation.

## Design decisions (approved)

1. **Attribution grain (affects Spec 2, recorded here):** hybrid — ridge always
   on 187 individual annotations; Spec 2 aggregates contributions to ~10
   families for headline/figure and keeps all 187 in a detail table. This spec
   curates the 187->family map and carries `family` in the weights table.
2. **Storage:** ONE parquet sorted by (CHR, BP), columns
   `CHR, BP, SNP, CM, <187 annotations:f32>`. Enables both a cheap streaming
   scan for the ridge and predicate-pushdown per-locus lookups in Spec 2.
3. **Tarball at rest:** download and store the 11GB tarball as a reference asset
   (reuse the tested downloader + md5 + build cache); mark it a `path_remap`
   candidate for `/mnt/d`. Not a bespoke stream-and-discard.
4. **alpha selection:** leave-one-chromosome-out over a small grid, refit on all;
   report held-out R^2. (R^2 is flat in alpha, so this is for defensibility.)
5. **Coefficient scale:** store `gamma_raw` (raw annotation scale, for Spec 2's
   contrast) AND `gamma_standardized` (for global importance ranking). `intercept`
   goes in the diagnostics sidecar. Per-column `mean_c`/`std_c` are NOT stored:
   they are redundant (`gamma_standardized = gamma_raw * std_c`) and nothing
   downstream needs them (Spec 2's contrast uses raw annotation values with
   `gamma_raw`).

## Components

### A. Tarball download — reuse `DownloadFileTask`

- `url = https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.tar.gz`
- `md5_hash` pinned (compute during implementation; S3 ETag is multipart so not a
  plain md5 — compute md5 after first download and pin it).
- `meta`: `ReferenceFileMeta(group="polyfun", sub_group="annotations",
  sub_folder="raw", id="baseline_lf_2.2_ukb_annotations_tarball", extension=".tar.gz")`.
- Location: `mecfs_bio/assets/reference_data/polyfun/annotations/`.

### B. `BuildBaselineLFAnnotationParquetTask` (new)

- Package: `mecfs_bio/build_system/task/annotation_weights/` (new subpackage).
- Deps: `[tarball_task]`.
- `execute`:
  1. Extract only `*.annot.gz` members from the tarball via `execute_command`
     (`tar`; determine the exact member paths by listing first, then extract by
     explicit name — the spike showed `--wildcards '*.annot.gz'` unreliable,
     exact names work).
  2. Build a polars LazyFrame = vertical concat, in chromosome order 1..22, of
     `pl.scan_csv(annot_gz, separator="\t", infer_schema_length=None)` per file,
     with the 187 annotation columns cast to Float32 and a dedup on `SNP` (drop
     the handful of multiallelic rsid dups). Chromosome order + per-file BP sort
     make the result globally (CHR, BP)-sorted.
  3. Write with `write_df_according_to_format(lazy, out_path, ParquetOutFormat())`
     — the default (no write options) path uses `df.sink_parquet`, a STREAMING
     write, so the ~19M x 191 frame is never fully materialized. (The
     `ParquetWriteOptions` path is non-streaming — `.collect().to_arrow()` — so
     we deliberately use the default and forgo byte-stream-split here.) If
     streaming a scan over gzip proves problematic, fall back to eager
     per-chromosome reads still sunk through the same writer.
- Output: single `FileAsset` parquet, `ReferenceFileMeta`
  `id="baseline_lf_2.2_ukb_annotations"`,
  `read_spec=DataFrameReadSpec(DataFrameParquetFormat())`. No json sidecar — the
  187 annotation names are read from the parquet schema where needed
  (family-map test, ridge task). Note the ~19M x 191 parquet is a `path_remap`
  candidate.
- `create()` DERIVES this meta from the tarball task's `ReferenceFileMeta`
  (reuse group/sub_group/sub_folder; raise on unknown meta) — the
  `CompressedCSVToParquetTask.create` pattern. Only id/extension/read_spec differ.

### C. Annotation -> family map (new curated asset)

- `mecfs_bio/assets/reference_data/polyfun/annotations/annotation_families.py`.
- A rule-based classifier `family_for_annotation(name) -> AnnotationFamily`
  assigning all 187 annotations to one of **11 published-grounded families**:
  `non_synonymous`, `coding`, `conserved`, `promoter_or_enhancer`,
  `histone_marks`, `repressed`, `open_chromatin`, `maf_bins`,
  `ld_related_continuous`, `molecular_qtl`, `other`.
- Taxonomy is grounded in the literature, not invented (checked: no paper ships
  a reusable per-annotation family table, and the polyfun GitHub repo has no
  grouping — only LDSC's per-annotation "category" sense):
  - The 7 functional-group names (non_synonymous, coding, conserved,
    promoter_or_enhancer, histone_marks, repressed, other) are the grouping the
    polyfun authors themselves use in the sub-additive simulation of their
    Supplementary Note (Weissbrod et al. 2020).
  - maf_bins and ld_related_continuous are the MAF-bin and LD-related continuous
    groups of Gazal et al. 2017 (the Continuous rows of Gazal 2018 Table S1).
  - molecular_qtl are the MaxCPP molecular-QTL annotations of Hormozdiari 2018.
  - open_chromatin (DHS_Trynka, DHS_peaks_Trynka, FetalDHS_Trynka, DGF_ENCODE) is
    the ONE deliberate refinement of polyfun's scheme (which lumps these in
    "others"), broken out for the explainability figure's accessibility panel;
    TFBS/CTCF/Transcribed/Intron stay in `other`. Documented as a deviation.
  - Per-annotation assignment follows annotation names + their source datasets
    (Gazal 2018 Table S1). MAF-split (`_lowfreq`/`_common`) and `.flanking.500`
    suffixes do not change the family.
- A `Literal` type alias `AnnotationFamily` for the family names; module docstring
  carries the four citations.
- **Test:** keys == the 187 annotation names (committed name-list constant,
  checkable against the parquet schema); every annotation resolves to a valid
  family; open_chromatin membership == exactly the DHS/DGF accessibility set.

### D. `RidgeAnnotationWeightsTask` (new) + weights asset

- Package: `mecfs_bio/build_system/task/annotation_weights/`.
- Deps: `[annotation_parquet_task, COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS]`.
- `execute` (single streaming pass, no full design matrix in memory):
  1. Load `snpvar_meta` (SNP, snpvar_bin), dedup rsid, into a polars frame (read
     natively, not via pandas, to keep RAM well under 16GB).
  2. Scan the annotation parquet **one chromosome at a time** (CHR filter with
     predicate pushdown; a chromosome can itself be read in sub-batches to cap
     memory). Accumulate cross-products into **per-chromosome partials** — for
     each chromosome store `n, Sx (p), Sxx_gram (p x p) = sum a a^T,
     Sxy (p) = sum a*y, Sy, Syy` over the 187 annotations `a`, after joining
     `snpvar_bin` on rsid. Per-chromosome partials are what make LOCO free (any
     train fold = sum of the other chromosomes' partials); the genome-wide fit is
     the sum of all 22. All 22 partials together are ~6MB.
  3. Center + standardize the Gram analytically:
     `mean = Sx/n`, `var = Sxx_diag/n - mean^2`, `sd = sqrt(var)`;
     `G_std = (Sxx_gram - n*mean mean^T) / (sd sd^T)`,
     `b_std = (Sxy - mean*Sy) / sd`.
  4. **alpha by LOCO:** for each chromosome held out, combine the other
     chromosomes' partials into `G_train, b_train`, solve for each alpha in a
     small grid, score R^2 on the held-out chromosome's partials (R^2 computable
     from sufficient stats: `SS_res = yty - 2 g.b + g.G.g`). Pick the alpha with
     best mean held-out R^2.
  5. Refit on all chromosomes at the chosen alpha -> `gamma_std`.
     `gamma_raw = gamma_std / sd`; `intercept = mean_y - sum(gamma_raw * mean)`.
  6. Write outputs into the scratch directory (returned as a `DirectoryAsset`):
     - `weights.parquet` via `write_df_according_to_format(..., ParquetOutFormat())`:
       rows = 187 annotations, columns `annotation, gamma_raw,
       gamma_standardized, family`.
     - `diagnostics.json`: `alpha, intercept, heldout_r2_per_chrom,
       mean_heldout_r2, n_variants`.
- Output: `DirectoryAsset`, `ReferenceDataDirectoryMeta`
  `id="baseline_lf_2.2_ukb_annotation_ridge_weights"`. (A directory, not a file,
  because the task emits both the weights parquet and the diagnostics json.)
  `create()` DERIVES group/sub_group/sub_folder from the annotation-parquet
  dependency's `ReferenceFileMeta` (raise on unknown meta), same pattern as
  Component B. Filename constants for the two members are module-level so Spec 2
  and tests reference them, not string literals.
- Weights-asset wiring instance lives in
  `mecfs_bio/assets/reference_data/polyfun/annotations/`.

## Memory / compute notes

- Annotation parquet build: streamed row-group writes; peak ~ one chromosome
  (~1.6M x 191 f32 ~= 1.2GB) + writer buffers.
- Ridge: peak ~ meta frame (~1-2GB native polars) + one batch + the p x p Gram
  (0.3MB). Comfortably < 16GB. The per-chromosome partial Grams for LOCO are 22 x
  (187x187 f64) ~= 6MB total — negligible.
- Both tasks are single-machine, no Docker, default pixi env.

## Testing (all via dependency injection; no monkeypatch; no skipif)

- **BuildBaselineLFAnnotationParquetTask:** synthetic `.tar.gz` fixture built in
  the test — 2 fake `baselineLF2.2.UKB.<chr>.annot.gz` (a few SNPs, a few
  annotation columns) + one decoy `.l2.ldscore.gz`. Inject a fake `WF`
  downloader that yields the local fixture. Assert: decoy skipped; output columns
  = keys + annotations; rows deduped and (CHR, BP)-sorted.
- **RidgeAnnotationWeightsTask:** synthetic annotation parquet + a snpvar built as
  a known linear combination of the annotations (+ small noise). Assert recovered
  `gamma_raw ~= truth` (tolerance) and mean held-out R^2 ~ 1. Task-level.
- **annotation_families:** assert keys == the 187 annotation names (from the
  committed name-list constant), and every value is a valid `AnnotationFamily`.
- Follow repo conventions: Task-level tests (Task is the public API); no
  assertions on log/error text; share column-name constants between build and
  assertion; prefer isinstance over re-deriving.

## Out of scope (Spec 2)

SUSIE emitting per-variant prior + full PIP; inner generator running SUSIE twice
(with/without prior); PIP-weighted-all-variants baseline
`baseline_c = sum_j PIP^noannot_j a_jc / sum_j PIP^noannot_j`; per-annotation and
family-aggregated contrast tables; display summary table; 8-panel stacked figure;
outer generator (8 SUSIE runs). Row set for attribution = union of both runs'
95% credible-set SNPs. These are designed in Spec 2.
