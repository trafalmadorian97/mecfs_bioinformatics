# Polyfun-explainability ridge surrogate — feasibility spike (THROWAWAY)

Date: 2026-08-25. Goal: can we fit ridge `snpvar_bin ~ 187 baseline-LF annotations`
to get per-annotation weights gamma_c, within the 16GB budget, with a usable R^2?

## Verdict: GO. Clean, high R^2, comfortably in budget.

## Data access (matters for the real download task)
- Raw annotations live in the **11GB** `baselineLF_v2.2.UKB.tar.gz`
  (NOT the 30GB `.polyfun.tar.gz`, which holds only `.l2.ldscore.parquet`).
- Inside: `baselineLF2.2.UKB.<chr>.annot.gz`, LDSC text, one per chromosome,
  ~20-80MB gz each (~1-1.5GB total for all 22). Columns: CHR, BP, SNP(rsid),
  CM, then **187 annotation columns** (coding, conserved, MAF-split
  `_lowfreq`/`_common`, flanking windows, chromatin/eQTL MaxCPP, ...).
- annot.gz files are **interleaved** with 600MB `.l2.ldscore.gz` members and in
  arbitrary chromosome order; gzip isn't seekable, so getting all 22 annot.gz
  means streaming the whole 11GB once and keeping only `*.annot.gz` (~1.5GB).
- No A1/A2 in annot.gz -> **join key to snpvar_meta is rsid (SNP)**.

## Join
- snpvar_meta (chr1_7 + chr8_22) = 19.48M rows; deduped on rsid = 19.44M.
- Join rate annot->meta = **1.000** (essentially every annotated SNP has a prior).

## Target
- `y = snpvar_bin` is the polyfun precomputed prior weight. It is **binned**:
  only ~13 distinct values, 26.7% nonzero. Values are tiny (mean ~7.8e-9) since
  per-SNP h^2 spread over ~19M SNPs — pure scale, not a problem.

## Model result (train chr 2,3,7 -> test chr 19; n_train=4.17M, p=187)
- **Cross-chromosome TEST R^2 = 0.869** (train 0.895). The linear annotation
  surrogate reproduces the polyfun prior well out-of-sample despite the binning
  step ChatGPT worried about. gamma_c decomposition is a compelling approximation.
- R^2 is **flat across alpha 0.1 -> 1e5** — 187 features vs 4M rows, no
  overfitting; regularization is nearly irrelevant (even OLS would do). alpha
  choice is not delicate.
- Top |gamma_c| are biologically sane: conservation (phastCons primate/mammal,
  GERP, LindbladToh), non-synonymous/coding, promoter/enhancer, eQTL & H3K27ac
  MaxCPP, TSS; sensible signs (synonymous negative).

## Memory / compute
- Peak RSS = **8.95GB**, dominated by the 19.4M-row meta table in RAM plus one
  chromosome. The Gram is 187x187 = 0.3MB. Genome-wide (22 chr) stays at the
  same peak because we stream one chromosome at a time -> fits 16GB.
- Production path: a **single** streaming pass suffices — accumulate raw
  cross-products (Sx x^T, Sx, Sx y, Sy, n) per chromosome, then center/standardize
  the 187x187 Gram analytically and solve `(G + alpha I) gamma = b`. (The probe
  used two passes for clarity.) Meta RAM can be cut by joining per-chromosome.

## Implication for the build
- We do NOT need the 30GB tarball. Download task streams the 11GB tarball,
  keeps `*.annot.gz`, converts to (per-chr or concatenated) parquet ~19M x 187.
- Ridge task: streaming Gram -> `gamma_c` weights asset (187 rows: annotation, weight).
- The high linear R^2 justifies the `C_c(i,j) = gamma_c (a_ic - a_jc)` local-contrast
  attribution as the explainability mechanism.
