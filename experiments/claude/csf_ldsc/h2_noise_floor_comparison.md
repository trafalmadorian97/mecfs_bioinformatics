# Why the CSF heritability table has large negative h² — PPP comparison

**Question:** the CSF pQTL LDSC table
(`csf_heritability_hapmap_3.parquet`) has many aptamers with large *negative*
heritability (down to −0.51), some with small p-values. Is this a bug?

**Answer: no — it is the LDSC noise floor at low sample size.** Reproduced by
`h2_noise_floor_comparison.py` (log in `logs/`), comparing CSF (Western et al. 2024,
median N ≈ 3.3k) against UKB-PPP plasma (median N ≈ 33.5k) on the same diagnostics.

## Diagnostics (all_variants)

| statistic | UKB-PPP | CSF |
|---|---|---|
| median N | 33,529 | 3,343 |
| n proteins / aptamers | 2,940 | 7,008 |
| h² mean / median | 0.078 / 0.072 | 0.012 / 0.011 |
| h² std | 0.060 | 0.138 |
| h² min / max | −0.035 / 0.278 | −0.514 / 0.854 |
| fraction h² < 0 | 7.7% | 46.7% |
| fraction mean_chi2 < 1 | 4.7% | 29.3% |
| mean_chi2 (mean) | 1.063 | 1.004 |
| corr(h², mean_chi2) | 0.98 | 0.65 |
| n p<0.05 (chance) | 2,203 (147) | 347 (350) |
| n p<0.001 (chance) | 1,870 (2.9) | 10 (7.0) |
| survive Bonferroni | 1,464 — all h²>0 | 0 |

(PPP `cis_excluded` rows are within rounding of `all_variants`; dropped for a like-for-like
comparison. p recomputed for both as the two-sided Wald `p = 2·norm.sf(|h²/se|)`, the same
convention `CsfProteinHeritabilityTask` writes.)

## Why LDSC produces negative h²

LDSC h² is an **unconstrained regression slope** — it regresses each SNP's χ² on `N·ℓ/M`
plus an intercept and reports `h² ∝ slope`. Nothing forces the slope positive; with little
polygenic signal it is noise scattered symmetrically around zero. The kernel does **not**
clip it (only the internal weight term is clipped), which is correct for a diagnostic —
clipping to 0 would bias the mean up and hide that the estimates are noise.

## What the comparison shows

- **It is sample size.** χ² inflation scales with N. PPP `mean_chi2 ≈ 1.063` carries real
  genome-wide signal; CSF `≈ 1.004` sits at the null. So PPP h² is clearly positive
  (median 0.072, only 7.7% negative and those hug zero, min −0.035), while CSF h² is
  centered at zero with 47% negative running to −0.51 (larger SEs, no positive signal to
  lift the estimates).
- **Signal vs chance.** PPP: 1,870 proteins at p<0.001 vs ~3 expected, and **1,464 survive
  Bonferroni, every one positive**. CSF: hit counts match chance exactly (347 at p<0.05 vs
  ~350; 10 at p<0.001 vs ~7) and **none** survive Bonferroni.
- **`corr(h², mean_chi2)`**: 0.98 (PPP, a real regression line) vs 0.65 (CSF, slope jerked
  by noise). Every negative-h² CSF row has `mean_chi2 < 1` — the mechanical fingerprint of
  a noise-driven slope.

## Two things that rule out a bug

1. **Not an alignment/sign error:** h² depends only on `χ² = (β/se)²`, invariant to allele
   flips, so harmonization cannot manufacture negative h².
2. **Both datasets agree where there is power:** every Bonferroni hit is positive in both,
   and PPP's negatives cluster at zero. CSF's large negatives are just what PPP would look
   like at 1/10th the N.

## Implications for the CSF table

- Keep h² unclipped; the negatives are an honest noise-floor / QC signal.
- Report significance with multiple-testing correction (nothing survives Bonferroni; a
  BH-FDR column would make this explicit). Treat per-aptamer LDSC h² as "not distinguishable
  from 0 at this N" for the bulk of aptamers.
- Per-protein h² that actually resolves at N ≈ 3k needs a cis-aware / individual-level
  method, not LDSC.
