# DecodeME outputs before and after the switch to the gnomAD panel

PR #1233 switched DecodeME's genome-reference harmonization from the 1000 Genomes EUR panel to
the gnomAD v2.1.1 `nfe_nwe` panel. It also added a palindrome AF-distance cap of 0.1. This page
compares key DecodeME outputs before and after the switch.

## Method

- `record_decode_me_baseline.py` (2026-10-07, before the merge) recorded the outputs that the switch
  would rebuild. The scope is the cis-pQTL MR, LDSC, S-LDSC and HBA/ensembl MAGMA analyses, plus the
  `with_palindromes`, `polyfun_explainability` and `polyfun_approach_2_l2_sldsc_prior` fine-mapping
  packages. Within that scope it took each DecodeME harmonization and everything downstream of it.
  It recorded 510 built assets, with per-file md5s. Outputs of up to 500 MiB were copied into
  `assets/baselines/decode_me_gnomad_panel_switch/`, which is gitignored.
- `rebuild_decode_me.py` (2026-10-09, on main after the merge) rebuilt exactly those 510 assets.
- `compare_decode_me_to_baseline.py` compares md5s, then compares the contents of the key output
  types, matching rows on a natural key.
- `harmonized_row_changes.py` classifies the rows that the switch added, dropped or re-oriented in
  the harmonized tables.

Each rebuilt task also picked up any other code changes merged to main, such as polars 2.0. A small
numeric difference therefore does not necessarily come from the panel.

## Summary

52 of the 510 assets are byte-identical and 458 changed. No asset is missing.

| Analysis | Result |
|---|---|
| S-LDSC (6 cell-type analyses, multiple-testing tables, scatter plots) | Byte-identical |
| HBA MAGMA (gene analysis, gene-set analysis, conditional analysis) | Numerically identical; only logs and intermediate SNP files changed |
| Ensembl specific-tissue MAGMA | Not rebuilt: it assigns rsIDs through dbSNP151 and does not use the harmonized data |
| LDSC h² | h2_liab 0.081358 → 0.081425 (+0.08%); intercept, λGC and mean χ² identical; cause not identified |
| cis-pQTL MR | All 1,506 exposure/method pairs present in both runs have identical b, SE and p; 95 nominal hits before and after; 8 exposures gained an instrument and TNN lost one; the 2 corrected-significant hits are unchanged |
| L2 S-LDSC prior | Regularization alphas identical; per-SNP snpvar correlation with the old prior is 0.9989 |
| SuSiE / PolyFun (107 runs) | Credible-set structure unchanged in 104 runs; most PIP changes are at most 0.02 |

## Harmonized tables

`decode_me_gwas_1_genome_reference_harmonized` grew from 8,371,657 to 8,628,839 rows.

| Change | Indels | Palindromic SNVs |
|---|---|---|
| Added | 248,601 | 29,981 |
| Dropped | 12,745 | 8,655 |
| BETA sign changed | 14 (effect allele changed) | 301 (alleles unchanged) |

- **Added rows:** these are mostly ambiguous indels and palindromes that the gnomAD panel resolves.
  The 1000 Genomes panel often had no record for them.
- **Dropped rows:** these include palindromes now cut by the AF-distance cap.
- **Sign changes without an allele change:** these are palindromes whose strand decision reversed.
  After REF orientation they keep the same alleles, but get −BETA and 1−EAF.

`decode_me_gwas_1_keep_ambiguous_genome_reference_harmonized` grew from 8,555,668 to 8,791,524 rows.
It has the same indel changes, and its 676 palindrome strand reversals change sign; it keeps every
palindrome either way.

## Fine mapping: the chr20 locus

At chr20:47.0–48.2 Mb, three uniform-prior L=10 runs lost a credible set:

- `decode_mechr20_47000000_48200000_palindromes_keep_susie_finemap` (the main fine-mapping run)
- `decode_me_polyfun_explainchr20_47000000_48200000_palindromes_keep_l10_susie_uniform`
- `decode_me_polyfun_explain_l2_sldsc_priorchr20_47000000_48200000_palindromes_keep_l10_susie_uniform`

| Credible set | Before | After |
|---|---|---|
| L2 | 1 variant at 47,743,125, PIP 1.0 | Unchanged |
| L3 | 95 variants, max PIP 0.057, minimum purity 0.94 | 86 variants (84 shared), max PIP 0.107, minimum purity 0.97 |
| L1 | 119 variants, max PIP 0.067, minimum purity 0.65, spanning 47.51–47.91 Mb | Gone |

In the old fit, L1 and L3 covered almost the same region and split one diffuse signal between
them. The new fit assigns that signal to a single, tighter set. These are the only runs with a
PIP change above 0.05; the largest is 0.066. Documentation built from the main run will show 2
credible sets at this locus instead of 3.

## Reproducing

```
pixi r python -m experiments.claude.decode_me_panel_switch.rebuild_decode_me
pixi r python -m experiments.claude.decode_me_panel_switch.compare_decode_me_to_baseline
pixi r python -m experiments.claude.decode_me_panel_switch.harmonized_row_changes
```

All three need the baseline under `assets/baselines/decode_me_gnomad_panel_switch/`, which
`record_decode_me_baseline.py` creates.
