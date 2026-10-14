# Palindrome AF distance predicts disagreement between panels

Evidence for `palindrome_max_af_distance` in `GenomeReferenceHarmonizationOptions`.

Source: `experiments/claude/gnomad_af_reference/palindrome_eaf_distribution.py`, run 2026-10-06 on
DecodeME (`DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED` pre-harmonization table, untrusted)
with four panel choices: 1000 Genomes `eur`, gnomAD v2.1.1 `nfe` and `nfe_nwe`, Pan-UKBB `ukb_eur`
(after the chr21/chr22/X swap-profile filter). The palindrome rule at the time was gwaslab's
side-of-0.5 rule with both MAF thresholds at 0.4 and no distance check.

## Definitions

- **Distance**: |EAF − panel AF| after harmonization, with EAF in the output orientation and
  panel AF for REF = NEA, ALT = EA. A keep decision gives |EAF − AF|; a strand flip gives
  |(1 − EAF) − AF|.
- **Opposite strand**: among palindromic SNVs that both panel choices resolve, the fraction for
  which one keeps the strand and the other flips it. At least one of the two is then wrong.

## Disagreement by distance (gnomAD `nfe` distance, compared with 1000 Genomes `eur`)

| distance from panel AF | palindromes both resolve | opposite strand | fraction |
|---|---|---|---|
| ≤ 0.02 | 977,176 | 197 | 0.02% |
| 0.02–0.05 | 72,447 | 63 | 0.09% |
| 0.05–0.1 | 1,851 | 40 | 2.2% |
| 0.1–0.2 | 922 | 88 | 9.5% |
| 0.2–0.3 | 345 | 47 | 14% |
| > 0.3 | 154 | 40 | 26% |

The other pairs show the same pattern:

| distance bin | 1000g vs gnomad_nfe | pan_ukbb vs gnomad_nfe | gnomad_nfe vs pan_ukbb |
|---|---|---|---|
| ≤ 0.02 | 0.016% | 0.015% | 0.014% |
| 0.02–0.05 | 0.05% | 2.3% | 0.06% |
| 0.05–0.1 | 0.8% | 4.2% | 1.8% |
| 0.1–0.2 | 10% | 13% | 9.9% |
| 0.2–0.3 | 21% | 21% | 14% |
| > 0.3 | 20% | 33% | 31% |

(The first-named choice supplies the distance.) `gnomad_nfe` and `gnomad_nfe_nwe` never disagree,
as they share records.

## Cost of a distance cap

Resolved palindromes dropped, out of about 1.06–1.09 million per choice:

| cap | 1000g_eur | gnomad_nfe | gnomad_nfe_nwe | pan_ukbb_eur |
|---|---|---|---|---|
| 0.1 | 1,679 | 2,832 | 2,777 | 895 |
| 0.05 | 9,750 | 5,744 | 5,067 | 1,991 |

A cap of 0.1 removes the bins where roughly one call in ten to one in four is contradicted by
another panel, at a cost of about 0.2% of palindromes.

## Panel AF exactly 0 or 1

The side-of-0.5 rule accepts a panel AF of 0 or 1 (MAF 0). DecodeME has no such palindrome with
MAF below 0.01, so these are never "rare in both" cases; many are common in DecodeME (for
`gnomad_nfe`, 110 of 470 have MAF above 0.2). The distance cap drops those with MAF above 0.1 and
keeps those with MAF at most 0.1.
