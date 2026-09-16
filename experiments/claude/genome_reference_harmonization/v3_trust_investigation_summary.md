# Ambiguous-indel orientation in "trusted" summary statistics

## Question

The genome-reference harmonization design trusts a table (keeps its ambiguous
indels in source orientation) when it is 100% consistent with the reference on
SNVs and indels. V3 tested that assumption on build-38 DecodeME, whose source
orientation we treat as truth. It found 11 ambiguous indels that the
frequency-based rules would swap. We then asked: is that specific to DecodeME,
or would other "fully trusted" summary statistics show the same thing?

## Method

An ambiguous indel is one where both the insertion and the deletion orientation
match the reference genome (gwaslab class indel_both). For the eight other
locally built gwaslab pickles, in their source build (no liftover), we:

1. picked hg19 or hg38 by whichever matched more chromosome-1 SNVs;
2. measured trust genome-wide with the default thresholds;
3. if trusted and EAF present, resolved every ambiguous indel with the
   stringent rules (distance 0.02, margin 0.3) and counted "swap" decisions.

A swap means the panel allele frequencies decisively contradict the source
orientation: the effect-allele frequency matches the complement of the opposite
orientation's panel frequency. On a table we are calling trusted, that is
evidence the source is NOT perfectly reference-aligned for indels.

## Which datasets could be tested

| Dataset | Trust check | Testable? |
|---|---|---|
| DecodeME GWAS-1 (build 38, regenie) | trusted | yes (truth set) |
| Kerrebijn fibromyalgia | trusted (hg38) | yes |
| MVP myocardial infarction | trusted (hg38) | yes |
| Bellenguez Alzheimer's | trusted (hg38) | yes |
| deCODE seropositive RA | 100% SNV but no indels | no indels |
| PGC 2022 schizophrenia | no EAF column | no |
| Yengo height | HapMap3 SNVs, 99.98% | not trusted |
| Johnston pain, Han asthma | ~25% consistent (build/format issue) | no |

All three testable non-DecodeME sets are GWAS Catalog harmonised files
(.h.tsv.gz); Bellenguez's EAF is the catalog's hm_effect_allele_frequency. On
chr1 SNVs their EAF tracks the 1KG EUR panel AF closely (median absolute
difference 0.002-0.009), so EAF is comparable to the panel.

## Result: contradicting swaps by dataset

At distance 0.02, margin 0.3:

| Dataset | Ambiguous indels | Swap ("suspicious") | Proportion |
|---|---|---|---|
| DecodeME build 38 | 769,705 | 11 | 0.0014% |
| MVP myocardial infarction | 944,565 | 249 | 0.026% |
| Bellenguez Alzheimer's | 1,218,978 | 26,829 | 2.2% |
| Kerrebijn fibromyalgia | 1,592,151 | 131,558 | 8.3% |

Requiring a panel record for both readings does not fix it: Kerrebijn still has
102,760 contradicting swaps, MVP 110.

## What the examples show

In most contradicting rows the frequencies favour the swap, and the source
orientation is what looks wrong.

- Kerrebijn 1:869598: source EA=T, NEA=TA, EAF 0.92; panel only has T>TA at
  0.078. T is the reference allele at 92%, but the source labels TA as NEA.
- MVP 1:8622779: source EA=AT, NEA=A, EAF 0.213; panel has AT>A at 0.791. AT is
  the reference at 21%: the "reference allele is minor, frequencies mirror"
  pattern also seen in DecodeME.
- Bellenguez: mostly EAF ~0.999 with panel records near AF 0; the frequency
  evidence still points to a swap but at near-monomorphic sites is weak.

## Conclusions

1. The failure is not DecodeME-specific. Passing the 100% SNV+indel consistency
   check does not guarantee ambiguous indels follow the reference convention.
   The NEA==REF check cannot distinguish orientations when both alleles match
   the genome, so a trusted table's ambiguous-indel labels are not proof of
   correct orientation.

2. There is a wide, clean gap between DecodeME (0.0014%) and every other trusted
   dataset (>= 0.026%, ~20x higher). A trust check that also bounds the
   proportion of contradicting/"suspicious" ambiguous indels would keep DecodeME
   on the trusted path and move the GWAS Catalog harmonised files to the
   untrusted path.

3. GWAS Catalog harmonised files reach 100% consistency through the catalog's
   own harmonisation, not through being raw reference-aligned analysis output.
   DecodeME's regenie output is the only tested table that is both trusted and
   raw.

## Open items

- Whether the trusted path should keep ambiguous indels in source orientation at
  all, versus applying frequency rules everywhere.
- Handling of ambiguous indels with 0/1 panel frequency, and whether to require
  both orientations in the panel (informed by the above).
- What to do when EAF is missing.

Source scripts and logs: v3_other_truth_sets.py / .log,
inspect_v3_wrong_choices.py / .log, v3_wrong_choices.md, v3_rule_variants.py /
.log, tune_ambiguous_indel_rules.py / .log.

## DecodeME preprint / data analysis plan evidence

We downloaded the DecodeME preprint (medRxiv 2025.08.06.25333109, via the
University of Edinburgh mirror) and the 2024 Data Analysis Plan, and read the
Methods with pdftotext. Extracts are in decodeme_preprint/*.txt (PDFs are
gitignored).

What the preprint says (main Methods):
- Cases genotyped on the UKB Axiom array; case and control genotypes merged and
  jointly imputed against "a reference panel from the whole-genome sequences of
  more than 200,000 gender-matched UKB samples".
- Association by REGENIE Firth logistic regression; variants on chr1-22 with
  INFO >= 0.90 and MAF >= 1%.
- Results are reported in GRCh38 coordinates.
- All imputation/QC detail is deferred to a Supplementary Methods file we could
  not retrieve (medRxiv is behind bot protection); the main text and DAP contain
  no statement on indel normalization or allele-orientation of the released
  summary statistics.

What the 2024 Data Analysis Plan says (superseded by the preprint pipeline):
- Planned a GRCh37 pipeline (HRC r1.1 + UK10K+1000Gp3, Sanger Imputation
  Service, PBWT), with pre-imputation QC by the HRC/1000G preparation tool that
  "checks the strand, alleles, position, ref/alt assignments ... [and] can be
  updated if discrepancy is [found]", dropping ambiguous A/T,G/C with MAF > 0.4.
- The released sumstats are GRCh38 from a UKB WGS panel, so the DAP's GRCh37
  HRC pipeline was not what produced them.

Bearing on the working assumption (DecodeME build 38 is reference-aligned):
- Supportive but indirect. Imputed-genotype summary statistics inherit the
  imputation panel's REF/ALT; the panel is a GRCh38 WGS panel, and DecodeME
  build 38 is 100% consistent with the hg38 FASTA on SNVs and indels. So its
  alleles follow a GRCh38 WGS-panel representation.
- The 11 suspicious ambiguous indels are most likely representation differences
  between the UKB WGS imputation panel and the 1000 Genomes EUR panel in repeat
  regions (different anchor/length for the same event), not source misalignment.
  This means even a perfectly aligned dataset yields a small non-zero suspicious
  rate, so the suspicious-indel threshold must be > 0 and is best expressed as a
  proportion.
- As anticipated, neither document discusses how ambiguous-indel orientation was
  chosen for the released summary statistics.
