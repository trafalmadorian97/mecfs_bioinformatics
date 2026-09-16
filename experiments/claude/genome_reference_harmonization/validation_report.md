# Genome-reference harmonization: validation report

Results of the validation experiments (V1-V5) for the in-repo genome-reference
harmonization Task. See the design spec and plan in
experiments/claude/design_specs/2026-09-14-genome-reference-harmonization-*.md.

## V3: stringent indel rules and the suspicious-indel trust gate

The original V3 goal -- a (distance, margin) grid cell with zero wrong swaps on
build-38 DecodeME -- is unreachable. No cell gives zero wrong, because a trusted,
reference-aligned table still carries a handful of ambiguous indels that differ
only in panel representation (UKB WGS vs 1000 Genomes) in repeat regions. That
finding (v3_rule_variants.py, inspect_v3_wrong_choices.py, v3_wrong_choices.md)
and the cross-dataset investigation (v3_other_truth_sets.py,
v3_trust_investigation_summary.md) motivated the suspicious-indel trust gate.

### Untrusted-path resolution defaults

Set to the most conservative sensible grid cell (fewest wrong swaps):

- indel_max_af_distance = 0.02
- indel_min_af_margin = 0.3

Untrusted tables are exactly those whose orientation is unreliable, so aggressive
dropping is appropriate.

### Suspicious-indel gate calibration (v3_suspicious_calibration.py)

Run of the production count_trust_evidence_genome_wide (with the monomorphic-panel
filter and the checkable denominator) against the hg38 FASTA and the 1000 Genomes
EUR hg38 panel, at the 0.02/0.3 resolution above. "checkable" = ambiguous (BOTH)
indels with EAF and a panel record after dropping monomorphic (AF 0 or 1) records;
"suspicious" = those decide_ambiguous_indels resolves to a swap.

| Dataset | suspicious | checkable | fraction | trusted @ 1e-4 |
|---|---|---|---|---|
| DecodeME build 38 (truth set) | 5 | 734,721 | 6.8e-6 | yes |
| MVP myocardial infarction | 242 | 744,851 | 3.25e-4 | no |
| Bellenguez Alzheimer's | 20,160 | 772,074 | 2.61e-2 | no |
| Kerrebijn fibromyalgia | 131,943 | 1,383,960 | 9.53e-2 | no |

All four pass 100% SNV+indel consistency, so the older trust check trusted all of
them. With the gate, DecodeME (raw imputed output) stays trusted and the three
GWAS Catalog harmonised files move to the untrusted path.

### Chosen defaults (confirmed)

- max_suspicious_indel_fraction = 1e-4 -- ~48x above DecodeME, ~3.25x below the
  nearest untrusted dataset.
- min_checkable_ambiguous_indels = 100 -- every dataset has far more than 100
  checkable ambiguous indels, so the noise floor never suppresses these.

The monomorphic-panel filter lowered DecodeME's suspicious count from 11 to 5
(and left MVP above the threshold), as anticipated when the filter was adopted.

## V1, V2: gwaslab parity (DecodeME build 37, Liu et al. 2023 IBD)

Pending (compare_with_gwaslab.py).

## V4: peak memory

Pending (memory_benchmark.py / measure_harmonization_memory.py).

## V5: trust decision across rsID-assignment chains

Pending (survey_trust.py, Task 12).
