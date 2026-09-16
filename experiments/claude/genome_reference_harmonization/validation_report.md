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

compare_with_gwaslab.py joins the new harmonized output with the cached
gwaslab-harmonized table on SNPID and buckets every row. Both tables are
untrusted (build-37 liftover leaves reference-difference blocks: DecodeME 17,092
inconsistent SNVs, Liu 15,339), so both take the stringent path.

| bucket | DecodeME build 37 | Liu 2023 IBD |
|---|---|---|
| identical -- SNV | 6,770,430 | 7,841,922 |
| identical -- palindromic SNV | 1,062,923 | 1,225,928 |
| identical -- indel | 502,156 | 316,386 |
| opposite_orientation (indel) | 128 | 8 |
| new_only (indel) | 36,020 | 18,016 |
| gwaslab_only (indel) | 109,937 | 2,547 |
| same_orientation_beta_differs | 0 | 0 |
| other | 0 | 0 |

Every difference is explained:

- **opposite_orientation**: every one of these indels carries gwaslab STATUS digit
  7 = 4, the indel-inference flip signature -- i.e. exactly the gwaslab bug this
  Task replaces. We orient them by the FASTA; gwaslab flipped them to a
  different variant at the repeat. This is the intended correction.
- **gwaslab_only** (rows gwaslab kept, we dropped), by our drop reason:
  DecodeME 108,249 ambiguous_indel_af_mismatch, 723 af_indecisive, 498
  indel_not_on_reference, 467 not_in_panel; Liu 1,612 / 5 / 179 / 751. These are
  the stringent ambiguous-indel rules dropping in the safe direction, plus
  reverse-strand indels we do not try to rescue.
- **new_only**: ambiguous indels we keep whose gwaslab STATUS put them in the
  drop set.
- **same_orientation_beta_differs = 0 and other = 0**: no unexplained orientation
  or statistic differences. (A1FREQ_CASES/CONTROLS, which gwaslab does not flip
  but we do, do not change orientation, so they stay in identical.)

Conclusion: genome-reference harmonization reproduces gwaslab's orientation
everywhere except the indel-inference bug, and its extra drops are the intended
stringent ambiguous-indel handling.

## V4: peak memory

memory_benchmark.py runs each scenario in its own process on DecodeME build 37
and samples RssAnon (anonymous memory; the memory-mapped FASTA is file-backed and
excluded). The acceptance criterion is that peak memory follows the largest
chromosome, not total rows: all-chromosomes peak <= 1.5x the chr1-and-2 peak.

| scenario | peak RssAnon (GiB) | wall time (s) |
|---|---|---|
| genome_reference chr1-2 | 1.56 | 6.1 |
| genome_reference all | 1.58 | 36.2 |
| gwaslab | 11.11 | 927.7 |

- all / chr1-2 peak ratio = 1.01 (<= 1.5): peak memory is bounded by the largest
  chromosome slice, independent of total rows. The OOM-safety goal is met.
- genome-reference harmonization uses ~7x less peak memory than gwaslab
  harmonization (1.58 vs 11.11 GiB) and runs ~25x faster (36 vs 928 s) on the
  same input.

## V5: trust decision across rsID-assignment chains

survey_trust.py imports the asset package, collects RSIDAssignmentTaskGroup
instances, and reports the trust decision of each group's pre-harmonization table
whose upstream gwaslab pickle is materialized. Only the DecodeME chains build
their groups through annovar_37_basic_rsid_assignment (the Liu chain constructs
its harmonize task directly, so it is not a group); two unrelated asset modules
fail to import at construction time and are skipped
(decode_me_filtered_gene_list_with_gene_metadata_drop_cols and
partitioned_model_allele_freq -- pre-existing, not touched by this work).

| chain | trusted | inconsistent SNVs | suspicious / checkable |
|---|---|---|---|
| decode_me_gwas_1_genome_reference_harmonized | no | 17,092 | 660 / 598,255 |
| decode_me_gwas_1_keep_ambiguous_...          | no | 17,092 | 660 / 598,255 |

Both DecodeME build-37 chains are untrusted, driven by the 17,092 inconsistent
SNVs in the GRCh37-vs-GRCh38 reference-difference blocks (the suspicious gate is
moot once consistency fails). So production rsID assignment runs on the stringent
path, which is the safe direction for a lifted table. Note this is the build-37
lifted DecodeME; the build-38 DecodeME used to calibrate the suspicious gate (V3)
is trusted.
