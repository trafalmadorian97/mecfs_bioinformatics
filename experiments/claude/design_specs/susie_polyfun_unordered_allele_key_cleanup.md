# Spec: retire the unordered allele key from the SUSIE-polyfun path

Status: high-level, for discussion (2026-09-17). Detailed plan to follow.

## Problem

`unordered_allele_key(a, b)` sorts the two alleles and joins them, so `{A,B} == {B,A}`.
For SNVs that is a correct, orientation-agnostic variant key. For indels it is **wrong**:
mirrored pairs like `T/TCA` and `TCA/T` are distinct variants but collide on one key. It
was introduced for the PPP database (SNV-only by construction) and then reused in the
polyfun annotation path and the SUSIE prior join, where indels are present (~9.8% of the
baseline-LF / polyfun-prior variants). This can silently mis-match or collapse indels in
fine-mapping.

## Target invariant

Every input to a fine-mapping join is harmonized to a single genome reference (the hg19
FASTA), so the non-effect allele equals the reference allele. Variant identity is then the
ordered tuple **(chrom, pos, ea, nea)**, and every join uses it directly. No runtime
orientation reconciliation (unordered key, allele-flip harmonizer) is needed on this path.

## Verification (done, gates this spec)

Script: `experiments/claude/orientation_check/verify_fasta_orientation.py`. Against hg19:

- baseline-LF annotations: **A1 == REF 99.9995%** (19.48M rows), 9.77% indels.
- Broad UKBB LD panel (local sample, 131k): **allele1 (=NEA) == REF 100%**.
- polyfun prior / snpvar meta (one asset, feeds both the SUSIE prior join and ridge
  weights): **A1 == REF 99.9995%** (19.48M rows), 9.77% indels.
- PPP HapMap3 index: **0 indels / 1.20M rows (100% SNV)**.
- Incidental: prior has exactly **548 more rows** than the annotation matrix — the mirrored
  indel pairs the annotation build's unordered-key dedup currently drops.

Conclusion: the reference-orientation assumption holds for gwas, LD panel, annotations, and
prior, so switching these joins to the ordered 4-tuple is orientation-safe and strictly more
correct for indels.

## Changes (high level)

Grouped by concern, not yet sequenced.

### 1. Make the key honest about its domain
- `ppp_database/allele_key.py`: document that it is SNV-only / unsafe for indels; consider a
  cheap SNV guard so misuse fails loudly rather than by comment.

### 2. Guard the PPP infrastructure (unordered key is load-bearing there, SNV-only)
- `construct_ppp_variant_index_task.py`: assert every index variant is an SNV. No-op for the
  HapMap3 index (already 100% SNV); documents the invariant the whole PPP DB depends on.
- Decision needed for the common-1kg mode (contains indels): filter-then-assert, or keep the
  assertion HapMap3-only. (See open question C.)
- `common_1kg_membership_task.py`, `build_slim_protein_parquet_task.py`: same family; confirm
  they inherit the SNV guarantee or add it.

### 3. Restrict the allele harmonizer to its only correct use
- `harmonize_gwas_with_reference_table_via_chrom_pos_alleles.py`: assert all variants are
  SNVs (its palindrome logic is SNV-only anyway). This is the last legitimate caller once the
  generators stop using it.

### 4. Convert the baseline-LF annotation join to the ordered 4-tuple (3 files, coupled)
The annotation matrix is a **producer** with two consumers that also key on the unordered key;
all three change together:
- `build_baseline_lf_annotation_parquet_task.py` (producer): dedup / key on (CHR, BP, EA, NEA)
  instead of the unordered key. This is indel-aware and stops dropping the 548 mirrored pairs.
- `ridge_annotation_weights_task.py` (consumer): join annotations↔snpvar on the 4-tuple.
- `polyfun_explain_contrast_task.py` (consumer): join run-variants↔annotations on the 4-tuple.
- Orientation is safe (annotations and snpvar are both A1==REF; map A1→nea, A2→ea).
- Cost: rebuilds the 19.5M-row matrix + ridge weights (cache invalidation). See open question B.

### 5. Remove the unordered key from the SUSIE finemap task
- `susie_r_finemap_task.py`: the main gwas↔LD join already uses the 4-tuple; only the **prior
  join** uses the unordered key. Switch it to (chrom, pos, ea, nea). Prior NEA==REF makes this
  orientation-safe; the existing "prior does not cover N variants" check stays as the guard.

### 6. Drop the allele harmonizer from the fine-mapping generators
- `fine_mapping_asset_generator.py` and `polyfun_explain_fine_mapping_asset_generator.py`:
  stop wrapping the sumstats in `HarmonizeGWASWithReferenceViaAlleles`; rely on the sumstats
  already being FASTA-harmonized (NEA==REF) so `align_data`'s existing 4-tuple join aligns to
  the (also REF-oriented) LD panel. Depends on open question A.

### 7. Defensive check that inputs really are harmonized
- Protect against someone feeding un-harmonized sumstats once the reconciler is gone. Options
  in open question D: assert once in the harmonization task vs. a NEA-vs-FASTA check in SUSIE
  vs. a join-coverage floor in SUSIE. Intent: a missing upstream harmonization becomes a loud
  failure, not silent variant loss.

### 8. Document the remaining/for-future uses
- Record which unordered-key sites are safe-by-SNV (PPP infra) vs. still indel-exposed
  (common-1kg mode), so future work has the full map. The grep-complete site list is in the
  conversation / memory.

## Open questions to resolve before the detailed plan

- **A. Are all fine-mapping sumstats FASTA-harmonized?** Item 6 assumes every
  `build_37_sumstats_task` passed to the generators has gone through genome-reference
  harmonization (NEA==REF). Need to audit the callers; any locus fed raw sumstats would break.
- **B. Rebuild the annotation matrix indel-aware, or keep the SNV-collapse?** Exact-key dedup
  restores the 548 mirrored indels and invalidates the 19.5M-row matrix + ridge weights. Is the
  correctness worth the rebuild, or do we accept dropping mirrored indels on this path?
- **C. common-1kg PPP index:** filter-then-assert (drop its 258 indels) or leave the SNV
  assertion HapMap3-only?
- **D. Where does the defensive harmonization check live**, and is it an assert, a FASTA
  re-check, or a coverage floor? (SUSIE runs many times per locus, so cost matters.)
- **E. Scope of the SNV assert helper:** one shared `assert_all_snv` free function reused by
  items 2/3, per the helper-free-function convention?

## Out of scope
- The known L=2 credible-set bug (`fine_mapping_asset_generator.py`, separate memory).
- Any change to the PPP database results themselves (HapMap3 index is already SNV-clean).
