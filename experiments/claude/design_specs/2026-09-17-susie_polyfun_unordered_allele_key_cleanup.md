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

## Decisions (resolved 2026-09-18)

- **A. All fine-mapping sumstats are FASTA-harmonized.** Traced: every DecodeME finemap locus
  (with/without palindromes and polyfun-explain) feeds
  `DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.join_task`. Its chain is
  `GenomeReferenceHarmonizationTask` (NEA==REF vs hg19 FASTA + 1kg EUR panel) ->
  `JoinDataFramesTask` (dbSNP150 rsid assignment, inner-join on CHR/POS/EA/NEA,
  orientation-preserving) -> `.join_task`. DecodeME is the only fine-mapped trait. Item 6 is safe.
- **B. Rebuild the annotation matrix indel-aware** (exact `(CHR,BP,EA,NEA)` key). Restores the
  548 mirrored-indel pairs; invalidates the 19.5M-row matrix + ridge weights (accepted cost).
- **C. Hard SNV assertion on the PPP index.** HapMap3 index is already 100% SNV (verified), so it
  passes with no data change. The common-1kg mode is not used now; revisit (filter-then-assert)
  when it is.
- **D. Shift-left harmonization provenance in metadata** (see below). The construction-time check
  verifies build equality; an optional execute()/audit data check can use the column fields.
- **E. One shared `assert_all_snv` free function**, reused by the PPP-index and harmonizer guards.

### D in detail: harmonization provenance types

New value type:

    @frozen
    class HarmonizationInfo:
        build: GenomeBuild          # "19" | "38"
        ref_allele_col: str         # the column asserted to equal the reference base
        pos_col: str                # the 1-based position column

Meta changes (all SIBLING leaves under the existing abstract `FileMeta`/`DirMeta` -- no
subclassing of the concrete reference metas, so no concrete class becomes a non-leaf):

- `FilteredGWASDataMeta` gains `harmonization_info: HarmonizationInfo | None = None` (a
  processed-GWAS asset can always be meaningfully asked "harmonized to which build?").
- `FASTAMeta(DirMeta)`: sibling of `ReferenceDataDirectoryMeta`, adds non-optional
  `build: GenomeBuild`. Used by `IndexedFastaTask`. `simple_meta_to_path` gets a branch
  returning the same relative path as `ReferenceDataDirectoryMeta`.
- `HarmonizableReferenceTableMeta(FileMeta)`: sibling of `ReferenceFileMeta`, adds
  `harmonization_info: HarmonizationInfo | None`. Used by the Broad LD-label task and the 1kg
  EUR panel task. `simple_meta_to_path` gets a branch returning the same relative path as
  `ReferenceFileMeta`.

Because these are siblings (not subclasses), the ~21 generic transformer `.create()` methods
that check `isinstance(source_meta, ReferenceFileMeta)` do NOT capture
`HarmonizableReferenceTableMeta` -> piping one through a generic transformer fails closed
(hits its `else: raise`). That is the desired safety and the reason `RenameColsTask` exists.

`RenameColsTask` (new; capability-restricted vs `PipeDataFrameTask` so it can update meta
accurately): applies a column rename map. Its `create()` derives the output meta from the
source; when the source is `HarmonizableReferenceTableMeta` (or a `FilteredGWASDataMeta`) with
`harmonization_info`, it remaps `ref_allele_col`/`pos_col` through the rename map (asserting
neither key is dropped/collided) and keeps `build`. Non-provenance sources fall back to standard
meta derivation. Replaces the `PipeDataFrameTask` used to rename the Broad LD labels.

Flow:
1. `UCSC_HG19_INDEXED_FASTA` -> `FASTAMeta(build="19")`; the 1kg EUR panel task and the raw Broad
   LD-label download -> `HarmonizableReferenceTableMeta(harmonization_info=...)` (build "19",
   ref/pos cols = that asset's reference-allele and position columns).
2. `GenomeReferenceHarmonizationTask.create` asserts `fasta.build == panel.build` and stamps its
   output `FilteredGWASDataMeta.harmonization_info = HarmonizationInfo(build, ref_allele_col="NEA",
   pos_col="POS")`.
3. Spot-2 propagation (the DecodeME path): `JoinDataFramesTask` (rsid assignment) is
   orientation-preserving, so its `create()` copies `harmonization_info` forward from the
   harmonized source. It is the ONLY intermediate task on the DecodeME finemap path.
4. LD labels are renamed via `RenameColsTask`, so the renamed labels task carries an updated
   `HarmonizableReferenceTableMeta`.
5. `SusieRFinemapTask.create` asserts (meta-to-meta, construction time):
   `gwas.harmonization_info is not None` and
   `gwas.harmonization_info.build == ld_labels.harmonization_info.build`.

## Decomposition into two plans

The work splits into two independently shippable subsystems (each its own plan):

- **Plan 1 -- Harmonization provenance plumbing (this spec's D, minus the SUSIE check).** Adds the
  types, `RenameColsTask`, asset tagging, stamping, spot-2 propagation, and LD-label tagging.
  Changes no results and enforces nothing yet; it makes `harmonization_info` reach the finemap
  generator input and the renamed LD labels. Independently testable via the asset metas.
- **Plan 2 -- Retire the unordered allele key + enforce the guard (items 1-3, 5, 6, 8).**
  `assert_all_snv`, the allele_key doc, PPP-index + harmonizer SNV asserts, the indel-aware
  annotation rebuild (producer + ridge weights + explain contrast), the SUSIE prior exact join,
  removing `HarmonizeGWASWithReferenceViaAlleles` from the two generators, AND the SUSIE
  build-match assertion. The assertion belongs here because `HarmonizeGWASWithReferenceViaAlleles`
  drops `harmonization_info`; only once it is removed (item 6) does the flag reach SUSIE.

Order: Plan 1 before Plan 2. Flow step 5 (the SUSIE `create()` assertion) is implemented in Plan 2.

## Out of scope
- The L=2 credible-set bug: VERIFIED FIXED 2026-09-18 (`fine_mapping_asset_generator.py` L=2 task
  uses `max_credible_sets=2`); stale memory deleted.
- Any change to the PPP database results themselves (HapMap3 index is already SNV-clean).
