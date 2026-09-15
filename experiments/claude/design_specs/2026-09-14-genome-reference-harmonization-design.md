# In-Repo Genome-Reference Harmonization --- Design

Date: 2026-09-14

## Goal

Replace gwaslab harmonization (gwaslab's Sumstats.harmonize, called through
HarmonizationOptions) with a Task implemented in this repo in polars. It
orients each variant's alleles and allele-specific statistics against a
genome FASTA, and resolves strand for palindromic SNVs with a reference
panel. It handles ambiguous indels with a conservative scheme that never
keeps an indel whose orientation it cannot establish.

Secondary goal: harmonization must not risk OOM on large summary statistics.
Peak memory is bounded by one chromosome of the input, not the whole table.

## Terminology

Three steps now carry the word "harmonization". Name them explicitly in code,
docstrings and prose:

- **gwaslab harmonization**: gwaslab's Sumstats.harmonize, which this design
  replaces.
- **genome-reference harmonization**: the new Task described here, which
  orients sumstats against a genome FASTA plus a reference panel.
- **task-reference harmonization**: HarmonizeGWASWithReferenceViaAlleles,
  which is unchanged and out of scope.

## Background

Measured findings are recorded in experiments/claude/polyfun_annot_dedup/.

### The gwaslab indel bug

Both alleles of an indel "match the genome" whenever the longer allele is
present in the reference. That holds for every deletion, and for insertions at
repeats. The short allele is a prefix of the long one, so it always matches.
gwaslab gives these rows STATUS digit 6 = 6.

Its indel-inference step (check_unkonwn_indel, hm/hm_harmonize_sumstats.py,
4.2.2 L1592-1633) then scans 1000 Genomes EUR VCF records at the position. It
flips the row on the first record whose REF and ALT are the string-swap of the
row's alleles, within a 0.2 allele-frequency tolerance. That record describes
a different variant: an insertion versus a deletion of the same bases at a
repeat.

- 1,598 DecodeME indels were relabelled this way.
- 176 of those matched a monomorphic (AF = 0) record.
- A one-row reproduction is in gwaslab_indel_flip_repro.py. The issue has
  been reported upstream.

### Mirrored indel pairs

UKB-derived tables contain distinct variants at the same CHR/POS with swapped
alleles (T/TCA and TCA/T). Mislabelling one as the other silently corrupts any
downstream allele-keyed join. See project memory
project_ukb_mirrored_indel_pairs.

### Trust signal (decodeme_pre_harmonization_trust_signal.py)

This is the fraction of non-palindromic SNVs whose NEA equals the FASTA
reference base.

| Table | NEA = REF | Swapped | Reverse-strand matches |
|---|---|---|---|
| DecodeME build 38 vs hg38 | 100.000% | 0 | 0 |
| DecodeME lifted to 37 vs hg19 | 99.748% | 0.160% | 0.092% |

Liftover never changes alleles. The build-37 mismatches are blocks where the
GRCh37 and GRCh38 references differ: 1.0% on chr22, 0.8% on chr19. So indels in
those blocks may carry the wrong orientation, and nothing short of the rules
below can safely identify them.

### What gwaslab harmonization does for us today

Two call sites use HarmonizationOptions:

- the annovar build-37 rsID-assignment generator, which every concrete
  standard analysis goes through, DecodeME included;
- Liu et al. 2023 IBD.

A system test, test_harmonize_drop_ambiguous.py, also uses it. Both production
call sites dump the result straight to parquet. Nothing reads STATUS except
our own drop logic in _do_harmonization.

The steps we rely on, as implemented in gwaslab 4.2.2, are these.

**1. Reference check (check_ref, _fast_check_status).** Alleles are compared
as prefixes of the reference starting at POS, case-insensitively. The
digit-6 outcomes are:

| Digit 6 | Condition |
|---|---|
| 0 | NEA matches, EA does not |
| 3 | EA matches, NEA does not |
| 6 | Both match and the lengths differ |
| 4 | Reverse complement of NEA matches |
| 5 | Reverse complement of EA matches |
| 8 | Nothing matches |

Gaps and quirks:

- Rows on contigs absent from the FASTA stay at 9 and survive the drop filter.
- Equal-length rows where both alleles match fall through the conditions.

**2. Statistic flipping (_flip_allele_stats, qc/qc_fix_sumstats.py
L1655-1795).**

- Digit 6 = 4 or 5: reverse-complement EA and NEA.
- Digit 6 = 3 or 5: swap EA and NEA, negate BETA, Z, T and the BETA CI, set
  EAF to 1 - EAF, invert OR and HR with their CIs, flip DIRECTION.
- Indel flip: pattern xxxx[123][67]6, the bug above.
- Palindromic minus strand (digit 7 = 5): flip the statistics only, leave the
  alleles.

Gaps and quirks:

- A1FREQ_CASES, A1FREQ_CONTROLS and NEAF are not flipped.
- BETA CI bounds are flipped wrongly. flip_by_sign assigns
  BETA_95U = -BETA_95L, then BETA_95L = -BETA_95U, reading the value it just
  overwrote, so BETA_95L comes back unchanged. Verified by reading the
  installed 4.2.2 source and gwaslab main on 2026-09-14.
  - This is the BETA analogue of the OR bug reported in gwaslab issue #209.
  - The OR and HR branches of flip_by_inverse now save both originals before
    assigning, so the OR case is fixed in code. The issue itself is still open
    on GitHub.
  - The BETA branch was not changed. Worth a follow-up comment on #209.

**3. Palindromic strand inference (_parallelize_infer_strand, mode "p").**
This applies only to standardized, normalized SNVs where {EA, NEA} is {A, T}
or {C, G}.

| Digit 7 | Condition |
|---|---|
| 7 | EAF missing, or 0.4 < EAF < 0.6 |
| 8 | No record with REF = NEA and ALT = EA, or reference MAF > 0.4 |
| 1 (keep) | EAF and AF are on the same side of 0.5 |
| 5 (flip) | Otherwise |

**4. Indel inference (mode "i").** This is the bug.

**5. Our drops.**

- drop_missing_from_ref_seq removes digit 6 = 8.
- drop_missing_from_ref_infer_or_ambiguous removes digit 7 in {7, 8}. For
  indels whose alleles both match, that means no panel record or MAF > 0.4.

Not used, so not reimplemented: rsID assignment, duplicate removal, sweep mode.

**Stays upstream in gwaslab** (GWASLabCreateSumstatsTask): basic_check and
liftover.

- basic_check includes allele normalization. That is parsimony trimming of
  shared prefix and suffix; it neither left-aligns nor uses the FASTA.
- Liftover does not change alleles.

## Design

### Components

1. **Pinned reference assets**, following the eur_hg38_30x_vcf.py precedent.
   - HG19_UCSC_FASTA_GZ: a DownloadFileTask of
     hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz, pinned by
     md5. The local copy is 806c02398f5ac5da8ffd6da2d1d5d1a9; confirm it
     against a fresh download.
   - HG19_UCSC_FASTA: an IndexedFastaTask that decompresses by streaming and
     writes the .fai with pysam.faidx. Output is a directory asset holding
     genome.fa and genome.fa.fai. pysam is already a pinned dependency.
   - THOUSAND_GENOMES_EUR_HG19_VCF: a DownloadFileTask of gwaslab's
     1kg_eur_hg19 URL, pinned by md5 2c78cb84cb1f90b576510decc45e5b9b. gwaslab
     has republished this panel in place before, so the checksum is what makes
     the content immutable. Rehost it if the URL changes.
   - ReferencePanelAlleleFrequencyTask converts a panel VCF to a parquet asset
     with columns CHR (gwaslab integer coding), POS, REF, ALT, AF.
     - Implementation: bcftools view -G (drop the 503 genotype columns), then
       bcftools query, as in FilterCommon1kgVariantsTask.
     - One row per (CHR, POS, REF, ALT). The panel is already split and
       normalized.
     - Exact duplicate keys with identical AF collapse to one row. Duplicate
       keys with differing AF are removed and counted in the log, so a lookup
       never sees two answers.
     - Rows are written sorted by CHR, POS so parquet row-group statistics let
       a per-chromosome filter read only that chromosome.
     - The hg19 panel indexes 1-22 and X only.
     - It also works on the hg38 30x panel, which the tuning experiment
       (V3 below) needs.
   - **Disk: store only derived assets.** The downloads are large: the hg19
     panel VCF is 3.5 GB (503 genotype columns) and hg19.fa.gz is 0.95 GB.
     - The panel parquet task is wrapped in DiscardDepsWrapper, so the VCF is
       materialized in a temporary directory and only the parquet (sites and AF
       only, far smaller) is kept.
     - IndexedFastaTask is wrapped the same way, so only the uncompressed
       genome.fa (3.2 GB, which memory-mapping requires) and its .fai are kept.
     - The wrapper appends "_wrapped_remove_deps" to the asset id and has no
       deps. A rebuild of the wrapped asset re-downloads the source, which is
       acceptable for reference data that is built once.

2. **Library module** mecfs_bio/build_system/task/genome_reference_harmonization/.
   These are pure functions on polars frames plus a FASTA reader; no gwaslab
   imports.
   - fasta.py
     - IndexedFasta: a memory-mapped uncompressed FASTA plus its .fai.
     - reference_matches(chrom, pos, allele): a boolean array, gathered
       vectorized per allele-length group.
     - This is the prototype in benchmark_fasta_ref_check.py: 4.6 s for 19.5M
       rows.
     - A group whose (rows x length) gather would exceed a byte budget is
       processed in slices. Very long alleles then never materialize a large
       matrix; this is the root of the gwaslab long-allele OOM.
   - classify.py: per-row allele classes (below).
   - trust.py: the trust decision from class counts.
   - indel_rules.py: trusted and stringent resolution of ambiguous indels.
   - palindromes.py: palindromic strand inference.
   - flip.py: the column flip registry and its application.

3. **Task** GenomeReferenceHarmonizationTask.
   - Dependencies: sumstats (any tabular asset with gwaslab column names,
     read with scan_dataframe_asset), fasta, reference panel.
   - Field pipe: DataProcessingPipe = IdentityPipe(), applied to the scanned
     narwhals LazyFrame before harmonization, following repo convention.
   - Output: one parquet file with FilteredGWASDataMeta, derived in .create()
     from the sumstats task's meta.

### Preconditions (asserted)

Violating any of these fails the Task with a message naming the offending
column or value.

- **Columns:** CHR (integer coding, X = 23, Y = 24, MT = 25), POS, EA, NEA are
  present. EAF is optional.
- **Alleles:** uppercase-able strings over ACGT, EA != NEA, and parsimonious
  (upstream basic_check normalization).
- **Contigs:** the fasta and panel dependencies cover every chromosome in the
  input.
  - Exception: rows on contigs listed in options.excluded_chromosomes are
    dropped with a reason, not asserted.
  - The default is {25}. UCSC hg19 chrM (16,571 bp) is not rCRS, so MT cannot
    be oriented against it.
- **Other columns:** every non-allele column has a registered flip rule or is
  registered invariant (flip registry below). An unregistered column fails
  loudly so a new statistic is never silently left unflipped.

### Per-row allele classes

These are computed on the plus strand, plus reverse complement where
meaningful. Equal-length variants (SNVs, MNPs) and different-length variants
(indels) are classified separately.

**Equal length (SNV or MNP):**

| Class | Condition |
|---|---|
| NEA_REF | NEA matches |
| EA_REF | EA matches |
| NEA_REF_RC | Reverse complement of NEA matches |
| EA_REF_RC | Reverse complement of EA matches |
| NONE | Nothing matches |

- Checks are made in that order: the reverse complement is only consulted when
  neither plus-strand allele matches, as in gwaslab.
- PALINDROMIC_SNV is a flag (EA is the complement of NEA) carried alongside the
  class.
- PALINDROMIC_MNP is the MNP analogue: the reverse complement of NEA equals
  EA. Strand cannot be inferred from the FASTA, and there is no frequency
  inference for MNPs. Handling mirrors palindromic SNVs in trusted mode:
  - **trusted:** oriented by the plus-strand class (NEA_REF keep, EA_REF swap)
    and kept, because a 100%-consistent source is plus-strand aligned;
  - **untrusted:** dropped (reason palindromic_mnp_untrusted).
  - Non-palindromic MNPs are oriented by the FASTA in both modes, with no
    strand inference, as gwaslab does.

**Different length (indel):**

| Class | Condition |
|---|---|
| NEA_ONLY | NEA matches, EA does not |
| EA_ONLY | EA matches, NEA does not |
| BOTH | Both match; ambiguous |
| NONE | Neither matches |

Reverse complements are deliberately not tried for indels. A reverse-complemented
left-anchored indel does not keep its anchor base, so gwaslab's digit 4 and 5
matches for indels are near-certainly spurious. This deviates from gwaslab; V1
below quantifies it.

### Trust decision (pass 1)

Computed over the whole input before any row is resolved.

- **Checkable SNVs:** non-palindromic SNVs on non-excluded chromosomes.
  - Consistent: class NEA_REF.
  - Inconsistent: any other class, including reverse-strand matches and NONE.
- **Checkable indels:** all indels.
  - Consistent: NEA_ONLY.
  - Inconsistent: EA_ONLY or NONE.
  - BOTH is excluded from the count.
- **trusted** is true iff all of the following hold:
  - checkable SNVs >= options.min_checkable_snvs (default 10,000), and all are
    consistent;
  - checkable indels >= options.min_checkable_indels (default 1,000), and all
    are consistent.

The bar is exactly 100% (user decision, 2026-09-14). Keeping one wrong
indel costs more than dropping many correct ones. A single stray mismatch
anywhere sends the table down the stringent path, which is the safe direction.

The trust counts, the decision and the per-chromosome inconsistency rates are
logged.

### Row resolution (pass 2)

**Equal-length variants:**

| Class | Action |
|---|---|
| NEA_REF | Keep |
| EA_REF | Swap alleles, flip statistics |
| NEA_REF_RC | Complement both alleles, statistics unchanged |
| EA_REF_RC | Complement and swap alleles, flip statistics |
| NONE | Drop (not_on_reference) |

**Palindromic SNVs, after the orientation above:**

Trusted, and the class is NEA_REF or EA_REF: keep the source strand. Skip
frequency-based strand inference.

- A 100%-consistent source is plus-strand aligned.
- Frequency-based inference against a EUR panel only adds risk there.

This deviates from gwaslab (confirmed by the user, 2026-09-14).

Untrusted: apply gwaslab's rule. Rows are oriented so EA = ALT and NEA = REF.

- **Unresolved** if any of these hold:
  - EAF is missing, or 0.4 < EAF < 0.6;
  - no panel record with REF = NEA and ALT = EA;
  - panel MAF > 0.4.
- **Keep** if EAF and AF are on the same side of 0.5.
- **Flip** otherwise: flip statistics, leave alleles unchanged.
- **Unresolved rows** are dropped (reason palindrome_unresolved), unless
  options.keep_unresolved_palindromes is set.
  - The option exists for the DecodeME with-palindromes runs.
  - It never applies to indels.

The 0.4 thresholds are options with gwaslab's defaults.

**Indels:**

| Class | Action |
|---|---|
| NEA_ONLY | Keep |
| EA_ONLY | Swap and flip |
| NONE | Drop (indel_not_on_reference) |
| BOTH, trusted | Keep the source orientation |
| BOTH, untrusted | Stringent rules below |

**Stringent rules for an ambiguous indel** with alleles (EA, NEA) and EAF. Two
candidate readings:

- **K (keep):** the variant is REF = NEA, ALT = EA. Its panel record has
  REF = NEA, ALT = EA. Predicted EAF = AF_K. Distance d_K = |EAF - AF_K|.
- **F (flip):** the variant is REF = EA, ALT = NEA. Its panel record has
  REF = EA, ALT = NEA. The row's EA is then the reference allele, so predicted
  EAF = 1 - AF_F. Distance d_F = |EAF - (1 - AF_F)|.

Rules, in order:

1. EAF missing: drop (ambiguous_indel_no_eaf).
2. Neither record in the panel: drop (ambiguous_indel_not_in_panel).
3. Exactly one record X:
   - choose X if d_X <= options.indel_max_af_distance;
   - otherwise drop (ambiguous_indel_af_mismatch).
4. Both records:
   - choose X if d_X <= options.indel_max_af_distance and
     d_other - d_X >= options.indel_min_af_margin;
   - otherwise drop (ambiguous_indel_af_indecisive).
5. Choosing F means swap and flip.

Defaults are set by V3 below; placeholders are 0.1 and 0.2. The rules drop in
the safe direction for non-EUR sumstats, where panel frequencies fit worse.

Worked example: the chr15:54698192 reproduction row in source orientation
(EA = T, NEA = TG, EAF = 0.845; the true variant is the deletion TG -> T).

- K reads the row as REF = TG, ALT = T. That record has AF 0.862, so
  d_K = |0.845 - 0.862| = 0.017.
- F reads the row as REF = T, ALT = TG. That record has AF 0.0, so the
  predicted EAF is 1.0 and d_F = 0.155.
- The margin is 0.138.
  - With margin 0.1, K is chosen: correct.
  - With the 0.2 placeholder, the row is dropped as indecisive: safe, but it
    loses a correct variant. V3 exists to settle this trade-off.

gwaslab's output row (EA = TG, NEA = T, EAF = 0.155) is the same variant in the
other orientation:

- K reads it as T -> TG: AF 0.0, so d_K = 0.155.
- F reads it as TG -> T: AF 0.862, so the predicted EAF is 0.138 and
  d_F = 0.017.

The distances simply swap, so F (flip back) is chosen under the same thresholds.
The rules are symmetric in orientation: the Task gives the same final variant
whichever orientation arrives.

### Flip registry (flip.py)

A module-level mapping from column name to a FlipRule. Every column name is a constant from gwaslab_constants.py or regenie_constants.py.

The registry covers two sets of columns:
- gwaslab's standard column names (the "gwaslab" entry of formatbook.json);
- the non-standard columns this repo's datasets carry: regenie binary-trait columns and N_EFF.

| Rule | Columns |
|---|---|
| negate | BETA, Z, T |
| complement (1 - x) | EAF, A1FREQ_CASES, A1FREQ_CONTROLS |
| invert (1 / x) | OR, HR |
| inverted bound pair | (OR_95L, OR_95U), (HR_95L, HR_95U): each bound becomes 1 / the other |
| invariant | SNPID, rsID, CHR, POS, SE, P, MLOG10P, CHISQ, F, P_HET, I2, SNPR2, DOF, N, N_CASE, N_CONTROL, N_CASES, N_CONTROLS, N_EFF, INFO, MAF, TEST, EXTRA |

Not registered (decided in plan review, 2026-09-15):
- **DIRECTION.** Only metal, mrmega and the auto formats map to it. This repo uses none of them, so an input carrying it fails the unregistered-column assertion.
- **NEAF and BETA_95L/BETA_95U.** Neither is a gwaslab standard column, and neither is used here.
- **REF and ALT.** After harmonization NEA and EA carry the reference orientation, which a source REF/ALT pair could contradict, so an input carrying them fails the unregistered-column assertion.

- STATUS is dropped from the output; gwaslab regenerates it when a later task builds a Sumstats object from the table.
- SNPID is left as-is, matching gwaslab.
- Task option extra_column_rules accepts additional name -> FlipRule entries,
  so a dataset-specific column can be declared without editing the module.
- The palindromic minus-strand flip applies every rule except allele swap.

### Uniqueness

After resolution, the Task asserts (CHR, POS, EA, NEA) is unique within each
chromosome. A mis-resolved mirrored pair would violate it.

### Streaming, per-chromosome execution

Input is pipe.process(scan_dataframe_asset(...)).to_native().

DataProcessingPipes may change a narwhals LazyFrame's backend, so to_native()
is not guaranteed to return polars. to_polars() would force collection, which
defeats streaming. The Task therefore asserts the result is a pl.LazyFrame,
with a message naming the backend actually received.

The list of chromosomes comes from a streaming unique over CHR.

**Pass 1 (trust).** For each chromosome:

- collect only CHR, POS, EA, NEA for that chromosome;
- classify;
- accumulate counts.

Peak memory is 4 columns of the largest chromosome. The FASTA is memory-mapped;
its pages are reclaimable page cache, not anonymous memory.

**Pass 2 (resolve).** For each chromosome:

- collect all columns for that chromosome;
- when untrusted, read only that chromosome's panel rows at positions
  holding a palindromic SNV or a BOTH indel in the slice. When trusted, the
  panel is not read at all, since trusted mode keeps palindromes and BOTH
  indels as-is;
- classify, resolve, flip, sort by POS;
- write scratch/parts/chr{CHR}.parquet;
- release the slice.

**Finish.** pl.scan_parquet(parts in CHR order).sink_parquet(output), which
streams.

Memory:

- Peak anonymous memory is roughly one chromosome of the input plus one
  chromosome of filtered panel rows plus the largest allele-length gather
  slice. Chromosome 1 or 2 is about 8% of a genome-wide table.
- Pass 1 re-reads the key columns, which costs one extra sequential scan.
  That is cheap for parquet (column projection) and slow but bounded for CSV.
  Large text inputs should be converted first (WhitespaceSepTextToParquetTask,
  per memory note).
- Parquet filter pushdown makes a per-chromosome collect read only the matching
  row groups when the input is CHR-sorted. gwaslab dumps are sorted. Unsorted
  inputs still stream through the filter with bounded memory.
- Chromosomes run sequentially. A process pool would multiply peak memory by
  the worker count, which defeats the purpose. Not planned.

### Logging

For each chromosome and in total:

- class counts;
- trust counts and decision;
- rows kept, swapped, complemented, palindrome-flipped;
- rows dropped by reason.

Drop reasons exist as a column inside the library functions, so experiments can
inspect every dropped row. The Task filters them out before writing.

## Validation (experiments, before switching call sites)

Scripts go in experiments/claude/genome_reference_harmonization/ with logs
captured via tee.

**V1: gwaslab parity, DecodeME build 37.**

- Inputs: the liftover-to-37 table, run through
  - the new Task;
  - gwaslab harmonization default (the cached decode_me_gwas_1__harmonized
    dump);
  - gwaslab harmonization with infer_strand mode "p".
- Join on SNPID and bucket every row:
  - identical;
  - dropped by one side only (by reason);
  - orientation differs;
  - statistics differ other than by a flip.
- Expected differences, each explained and counted:
  - ambiguous indels (stringent path, since trust is 99.83%);
  - reverse-strand indels;
  - A1FREQ_CASES/CONTROLS flipped;
  - unknown-contig rows gwaslab kept at digit 6 = 9.
- Any other bucket is a missed gwaslab behaviour to copy or consciously reject.

**V2: gwaslab parity, Liu et al. 2023 IBD.** Same as V1, and it exercises OR
columns.

**V3: tuning the stringent indel rules on a truth set.**

- Data: DecodeME build 38 against the hg38 FASTA and the pinned hg38 30x EUR
  panel. Trust is 100% and the source orientation is ground truth.
- Call indel_rules on every BOTH indel directly as if untrusted, over a grid of
  (indel_max_af_distance, indel_min_af_margin).
- Report per grid cell:
  - wrong choices (chose F for a trusted row);
  - drops by reason;
  - kept.
- Choose defaults with zero wrong choices that keep the most indels. Record the
  grid table in the log.
- The build-37 source cannot be used as truth: its reference-difference blocks
  contaminate it.

**V4: memory.**

- Run the Task on the largest local sumstats table (at least DecodeME, 9M
  rows) and, if available, a 20M+ row table.
- Sample RssAnon from /proc/self/status in a background thread and report its
  peak, alongside wall time.
- Compare against gwaslab harmonization's peak on the same input.
- Acceptance: peak RssAnon grows with the largest chromosome slice, not with
  total rows (compare the two inputs).

**V5: trust decision on other sources.** Run pass 1 alone on the concrete
standard analysis inputs. List each source's trust counts and decision, so we
know which sources get the looser path.

## Unit tests (Task level)

Use a synthetic FASTA (a few contigs with engineered repeats, written with its
.fai) and a synthetic panel parquet. No mocks.

**SNVs and palindromes:**

- SNV classes: keep, swap, complement, complement + swap, dropped.
- Palindromic SNV:
  - trusted: kept as-is;
  - untrusted: keep, strand-flip, unresolved-dropped, unresolved-kept when the
    option is set.

**Ambiguous indels:**

- Trusted: kept.
- Untrusted: every stringent branch (no EAF, not in panel, one record
  match/mismatch, both records decisive/indecisive).
- The mirrored-orientation symmetry: the same variant supplied in either
  orientation resolves identically.

**Trust decision:**

- One inconsistent SNV flips the decision.
- Below-minimum counts give untrusted.

**Flip registry:**

- A1FREQ columns and OR CI bounds flip correctly.
- An unregistered column fails.
- extra_column_rules is honoured.

**Execution:**

- Output is identical whether the input has one chromosome or rows
  interleaved across several.
- Excluded chromosome rows are dropped.
- Palindromic MNP: kept when trusted, dropped when untrusted.
- A pipe that switches the backend away from polars fails the backend
  assertion.
- A duplicate (CHR, POS, EA, NEA) after resolution fails.
- A very long allele (e.g. 5 kb) is classified correctly with a small gather
  byte budget.

**Reference assets:**

- ReferencePanelAlleleFrequencyTask on a tiny bgzipped VCF.
- IndexedFastaTask on a tiny gzipped FASTA.

## Rollout

1. **Reference assets and tasks** (FASTA download and index, panel download and
   parquet), with unit tests.
2. **Library module and GenomeReferenceHarmonizationTask**, with unit tests.
   Run pixi r invoke green (logged to file).
3. **Experiments V1-V5.** Set the indel defaults from V3 and write the results
   into this spec's directory as a short report.
4. **Switch call sites.**
   - Annovar generator: GWASLabTransformSumstatsTask(harmonize) plus
     GwasLabSumstatsToTableTask become
     GwasLabSumstatsToTableTask(pre-harmonization pickle) ->
     GenomeReferenceHarmonizationTask -> JoinDataFramesTask (unchanged).
     drop_palindromic_ambiguous maps to keep_unresolved_palindromes = not
     drop_palindromic_ambiguous.
     filter_indels_in_harmonized stays upstream in GWASLabCreateSumstatsTask
     semantics (applied before the dump).
   - Liu IBD harmonized and its dump are replaced the same way.
   - Rewrite test_harmonize_drop_ambiguous.py against the new Task.
5. **Remove dead code** once nothing uses it: HarmonizationOptions,
   _do_harmonization, the harmonize_options fields.
6. **Rebuild.** Changed task graphs rebuild downstream assets automatically:
   the asset IDs and dependencies change, and the cache does not key on code.
   Note in the PR which analyses will shift (indel sets, palindrome handling).

## Out of scope (tracked separately)

- **Liftover and basic_check** remain in gwaslab. Whole-table memory in
  GWASLabCreateSumstatsTask and in the post-rsID gwaslab Sumstats creation is
  unchanged by this work, so end-to-end OOM safety of the chain needs follow-up.
- **Mirrored-pair handling downstream:**
  - unordered_allele_key joins;
  - baseline-LF annotation dedup (the failing PolyFun SUSIE system test);
  - SUSIE prior and contrast join row-multiplication checks.
- **Task-reference harmonization fixes:** double-match refusal, flipping
  A1FREQ_CASES/CONTROLS.
- **The uncommitted unit-test edit** in
  test_build_baseline_lf_annotation_parquet_task.py belongs to the
  annotation-dedup fix, not this work.

## Decisions from review (2026-09-14)

1. **Trusted-mode palindromic SNVs** keep the source strand, with no frequency
   inference.
2. **Small tables** below the minimum trust counts (10,000 SNVs, 1,000 indels)
   take the stringent path. Accepted.
3. **MNPs** are oriented by the FASTA with no strand inference. Palindromic
   MNPs are kept at 100% trust and dropped otherwise.
4. **Large downloads** (panel VCF, FASTA gz) are discarded via
   DiscardDepsWrapper; only derived assets are stored.
5. **The pipe output** is asserted to be a polars LazyFrame after to_native().
