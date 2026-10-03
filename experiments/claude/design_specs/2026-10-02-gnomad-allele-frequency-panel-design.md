# gnomAD and Pan-UKBB Allele-Frequency Panels for Genome-Reference Harmonization --- Design

Date: 2026-10-02 (Pan-UKBB added 2026-10-03)

## Goal

Give genome-reference harmonization (GenomeReferenceHarmonizationTask) allele-frequency
panels built from large population resources, with every ancestry group's frequencies
available, and let each harmonization choose the panel and ancestry that match its GWAS.

- hg19, gnomAD: built now, from gnomAD v2.1.1 genomes (GRCh37).
- hg19, Pan-UKBB: built now, from the Pan-UK Biobank variant manifest (GRCh37, imputed UK
  Biobank). It complements gnomAD: far larger European sample, but imputed and limited to
  imputation-panel sites.
- hg38, gnomAD: same code, from gnomAD v4.1 genomes (GRCh38). The asset is defined in this
  project but built only on request.

The 1000 Genomes EUR panels stay. Existing harmonized outputs must not change.

Non-goals: changing any resolution rule (palindromes, ambiguous indels, trust); using AN in
the harmonizer; lifting v4 down to hg19 or Pan-UKBB up to hg38; adding an ancestry field to
GWAS metadata; combining panels into one table.

## Why

The panel resolves palindromic SNVs and ambiguous indels in untrusted tables, and drives the
suspicious-indel trust gate (palindromes.py, ambiguous_indels.py, trust.py). The 1000 Genomes
EUR panel has about 1,006 chromosomes. At AF 0.3 its sampling SE is about 0.014, which uses
up most of the ambiguous-indel tolerance (indel_max_af_distance = 0.02). gnomAD NFE has
15,436 chromosomes in v2.1.1 and 68,058 in v4.1, with SEs of about 0.004 and 0.002.
Pan-UKBB EUR has 917,874 (SE about 0.0005), and its population (mostly British) is the
closest match to DecodeME. Its frequencies are imputed dosages rather than sequenced calls,
so their accuracy rests on imputation quality, which the manifest bounds at INFO > 0.8. The
expected gain is mainly in fewer correct ambiguous indels dropped as AF_MISMATCH. That is
a prediction; validation step 3 measures it for each panel.

A panel restricted to SNVs would silently disable ambiguous-indel resolution and the trust
gate, so the panel keeps indels.

## Measured facts

All measurements are on chr21, from experiments/claude/gnomad_af_reference/ (scripts are
committed; logs are gitignored, so the numbers are recorded here).

Sources, from the public bucket gnomad-public-us-east-1, which allows anonymous access and
has no egress charge to us:

| | v2.1.1 genomes | v4.1 genomes |
|---|---|---|
| Path | release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.{chrom}.vcf.bgz | release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz |
| Build, contig naming | GRCh37, "21" | GRCh38, "chr21" |
| Chromosomes | 1-22, X | 1-22, X, Y |
| Total VCF size | 451 GiB (per-chromosome files) | 524 GiB |
| chr21 records, PASS | 3.48 M, 2.91 M | 10.96 M, 8.20 M |
| Multiallelics | already split (one ALT per record) | already split |
| FILTER values | PASS, RF, AC0, InbreedingCoeff | PASS, AS_VQSR, AC0, InbreedingCoeff |

Ancestry groups (approximate samples = max AN / 2, PASS records):

- v2.1.1: afr 4,359; amr 424; asj 145; eas 780; fin 1,738; nfe 7,718; oth 544; and the NFE
  subgroups nfe_nwe 4,299, nfe_est 2,297, nfe_onf 1,069, nfe_seu 53. There is no sas group.
- v4.1: afr 20,805; amr 7,657; asj 1,736; eas 2,604; fin 5,316; nfe 34,029; sas 2,419;
  ami 456; mid 158; remaining 1,058.

Within every group, median AN is within about 0.3% of max AN. In both releases AF is missing
exactly where AN = 0. None of the v4.1 genomes are UK Biobank samples, so v4's non_ukb strata
(exome-only) are irrelevant here.

Streaming chr21 over HTTPS with bcftools took 11.5 min (v2) and 11.8 min (v4), about
10 MiB/s, with bcftools using 40-65% of one core. A serial genome-wide pass is therefore about
14 h per release.

Panel size at the chosen filter (PASS, and AF > 0 in at least one main group), storing AF and
AN for every group:

| | chr21 rows | genome-wide rows (extrapolated) | genome-wide size, chosen encoding |
|---|---|---|---|
| v2.1.1 | 2.84 M | ~210 M | ~3.7 GiB |
| v4.1 | 8.10 M | ~590 M | ~10 GiB (int32 byte-stream-split for AN was measured on v2 only; v4 is extrapolated from its 10.9 GiB at default encoding) |

Encoding measurements (v2 chr21 at the chosen filter; v4 agrees):

- Byte-stream-split on the AF columns makes them 50-75% larger and slows reads. AFs are
  quantized (AC / AN, about 6 significant digits) and mostly exactly 0, so dictionary + RLE
  wins. Byte-stream-split requires disabling the dictionary.
- AN columns at default encoding cost about as much as all the AF columns together. AN takes
  almost every even value from 0 to max (7,712 distinct values for v2 nfe). Its order-0
  entropy sums to 42.5 bits per row over 11 groups, a floor of 14.4 MiB for chr21.
  Int32 byte-stream-split reaches 16.0 MiB at zstd default and 14.3 MiB at zstd 19, against
  19.9 MiB for dictionary. Float32 byte-stream-split reaches 18.8 MiB and delta encoding
  28.1 MiB.
- zstd level: default (3) writes chr21 in 3.1 s at 55.0 MiB, level 9 in 6.0 s at 52.3 MiB,
  and level 19 in 71 s at 46.1 MiB. pyarrow compresses on one core.

REF vs FASTA: all 2,844,308 v2 chr21 REF alleles match the UCSC hg19 FASTA exactly. No
record lies over an N or IUPAC base, and no REF or ALT has a non-ACGT character.

### Pan-UKBB variant manifest

Measured genome-wide by experiments/claude/pan_ukbb_manifest/ (inspect_manifest.py,
inspect_swaps.py).

- Source: https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz,
  public S3. One bgzipped TSV of 2.7 GB (9.0 GB uncompressed) with a tabix index, last
  modified 2020-08-28, md5 e70ebc8289f762dd8d5086f54e766654 (equal to its S3 ETag).
- Content: UK Biobank imputed variants (HRC + UK10K/1000 Genomes imputation, GRCh37) with
  INFO > 0.8. 28,987,534 rows on 1-22 and X, in gwaslab chromosome order; 26.99 M SNVs and
  2.00 M indels; no duplicate keys; all alleles ACGT; one ALT per row.
- Columns used: chrom ("1".."22", "X"), pos, ref, alt, and af_{pop} for AFR, AMR, CSA, EAS,
  EUR, MID. af is the frequency of the listed alt. It has no nulls, and af = 0 occurs (EUR
  3,521 rows, EAS 1.16 M).
- ac is a fractional dosage sum and an is constant within a contig (EUR 917,874 on the
  autosomes, 916,634 on X), so AN carries no per-site information. Pan-UKBB reports 420,531
  EUR individuals against the 458,937 implied by an; the difference is unexplained and does
  not matter for orientation.
- Other columns, not used: rsid, info, gnomAD v2 genome frequencies joined in, and
  high_quality (passes gnomAD and agrees with gnomAD frequency in all four shared groups).
- REF vs the UCSC hg19 FASTA: 28,987,239 rows match. 295 do not (chr21 47, chr22 86, X 162).
  They are directly genotyped sites (info = 1.0), and every one is a ref/alt swap: SNVs
  with the FASTA base as alt, and insertions written as deletions (manifest CT>C where
  gnomAD has C>CT). All have high_quality = false and no gnomAD match, because Pan-UKBB's
  gnomAD join is on ordered alleles. No row lies over an N or IUPAC base.
- The af of a swapped row describes its listed alt, which is the reference allele. On chr21
  (swaps_vs_gnomad_chr21.py), 45 of the 47 have a gnomAD v2 record with the same alleles in
  FASTA orientation, and af_EUR is within 0.05 of 1 - gnomAD AF_nfe for 44 of them.
- Origin: the Pan-UKBB pipeline (github.com/atgu/ukbb_pan_ancestry, resources/genotypes.py)
  takes alleles straight from UK Biobank's v3 imputed BGEN files with hl.import_bgen and
  never checks them against a reference. The swap is therefore in UK Biobank's own allele
  order for these sites. Why only chr21, chr22 and X (with 89 of the 162 X rows at 88-92
  Mb) is not determinable from public data.

## Decisions

1. **Builds.** v2.1.1 genomes for hg19 now, v4.1 genomes for hg38 with the same code.
2. **One streaming pass per release**, storing AF_g and AN_g for every ancestry group of the
   release.
3. **Filter**: FILTER == PASS and AF_g > 0 for at least one main group. Main groups are the
   large continental groups: v2 afr, amr, asj, eas, fin, nfe; v4 those plus sas. Extra groups
   (NFE subgroups, oth, ami, mid, remaining) are stored but do not admit records, since a
   single allele in a small group would admit noise. Any nonzero AF threshold would mostly
   admit records through small groups: at 1e-3 on v2, 20% of rows enter only through asj or
   amr. Size does not justify a threshold.
4. **Encoding.** AF_g Float32 with dictionary encoding; AN_g Int32 with byte-stream-split;
   zstd default level. Keys CHR (Int32 gwaslab code), POS (Int32), REF, ALT (String).
5. **Ancestry is chosen at harmonization time**, through a required argument.
6. **AN is stored but not read by the harmonizer** in this project. It preserves the option of
   sample-size-aware rules later.
7. **gnomAD REF must match the FASTA**, enforced while building.
8. **Pan-UKBB is a separate hg19 panel**, not merged with gnomAD. Each harmonization picks
   one panel and one ancestry; validation step 3 compares them.
9. **Pan-UKBB keeps every manifest row except REF/FASTA mismatches.** No frequency filter
   (the manifest is already filtered) and no high_quality filter. That filter would drop
   2.3 M further rows (1.39 M absent from gnomAD, 0.73 M gnomAD-PASS but frequency-
   discordant, 0.18 M failing gnomAD filters) and would couple the panel to gnomAD v2.
10. **The 295 Pan-UKBB mismatches are dropped, and their count is asserted exactly.** The
    download is pinned by md5, so the count is deterministic. Asserting it equal to 295,
    rather than tolerating some, keeps the zero-surprise guarantee: any change in the file or
    in our parsing fails the build. They are dropped rather than swapped back (with af
    replaced by 1 - af): gnomAD supports swapping back for 44 of 47 chr21 rows but not all,
    confirming the rest would couple this panel to gnomAD, and 295 rows are 0.001% of the
    panel.
11. **Pan-UKBB stores AF only.** Its AN is constant per contig and would suggest per-site
    information that does not exist.

## Architecture

Classes live under mecfs_bio/build_system. Concrete releases, URLs and asset instances live
under mecfs_bio/assets. Vocabularies that build_system code types against live in
mecfs_bio/constants.

### Vocabulary (mecfs_bio/constants/)

- GnomadGroup = Literal["afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "nfe_est",
  "nfe_nwe", "nfe_onf", "nfe_seu", "oth", "remaining", "sas"]
- PanUkbbGroup = Literal["ukb_afr", "ukb_amr", "ukb_csa", "ukb_eas", "ukb_eur", "ukb_mid"].
  The prefix keeps Pan-UKBB's genetically assigned groups distinct from gnomAD's groups and
  from 1000 Genomes "eur".
- PanelAncestry = GnomadGroup | PanUkbbGroup | Literal["eur"]. "eur" is the 1000 Genomes EUR
  super-population and is deliberately not a synonym for gnomAD nfe or ukb_eur.

### Panel column declaration (mecfs_bio/build_system/meta/)

- PanelAlleleFrequencyColumns: frozen. Holds (ancestry, column) pairs as a tuple of tuples, so
  it stays hashable on a frozen Task. __attrs_post_init__ asserts the ancestries are unique.
  column_for(ancestry) returns the column, and the assertion message lists the available
  ancestries.
- HarmonizableReferenceTableMeta gains allele_frequency_columns: PanelAlleleFrequencyColumns |
  None = None. The default keeps the other users (UKBB LD labels, RenameColsTask) unchanged.

### gnomAD extraction (mecfs_bio/build_system/task/genome_reference_harmonization/gnomad/)

**GnomadRelease** (gnomad_release.py): a frozen descriptor holding

- name;
- build;
- vcf_url_template, with a chromosome placeholder;
- contig naming (the contig name for a gwaslab code);
- chromosomes, as gwaslab codes;
- main_groups and extra_groups (GnomadGroup tuples);
- header_assembly, the header's contig assembly string.

__attrs_post_init__ asserts that the groups are disjoint and non-empty, that the chromosome
list is non-empty and has no MT, and that the template contains the placeholder.

**GnomadChromosomeAlleleFrequencyTask** (gnomad_chromosome_allele_frequency_task.py)

- Fields: meta, release, chrom (gwaslab code), fasta_task, read_block_size_bytes (production
  default), and the retry policy (max attempts and sleep; production defaults from
  execute_command_with_retries).
- Deps: fasta_task. create() asserts the FASTAMeta build equals release.build, and that chrom
  is in release.chromosomes.

execute():

1. **Header check.** Run bcftools view -h on the URL. Assert every AF_g and AN_g for the
   release's groups is declared, and that the header's contig assembly equals
   release.header_assembly. This fails in seconds, before a long stream.
2. **Stream, as one process**: bcftools query -i 'FILTER="PASS" && (AF_g1>0 || AF_g2>0 ...)'
   -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF_g\t%INFO/AN_g...\n' URL -o scratch TSV. Run through
   execute_command_with_retries; each attempt overwrites the TSV. It must be one process, not
   a view | query pipe. execute_command runs through sh, where a pipeline's exit status is the
   last command's, so a network failure upstream could yield a truncated TSV and exit 0. A
   single process lets htslib's BGZF CRC and EOF checks fail the Task.
3. **TSV to parquet, memory-bounded.** pyarrow's streaming CSV reader (explicit schema, "." as
   null, block size read_block_size_bytes) feeds a ParquetWriter. Each batch passes through
   the shared panel batch checks (below), with the release's contig naming as the only
   allowed contig. Then, for every group, AF must be null exactly where AN == 0.
   **REF vs FASTA**: any mismatch over a pure-ACGT span fails the Task (the batch check's
   mismatch count must be zero).
4. The output must have at least one row. Log row counts, FASTA-ambiguous drops and per-group
   max AN with structlog.

### Shared panel batch checks (mecfs_bio/build_system/task/genome_reference_harmonization/panel_batch_checks.py)

Both the gnomAD and the Pan-UKBB extraction stream large sorted tables into parquet and need
the same invariants, so one helper holds them. PanelBatchChecker is a small stateful object
(it carries the previous row across batches) constructed with the IndexedFasta and the map
from source contig name to gwaslab code. Its check(batch) returns the cleaned batch and
per-batch counts as a frozen result object:

- The contig column must hold only mapped names; it is replaced by an Int32 CHR gwaslab code.
- REF and ALT must be non-null and ACGT-only.
- (CHR, POS) must be non-decreasing, and there must be no duplicate (CHR, POS, REF, ALT).
  Duplicates are adjacent in sorted input, so each row is compared with its predecessor, and
  the last row is carried across batch boundaries.
- **REF vs FASTA.** Records whose FASTA span contains N or another IUPAC code are dropped and
  counted as fasta_ambiguous: their orientation is unverifiable, and the harmonizer's own
  FASTA classification would reject them anyway. Records whose REF differs from a pure-ACGT
  span are dropped and counted as ref_mismatch. The calling Task decides what count is
  acceptable: gnomAD requires zero, Pan-UKBB requires exactly its expected count.

At the end the caller asserts that every expected contig appeared.

Writer settings: the AN columns are passed as byte_stream_split columns, AF columns keep the
dictionary, and compression is zstd at the default level. These settings come from one helper
in dataframe_output.py that returns the pyarrow writer keyword arguments. write_parquet_table
and both streaming writers use that helper, so the dictionary-versus-split rule lives in one
place. Its docstring is updated: the split applies to any fixed-width column (pyarrow 25
supports Int32), not only floats.

FASTA gather: add a sibling of reference_matches in fasta.py that reports, per record, whether
the reference span is pure ACGT. It reuses the same memory-bounded, case-insensitive gather.
The Task combines it with reference_matches.

**GnomadAlleleFrequencyPanelTask** (gnomad_allele_frequency_panel_task.py)

- Deps: one GnomadChromosomeAlleleFrequencyTask per release chromosome.
- execute(): reads the parts in gwaslab code order (1-22, X=23, Y=24) with
  ParquetFile.iter_batches. Asserts every part has the identical arrow schema. Writes one
  FileAsset with the same writer settings. The output is sorted by (CHR, POS), so
  ParquetPanelLoader's per-chromosome filter prunes row groups.
- create() builds HarmonizableReferenceTableMeta with:
  - harmonization_info: build = release.build, ref_allele_col REF, pos_col POS;
  - allele_frequency_columns = {g: "AF_g"} for every group of the release.
- **Do not wrap this Task in DiscardDepsWrapper.** That wrapper rebuilds the inner Task's
  whole dependency graph inside a temporary store. The 14 h build would become all-or-nothing
  again, and the FASTA would be rebuilt in the temporary store. The per-chromosome parts
  therefore stay in the asset store, costing roughly one extra panel's worth of disk; they
  can be deleted by hand once the panel exists.

### Pan-UKBB extraction (mecfs_bio/build_system/task/genome_reference_harmonization/pan_ukbb/)

**PanUkbbAlleleFrequencyPanelTask** (pan_ukbb_allele_frequency_panel_task.py)

- Fields: meta, manifest_task (the downloaded manifest), fasta_task, groups (the
  PanUkbbGroup values to extract, mapped to manifest population codes), chromosomes (gwaslab
  codes expected), expected_ref_mismatches (int), read_block_size_bytes (production
  default).
- Deps: manifest_task and fasta_task. create() asserts the FASTAMeta build is GRCh37 (the
  manifest is GRCh37-only, so this is a constant of the class, not a field) and builds
  HarmonizableReferenceTableMeta with harmonization_info build GRCh37, ref_allele_col REF,
  pos_col POS, and allele_frequency_columns {ukb_g: "AF_ukb_g"} for every group.

execute():

1. **Header check.** Read the first line of the manifest and assert that chrom, pos, ref, alt
   and af_{pop} for every configured group are present.
2. **Stream the bgzipped manifest directly** with pyarrow's streaming CSV reader over a gzip
   input stream (tab separator, "NA" as null, explicit schema for the selected columns,
   other columns skipped). BGZF is multi-member gzip, and a reader that stopped after the
   first member would silently truncate the panel. The plan's first step confirms pyarrow
   consumes every member, and a unit test with a fixture spanning several BGZF blocks guards
   it. Fallback if not: decompress to scratch with bgzip -dc (9 GB) first.
3. Per batch: rename af_{pop} to AF_ukb_{pop} as Float32, assert no null AF, then run the
   shared panel batch checks with contigs "1".."22", "X" mapped to gwaslab codes.
4. At the end: assert the total ref_mismatch count equals expected_ref_mismatches, that no
   row was FASTA-ambiguous (the measurement found none; a change is a surprise), that every
   chromosome appeared, and that the output is non-empty. Log counts with structlog.

The output is one parquet sorted by (CHR, POS), because the manifest is already in gwaslab
order and the batch checks enforce it. Writer settings are the shared helper's, with no AN
columns. Estimated size is well under 1 GiB; validation step 2 records it. The whole build
reads 2.7 GB locally and should take minutes, so no per-chromosome split is needed.

### Assets (mecfs_bio/assets/reference_data/gnomad/ and mecfs_bio/assets/reference_data/pan_ukbb/)

- GNOMAD_V2_1_1_GENOMES and GNOMAD_V4_1_GENOMES: GnomadRelease instances with the values in
  Measured facts.
- GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES (with UCSC_HG19_INDEXED_FASTA) and
  GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES (with the hg38 indexed FASTA).
- PAN_UKBB_VARIANT_MANIFEST: a DownloadFileTask with the manifest URL and its md5.
- PAN_UKBB_HG19_ALLELE_FREQUENCIES: PanUkbbAlleleFrequencyPanelTask with all six groups,
  chromosomes 1-22 and X, expected_ref_mismatches = 295, and UCSC_HG19_INDEXED_FASTA. Like
  the gnomAD panel it is not wrapped in DiscardDepsWrapper, which would rebuild the FASTA in a
  temporary store; the 2.7 GB download stays in the asset store and can be deleted by hand.

### Harmonizer changes

- ReferencePanelAlleleFrequencyTask.create gains a required ancestry: PanelAncestry, and
  declares {ancestry: "AF"}. The 1000 Genomes assets pass "eur". Their parquet files do not
  change.
- GenomeReferenceHarmonizationTask.create gains a required panel_ancestry: PanelAncestry. It
  is required rather than defaulted, because a silent default is how the wrong ancestry slips
  in. create() asserts the panel meta has allele_frequency_columns that include
  panel_ancestry, then stores the resolved column as the field panel_af_col.
- ParquetPanelLoader gains af_col. It selects pl.col(af_col).alias(PANEL_AF_COL) and drops
  rows where it is null (AN == 0: nobody in that group was called, so the record is treated
  as absent). Both construction sites (trust pass, per-chromosome pass) receive the column.
- palindromes.py, ambiguous_indels.py, trust.py and GenomeReferenceHarmonizationOptions do not
  change.
- Call sites: the annovar_37_basic_rsid_assignment generator, the IBD liu_et_al_2023 asset and
  the test_harmonize_drop_ambiguous system test pass panel_ancestry="eur". The unit-test
  fixture's FakeTask panel meta declares {"eur": "AF"}.
- ReferencePanelAlleleFrequencyTask gets a comment noting that its bcftools view | query pipe
  has an unreliable exit status under sh, and that a single bcftools query -i would fix it.
  The fix itself is out of scope.

## Data flow

gnomAD VCF URL (per chromosome)
-> bcftools query -i PASS-and-polymorphic (one process, retried)
-> scratch TSV
-> pyarrow streaming reader, with invariants and REF/FASTA checks per batch
-> chromosome parquet (per-chromosome Task, cached, so a failure costs one chromosome)
-> GnomadAlleleFrequencyPanelTask concatenation
-> panel parquet with CHR, POS, REF, ALT, AF_g..., AN_g...
-> ParquetPanelLoader selects AF_<panel_ancestry> as AF and drops null AF
-> unchanged resolution rules.

Pan-UKBB manifest URL (md5-pinned DownloadFileTask)
-> local bgzipped TSV
-> pyarrow streaming reader over gzip, shared batch checks, 295 REF mismatches dropped
-> panel parquet with CHR, POS, REF, ALT, AF_ukb_g...
-> ParquetPanelLoader selects AF_ukb_<group> (never null here)
-> unchanged resolution rules.

## Error handling

Fail the Task on any of the following:

- missing INFO field, or wrong header assembly;
- bcftools non-zero exit after retries, which includes BGZF CRC/EOF failures on truncated
  or corrupt streams;
- unexpected contig;
- null or non-ACGT allele;
- decreasing POS or duplicate key;
- AF/AN null inconsistency;
- REF mismatch over pure-ACGT reference (gnomAD), or a Pan-UKBB mismatch count other than
  expected_ref_mismatches;
- Pan-UKBB: missing manifest column, null AF, any FASTA-ambiguous row, or a missing
  chromosome;
- md5 mismatch on the manifest download (DownloadFileTask);
- empty output;
- schema mismatch between parts;
- FASTA/release build mismatch, or an undeclared panel ancestry (both fail at graph
  construction).

Records over N/IUPAC reference bases are dropped and logged, not fatal. If v4 turns out to
contain non-ACGT alleles such as spanning-deletion "*", the assertion reports it and the
handling is decided then.

## Testing

Tests are Task-level, use no mocks, and do not match on error text.

Fixtures: a tiny bgzipped VCF written in tmp_path with the release's INFO fields, a test
GnomadRelease whose URL template points to that local file (bcftools reads local paths and
URLs the same way), and the synthetic FASTA from genome_reference_fixtures.

GnomadChromosomeAlleleFrequencyTask:

- **Happy path.**
  - Kept: PASS records polymorphic in a main group.
  - Dropped: non-PASS records, and records whose only nonzero AF is in an extra group.
  - Preserved: AF null where AN = 0.
  - Output: expected columns, dtypes and CHR code; parquet metadata shows BYTE_STREAM_SPLIT
    on the AN columns and dictionary encoding on AF.
- REF mismatch over ACGT fails.
- REF over an N or IUPAC FASTA base is dropped, not fatal.
- A header missing an expected field fails.
- A duplicate key fails, including a pair split across a batch boundary (forced with a tiny
  read_block_size_bytes). A decreasing POS fails.
- A truncated VCF fails (one retry attempt, no sleep).

GnomadAlleleFrequencyPanelTask:

- Parts are concatenated in gwaslab chromosome order (X after 22), and the encodings are
  preserved.
- Mismatched part schemas fail.

PanUkbbAlleleFrequencyPanelTask (fixture: a synthetic manifest with the real header, bgzipped
in tmp_path, large enough to span several BGZF blocks, and the synthetic FASTA):

- **Happy path.** Every row arrives (guards multi-member gzip reading); columns AF_ukb_g are
  Float32; CHR codes and (CHR, POS) order are right; the meta declares every group.
- Rows whose ref is swapped relative to the FASTA are dropped when their count equals
  expected_ref_mismatches, and the Task fails when it does not.
- A row over an N FASTA base fails.
- A missing af column fails, and out-of-order chromosomes fail.

The shared batch checks are covered through these two Tasks' tests, not tested directly.

GenomeReferenceHarmonizationTask:

- Existing tests pass with panel_ancestry="eur" and unchanged expectations.
- A panel with AF_nfe and AF_nfe_nwe set so that one palindrome resolves KEEP under one
  ancestry and STRAND_FLIP under the other, proving the chosen column is read.
- A null chosen AF makes an ambiguous indel NOT_IN_PANEL.
- create() rejects an undeclared ancestry, and a panel without allele_frequency_columns.

Not tested: GnomadRelease and PanelAlleleFrequencyColumns as data holders, whose invariants
run at import when the asset constants are built; asset definitions, which are wiring only;
a CI test against live gnomAD, which would need a test-only region field and network access,
while the real build enforces every invariant genome-wide.

pixi r invoke green must pass, including import-linter.

## Validation on real data (experiments/claude/gnomad_af_reference/)

1. **Byte identity.** Force-rebuild one existing 1000 Genomes harmonized asset (DecodeME or
   IBD) after the refactor and confirm the parquet output is identical. The build cache does
   not track code changes, so the rebuild must be forced.
2. **Full builds.** gnomAD v2.1.1 (about 14 h) doubles as the genome-wide REF/FASTA and
   invariant check; record its size and FASTA-ambiguous drop count. Pan-UKBB (minutes): record
   its size and confirm the 295-mismatch assertion holds through the Task's own parsing.
3. **DecodeME comparison.** Re-harmonize with four panel/ancestry choices: 1000 Genomes "eur",
   gnomAD v2.1.1 "nfe" and "nfe_nwe", and Pan-UKBB "ukb_eur". For each, report:
   - drop-reason counts, especially AF_MISMATCH for ambiguous indels and NOT_IN_PANEL
     (panel coverage: Pan-UKBB lacks sites outside the imputation panels);
   - palindromes and ambiguous indels resolved, and how often the four panels agree on each
     resolution (a disagreement between gnomAD and Pan-UKBB is a variant worth inspecting);
   - palindromes decided by a chosen-group AF of exactly 0.
   The script lives in experiments/claude/gnomad_af_reference/ and tees its log.

## Deferred

- **Palindromes decided by AF = 0.** At the chosen filter, a gnomAD panel has many rows where
  the chosen group's AF is exactly 0, for example a variant seen in afr but never in nfe. The
  palindrome rules treat panel AF 0 as decisive and call KEEP when EAF < 0.5. The ambiguous-indel
  rules already treat monomorphic records as absent. 1000 Genomes EUR has the same behaviour
  today, but gnomAD will make it more frequent. Decide from validation step 3 whether
  palindromes should also treat monomorphic records as absent.
- **Sample-size-aware rules using AN.**
- **The pipe fix in ReferencePanelAlleleFrequencyTask** (a comment only, in this project).
- **Building the hg38 panel** (about 14 h, about 10 GiB).
- **Choosing a default panel** for new harmonizations, once step 3 has been reviewed.
- **Pan-UKBB on hg38.** No GRCh38 manifest is known; liftover is out of scope.
- **Recovering the 295 swapped Pan-UKBB rows** by swapping alleles back and using 1 - af.
  The chr21 evidence supports it; it is deferred only because the rows are too few to
  matter.
