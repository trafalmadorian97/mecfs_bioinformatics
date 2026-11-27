# deCODE SomaScan database: probes and findings

What the exploratory probes behind the deCODE SomaScan database established (Nov 2026). The
probe scripts were removed once their results were recorded here. Access, bucket layout, file
format, download throughput and the first end-to-end mockup are in the design spec's "What the
survey established" (`../design_specs/2026-11-21-decode-somascan-database-design.md`). Every
measurement below used one SMP file (AP4M1, 10076-1) unless it says otherwise.

## Pipeline design probes

- **Decompression with a seek index.** rapidgzip builds an indexed_gzip-compatible index
  during the one sequential pass. One 3.32 GB file took 4.3 s at 8 threads (16 threads 3.9 s,
  1 thread 13.5 s). indexed_gzip took 26–29 s, zlib 13–19 s, and ISA-L 4.1 s but with no index.
  Every exported index was checked at 20 random offsets.
- **Where a pass spends its time.** Parsing all 33.4M lines cost 2.8 s, the newline scan 1.6 s,
  and a HapMap3 semi-join on all rows 5–7 s. Parsing only the needed lines (cut by the newline
  scan) gave the same result as the full parse plus join. A position-first join and an integer
  key gave no gain over the string-pair join.
- **Significant rows.** AP4M1 has 849 rows at p ≤ 1e-5, 377 at 1e-6 and 192 at 5e-8 (none on
  chrX). Filtering must use minus_log10_pval, because Pval rounds to 0 for strong hits. This
  sized the p ≤ 5e-8 supplement at about 10 KB per protein.
- **HapMap3-only screen.** In 25 proteins, 320 of 471 loci at p < 5e-8 have no significant
  HapMap3 variant. Most are weak single-variant loci; about 60 look real. This is why the
  per-protein task keeps significant_rows.parquet over all variants.
- **Downloads.** aria2 caps itself at 16 connections per server. With 128 connections on one
  972 MB file, s5cmd took 50 s, rclone 45 s, and plain Python threads with range requests 40 s
  (24 MB/s, near the server's ceiling), so the database uses its own parallel range downloader.
  Chunks of 4 MiB were as fast as 2 MiB, and 8 MiB was slower.

## Row order across files (the RAF1 finding)

- The gzip trailer's uncompressed size differs for every file (the numbers' text lengths vary),
  so it cannot show which files have different rows.
- RAF1 (10001-7) is the template minus 6 rows, none of them a HapMap3 index row; every other row
  is identical and in the same order. Gathering by template row number therefore failed. Rows
  are now matched to the index by Name.
- Name equals `Chrom:Pos:effectAllele:otherAllele` in all 33,449,999 AP4M1 rows. Index NAMEs are
  unique, but the file holds 133,250 duplicated Names, so the join is validated one-to-one.
- **Join cost.** A position pre-filter keeps about 1.4M candidate lines per file, and only those
  are parsed. The Name join then costs 0.41 s, the same as joining on (CHR, POS, EA, NEA).
  Parsing Name/Beta/SE/N for every row instead would add about 31 s per protein.
- **24-protein check.** The first 24 SMP proteins all built with 0 of 1,171,498 index rows
  missing.

## Allele semantics and multi-allelic sites

- **The two alleles.** otherAllele is the GRCh38 reference in 100% of rows. A row tests
  effectAllele against every other allele at the site.
  - Rows whose two alleles are equal (X:X, 3.0M rows) test the reference allele itself.
  - `*` (958k rows) is the spanning-deletion allele.
  - Some per-base rows copy pieces of a nearby complex allele.
- **Multi-row positions.** 27.0% of rows (9.03M) sit at 2.74M positions holding more than one
  row: 6.4% of common SNVs and 72% of common indels.
- **variant_info_iceland.**
  - A two-allele record lists only its ALT allele.
  - A multi-allelic record lists every allele, REF and `*` included, and its frequencies sum
    to 1.
  - Records overlapping one position are interleaved, so they cannot be separated.
- **Why multi-allelic rows are a problem.** At a multi-allelic site, a deCODE row is one allele
  against the pooled rest. A trait GWAS row is ALT against REF; after `bcftools norm -m-`
  (checked on a demo VCF), the other ALT alleles' carriers are counted as REF. So the two rows
  measure different contrasts.
- **The strict cleanliness rule.** A row is clean when variant_info lists only its effect allele
  at its position. The rule keeps:
  - 93.2% of common SNVs (dropping multi-row positions alone would keep 93.6%);
  - 99.5% of SNVs alone at their position;
  - 23.9% of common indels.

  It also covers the X:X, `*` and multi-row cases. On real data:
  - the clean table holds 23.9M of 27.2M positions;
  - 1,171,498 of the 1,197,450 HapMap3 index rows are clean.
- **An abandoned measure.** A "pooled other-allele mass" (q) was explored and dropped, because
  variant_info's interleaved records make it unreliable.
- **EAF.** In the files the effect allele is always the ALT; ImpMAF does not say which allele
  is minor. EAF is the variant_info frequency of the effect allele, checked against ImpMAF
  within 1e-5:
  - 93.0% of rows have a single variant_info candidate, and it always agrees;
  - 7.0% have several candidates, and in 6.7% only one agrees;
  - no row has candidates of which none agrees.

## DecodeME (trait side)

- **What the papers say** (`decodeme_preprint.pdf`, `decode_me_supplement.pdf` pp. 11–12):
  - association used REGENIE with Firth correction;
  - imputation used IMPUTE5 with the UK Biobank 200k WGS reference panel;
  - DecodeME "discarded … any variant that was monomorphic, multi-allelic, or very rare".
- **The raw GWAS-1 file** has 8,941,712 records at 8,941,712 unique positions. The harmonized
  GRCh38 table has 8,902,782 rows, all with MAF ≥ 1%. So DecodeME has no multi-record positions,
  and dropping them for deCODE coloc runs costs nothing here. Other sources keep such positions.
- **Effect scale.** deCODE's Beta is in protein-SD units: SE·√(2Nf(1−f)) has a median of 1.05
  (AP4M1) and 1.02 (CRYBB2), consistent with coloc's sd_y = 1.
- **EAF concordance with deCODE.** 8,284,093 clean AP4M1 rows were joined exactly to DecodeME.
  - The EAF correlation is 0.9967.
  - The 99.99th percentile of |ΔEAF| is ≤ 0.14 for common SNVs, and the distribution's body
    ends near 0.15.
  - 224 variants exceed 0.2 (87 of them indels), and few of those fit a flipped allele.

  This sets coloc's max_eaf_difference = 0.2.

## Smoke run and tuning (plan Task 17, real bucket, 16 GB machine)

| stage | wall | peak RSS |
|---|---|---|
| variant_info parquet | 60 s | 11.0 GB |
| clean variant_info | 8 s | 6.9 GB |
| HapMap3 index (template CRYBB2) | 64.5 s | 4.4 GB |
| protein (download ~35–40 s + one pass ~10 s) | ~51 s | 5.4 GB |

- **Correctness.** The AP4M1 slim equals a brute-force Name join. Three remote windows equal an
  awk scan restricted to clean rows, with no null EAF. They hold 4, 4,678 and 8,829 rows (chr1,
  chr11, chrX).
- **N per variant.** p05 is 35,662, the median 35,690 and p95 35,709, so the median-N
  heritability is sound.
- **Proteins fetched at once.** One 1 Mb window from each of 24 proteins took about 20 s at 4
  proteins at once, 14 s at 6, 12 s at 8, 9–12 s at 12 and 11–13 s at 16. Peak RSS was 1.7 GB.
  DECODE_MAX_PARALLEL_FETCHES = 8.
