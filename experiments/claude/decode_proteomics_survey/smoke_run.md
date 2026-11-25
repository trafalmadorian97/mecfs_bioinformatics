# deCODE SomaScan smoke run (plan Task 17), 2026-11-24

Real-bucket runs on the 16 GB WSL machine. Scripts: smoke_run.py (one stage per process under
/usr/bin/time -v), verify_smoke.py (brute-force checks against a local AP4M1 copy),
tune_defaults.py, diff_row_sets.py / shell diff. Logs are not committed.

## Stages (step 1)

| stage | wall | peak RSS | notes |
|---|---|---|---|
| variant_info | 60 s | 11.0 GB | 442 MB download; scan/unique/sort/sink in polars |
| SMP HapMap3 index (template CRYBB2 10000-28) | 75 s | 4.4 GB | 1,197,450 rows (as expected) |
| protein CRYBB2 10000-28 | 57 s | 5.3 GB | download 39.9 s, one pass 8.4 s, write 0.2 s |
| protein AP4M1 10076-1 | 53 s | 5.2 GB | download 35.4 s, one pass 9.0 s, write 0.2 s |

- Per protein: ~50-55 s, as expected (~54 s). Peak RSS 5.2-5.3 GB, slightly above the 5 GB
  target, with the 16 GB machine far from OOM.
- variant_info peaks at 11 GB: a one-off task, but the closest the build comes to the limit.
- AP4M1 slim == brute force (gzip -dc | cut, rows taken at the index's SOURCE_ROW;
  assert_frame_equal on float32 BETA, SE, N; 1,197,450 rows).
- significant_rows.parquet (AP4M1): 192 rows, 10.7 KB. The screen costs 3.1 s of a 7.4 s
  decompress-and-block pass (expected under 1 s): about 6% of a protein's 50 s.

## Window reads (step 3)

Three AP4M1 windows through DecodeWindowReader (remote, 7.2 s for all three, including the
variant_info scans) equal an awk scan of the local file, EAF aside:

| window | rows | null EAF |
|---|---|---|
| chr1:1-500,000 | 4 | 0 |
| chr11:116.5-117.0 Mb | 6,955 | 5 |
| chrX:154.0-156.1 Mb | 12,202 | 4 |

## N spread (step 4)

AP4M1 slim N: p05 35,662, median 35,690, p95 35,709 (min 34,839, max 35,781). Narrow, inside
the expected 30.4k-35.7k range: the median-N heritability is sound.

## Download chunk size (step 2)

AP4M1, 128 connections, two runs each: 2 MiB 37.8 / 38.9 s; 4 MiB 37.3 / 38.9 s; 8 MiB
43.0 / 44.1 s. 4 MiB (the default) stays.

## Blocker: files do not all share the template's rows

Building the first 24 SMP proteins stopped at RAF1 (10001-7): 33,449,993 rows against the
template's 33,449,999. A diff of RAF1 against AP4M1 (which matches the template) shows RAF1 is
the template minus 6 rows (chr16:81416170:G:T, chr16:82140032:T:C, chr18:5778560:T:C,
chr20:24684085:G:C, chr6:98465792:G:T, chr6:134629413:G:A), all other rows identical and in
order. None of the six is a HapMap3 index row, but every index row after the first deletion
(760,615 of them) sits at a different row number, so gathering by SOURCE_ROW cannot build
RAF1. The task's Name check caught it as designed.

The gzip trailer's uncompressed size differs for every file (the numbers' text lengths vary),
so it cannot count rows; how many files have deletions is unknown until each is decompressed.

## Resolution

- Rows are matched to the index by Name (commit 55409f3f). Each block parses Pos together with
  the significance screen's field, keeps the lines at an index position, parses only those, and
  left-joins the index NAME on Name (validated one-to-one). Index rows a file lacks are null,
  counted in the slim's metadata (index_rows_missing). AP4M1 one pass 10.3 s, 2.2 GB peak
  (check_one_pass_name_join.py); slim == brute-force join. RAF1 builds: 0 index rows missing,
  51 s end to end, 5.4 GB peak.
- Only clean rows (two-allele tests) enter the index and coloc windows (commit 22659a28). The
  clean variant_info table keeps entries at positions with one entry: 23.9M of 27.2M positions
  (8 s, 6.9 GB peak). HapMap3 index: 1,171,498 of 1,197,450 rows clean, every EAF determined
  (64.5 s, 4.4 GB). verify_smoke.py: slim == brute force; the three windows equal the awk scan
  restricted to clean rows, with no null EAF: chr1 4, chr11 4,678, chrX 8,829 rows (were 4,
  6,955, 12,202 with 0, 5, 4 null EAF).
- Coloc drops aligned variants whose EAFs differ by more than 0.2 (commit 7100c8d6). On
  8,284,093 clean AP4M1 rows joined exactly to DecodeME the correlation of EAF is 0.9967; the
  99.99th percentile of |difference| is 0.14 for SNVs with MAF > 20%; 224 exceed 0.2 (103
  SNVs, 34 palindromic SNVs, 87 indels), few of them consistent with a flipped allele
  (eaf_concordance_decode_vs_decodeme.py).

## Parallel window fetches (step 2, after the fix)

The first 24 SMP proteins built in 1,116 s; every one has 0 of 1,171,498 index rows missing.
One 1 Mb window (chr11:116.3-117.3 Mb, 220,632 rows in all) from each of the 24 through
ParallelWindowFetcher, two runs each: max_parallel 4 20.1 / 18.8 / 19.7 / 20.7 s; 6 13.9 / 12.9
/ 14.2 / 17.1 s; 8 12.9 / 11.8 / 11.6 / 12.0 s; 12 8.7 / 11.7 s; 16 10.9 / 12.7 s. Peak RSS 1.7
GB. The gain stops at about 8: DECODE_MAX_PARALLEL_FETCHES = 8.
