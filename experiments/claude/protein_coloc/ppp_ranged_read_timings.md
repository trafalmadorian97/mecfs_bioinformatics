# UKB-PPP ranged reads: timing results

Measured 2026-11-14 from a WSL2 laptop on a home connection, against the pre-signed S3 URLs
that Synapse serves for UKB-PPP protein tars (EUR discovery manifest). Every number is one
run, and repeat runs vary by up to a factor of 3, so treat them as rough.

Scripts (in this directory):

- probe_ppp_tar_range_layout.py, in the main checkout only and not committed: the original
  layout probe.
- verify_ppp_ranged_window_fetch.py: compares the ranged PppWindowFetcher with the full
  download, end to end.
- time_ppp_ranged_read_phases.py: times each phase.

## Tar layout (A1BG, OID30771)

- An uncompressed tar of about 549 MB, holding a directory entry and 23 gzipped regenie
  files, one per chromosome.
- Chromosomes are not stored in numeric order (the tar starts chr18, chr19, chr9, chr11, ...).
- Chromosome files range from 7.2 MB (chr21) to 45.3 MB (chr2), so 1-8% of the tar each.
- The files are plain gzip, not BGZF, so a chromosome file cannot be read partway; it is
  always read whole.
- The tars are GNU format. Names over 100 characters go in a separate ././@LongLink entry
  before the file. 573 of the 2,940 manifest proteins have such names (complexes, e.g.
  AMY1A_AMY1B_AMY1C).
- The pre-signed URL answers Range requests with HTTP 206.

## End to end: ranged fetch vs full download

AMY1A_AMY1B_AMY1C, a 1 Mb window each on chr1 and chr22. Both methods gave identical rows
(5,431 and 5,499).

| Version | Ranged fetch | Full download + extract |
|---|---|---|
| One new connection per range read | 10.0 s | 14.5 s |
| One session per thread (current) | 6.5 s | 16.9 s |

The full download is fast on this connection (550 MB in about 12-15 s, about 40 MB/s), so the
ranged fetch saves more network traffic (about 15 times less) and disk space than time.

## Header walk: new connection per read vs keep-alive session

A1BG, 23 files plus the directory entry, about 25 reads that must run one after another
because each header gives the position of the next.

| Run | New connection per read | One requests.Session |
|---|---|---|
| 1 | 5.50 s | 2.18 s |
| 2 | 5.14 s | 2.50 s |
| 3 | 5.31 s | 2.45 s |

Opening an HTTPS connection costs several network round trips. Reusing the connection leaves
roughly one round trip per read.

Production reader (SharedClientSynapseDownloader, one session per thread), three other
proteins in turn through one downloader:

| Protein | Open (Synapse lookup of the URL) | Header walk |
|---|---|---|
| OID30236 (first; includes Synapse login) | 0.88 s | 3.14 s |
| OID21311 | 0.17 s | 2.75 s |
| OID30146 | 0.21 s | 3.18 s |

Reusing connections from one protein to the next saves at most one connection setup per
protein, which is too small to see in this noise.

## Reading one chromosome file: 1 vs 4 vs 8 parallel range requests

A1BG chr2 (45.3 MB), split into equal parts read concurrently. Runs 2 and 3 alternate the
order of the connection counts.

| Connections | Run 1 | Run 2 | Run 3 |
|---|---|---|---|
| 1 | 2.71 s | 4.17, 1.48, 3.66 s | 1.84, 1.12, 1.50 s |
| 4 | 4.01 s | 1.46, 2.80, 0.94 s | 1.70, 1.66, 1.26 s |
| 8 | 2.66 s | 2.55, 2.99, 1.37 s | 1.33, 1.62, 1.79 s |

Speeds ranged from 11 to 48 MB/s whatever the number of connections, and splitting the file
showed no consistent gain. At best it saves about a second on the largest chromosome.

## Conclusions

- Ranged reads replace the 550 MB tar download with about 25 small header reads plus the
  needed chromosome files. On this connection that is faster, and it uses far less network
  traffic and no scratch disk.
- Connection reuse was the worthwhile tuning, and it is implemented as one requests.Session
  per thread in SharedClientSynapseDownloader. It cut the header walk by about 2.5-3 s per
  protein.
- Splitting one file across parallel connections is not worth the complexity.
- The bigger remaining gain is outside the reader: TraitProteinColocTask fetches proteins one
  at a time (_test_selected_pairs). Most of each fetch is spent waiting on the network, so
  fetching several proteins at once would overlap those waits.

## Fetching several proteins at once

probe_parallel_ppp_window_fetch.py: PppWindowFetcher.fetch_protein (ranged reads, one session per
thread) over batches of 12 proteins, one 1 Mb window each on a fixed mix of 12 chromosomes,
with a thread pool of 1, 2, 4 or 8 workers. Every batch uses different proteins, and each
worker count runs twice in the order 1, 4, 2, 8, 8, 2, 4, 1.

| Workers | Wall time per 12 proteins (two runs) | Per protein | Speedup | CPU time / wall time | Peak RSS |
|---|---|---|---|---|---|
| 1 | 77.0, 65.4 s | 5.93 s | 1.00x | 0.20, 0.21 | 546, 847 MB |
| 2 | 34.9, 34.9 s | 2.91 s | 2.04x | 0.35, 0.38 | 662, 979 MB |
| 4 | 21.3, 21.2 s | 1.77 s | 3.35x | 0.76, 0.63 | 796, 908 MB |
| 8 | 17.9, 13.4 s | 1.30 s | 4.56x | 0.83, 1.11 | 932, 1042 MB |

The process RSS was 247 MB before the first batch. RSS rarely shrinks once memory is
allocated, so later batches inherit earlier peaks. Even so, 8 workers stay near 1 GB.

- One worker at a time spends about 80% of each fetch waiting on the network (CPU time /
  wall time 0.2).
- Speedup is near-linear to 2 workers and still 3.4x at 4. By 8 workers, CPU time / wall
  time reaches about 1, so the work becomes CPU-bound (gzip decompression and CSV parsing
  under Python's GIL) and further workers add little.

## How much of a coloc run is fetching

The coloc task logs "fetched dense protein windows" once per protein, but nothing else in a
protein's cycle carries a timestamp. So the gap between consecutive messages is the previous
protein's pair tests plus the next protein's fetch. All three runs below used the old
full-download fetcher.

| Log (experiments/claude/logs) | Proteins | Median gap | Min gap | Max gap |
|---|---|---|---|---|
| coloc_gnomad_v4_1_rerun_20261104_1617.log | 98 | 17 s | 14 s | 149 s (between traits) |
| willer_ldl_ppp_coloc_run.log | 30 | 24 s | 19 s | 29 s |
| decode_me_ppp_coloc_run.log | 8 | 19 s | 18 s | 26 s |

The old full download plus extraction measured 14.5-16.9 s today. If that held for these
runs, the old fetch was roughly 70-85% of each protein's time, and testing was about 3-8 s
per protein. This is an inference from today's timings, not a measurement from the runs.

With the ranged fetcher run serially (about 6 s per protein), fetching should fall to about
45-65% of each protein's time. Fetching 4 proteins at once (about 1.8 s per protein) would
make the per-protein test time, not the fetch, the larger cost.

Earlier probes in experiments/claude/logs (in the main checkout) put this connection's
bandwidth cap at roughly 30-60 MB/s:

- synapse_max_threads_sweep.log: 30-42 MB/s for single PPP tar downloads.
- parallel_synapse_probe_n16.log: 49 MB/s across 16 concurrent PPP tar downloads.
- s3_concurrency_sweep.log: 30-60 MB/s across many connections to a 1000 Genomes S3 object.

At about 25 MB of compressed data per protein, 8 parallel fetches every 1.3 s come to about
19 MB/s, below those caps. So bandwidth is not yet the limit; CPU is.

## Implemented: fetch_all

The PPP coloc fetcher is now a ParallelWindowFetcher (dense_window_fetch.py) over
PppWindowSource (formerly PppWindowFetcher), fetching PPP_MAX_PARALLEL_FETCHES = 4 proteins at
once; the coloc task no longer fetches one protein at a time. verify_ppp_ranged_window_fetch.py checks it live: 6 proteins took 16.1 s
through fetch_all versus 42.7 s one at a time, with identical frames in request order.
