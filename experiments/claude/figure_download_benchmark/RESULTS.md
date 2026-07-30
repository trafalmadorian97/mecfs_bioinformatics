# Figure-download benchmark results

Run: 2026-07-24, 137 unique manifest assets, 3 reps each, into a fresh temp
dir per run. Wall-clock seconds (log: `experiments/claude/logs/figure_download_benchmark.log`).

| config                | median | min  | max  |
|-----------------------|-------:|-----:|-----:|
| per_asset_parallel_w8 |   18.1 | 16.4 | 19.7 |  ← current production path
| chunked_c20_w4        |    3.0 |  2.9 |  3.0 |  ← ~6x faster
| chunked_c20_w1        |    9.8 |  9.5 | 11.7 |
| chunked_all_w1_bulk   |    5.1 |  4.9 |  5.7 |

## Findings

1. **The chunked batch is far faster, not slower.** My earlier prediction
   (bulk/serial would lose to 8-way parallelism) was wrong in effect. Even the
   *fully serial* chunked run (c20_w1, 7 sequential batches) beats the 8-way
   parallel per-asset path (9.8s vs 18.1s).

2. **Reason: per-invocation overhead dominates.** Each `gh release download
   --pattern <one>` pays `gh` process startup + a full release-asset-list
   lookup just to locate one blob. The current path pays that 137 times.
   Chunking amortizes it: c20 pays it ~7 times, bulk pays it once.

3. **gh has no internal download parallelism.** The single bulk call (5.1s) is
   slower than 4 parallel chunks of 20 (3.0s), so within one invocation gh
   downloads its matched assets serially. Splitting into a few parallel chunks
   recovers concurrency on the actual transfers while keeping listing overhead
   near-zero.

4. **API pressure — the thing that provokes the HTTP 500s — drops sharply.**
   ~7 asset-list lookups (c20_w4) instead of 137, so far less burst against
   GitHub.

## Recommendation

`chunked_c20_w4` is the clear winner on every axis: ~6x faster wall-clock,
~20x fewer asset-list API calls, and with per-chunk retry the failure blast
radius is one ~20-asset chunk (re-downloaded with `--clobber`) rather than the
whole run. Manifest-membership is still exact (one `--pattern <hash>` per
manifest asset), so the CI consistency check is preserved.

Chunk size / worker count are mild knobs; c20/w4 is a good default. Larger
`chunk_workers` reintroduces the burst we're trying to avoid, so keep it small.
