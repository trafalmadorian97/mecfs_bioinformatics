# Build-system overhead on large targets: profiling results (2026-11-11)

## Question

Iterating on targets with thousands of transitive dependencies, such as the Willer LDL CSF
colocalization screen (`mecfs_bio/assets/gwas/ldl/willer_et_al/analysis/willer_ldl_csf_coloc.py`),
was slow even when nothing needed rebuilding.  Was the cost in verifying the assets (hashing),
or somewhere else?

## Method

`mecfs_bio.analysis.csf_coloc_analysis.run_willer_ldl_csf_coloc()` was run under `cProfile`
against the production `DEFAULT_RUNNER` (ImoHasher xxh3-128 tracer, 9.2 MB YAML info store).
Every asset was already built, so each run was a pure verification pass: 7,047 assets
verified, 0 materialized, in every run below.  Timings are under cProfile, which inflates
code that makes many small Python calls (notably `copy.deepcopy`), so absolute numbers are
higher than in a normal run, but the ranking is unambiguous.

## Results

| Stage | Baseline | + Kahn scheduler | + CSafeLoader | + conditional save |
|---|---|---|---|---|
| Total run | 872 s | 30.8 s | 13.1 s | **7.8 s** |
| Scheduler loop (verify 7,047 assets) | 843 s | 4.6 s | 7.0 s | 4.0 s |
| Info-store YAML load | 24 s | 23 s | 3.0 s | 2.8 s |
| Final info-store save | ~3 s | 1.9 s | 1.9 s | skipped |

The scheduler-loop time varies by a few seconds between runs that did not change it
(most likely filesystem caching).

### Baseline breakdown

- `_update_frontier_and_G` in `topological_scheduler.py`: **817 s (94%)**.  It
  `deepcopy`-ed the whole dependency graph after every node, which is O(N²) for N tasks
  (~274 M recursive `deepcopy` calls).
- PyYAML `FullLoader` (pure Python) loading the info store: 24 s.
- ImoHasher, 21,127 calls: 13 s total (~0.6 ms per call).  Hashing was **not** the
  bottleneck: ImoHasher samples three 16 KB blocks per file rather than reading whole files.
- `stat` / asset-existence checks: ~6 s.
- Unconditional final YAML save: ~3 s.

## Fixes (branch scheduler-kahn-frontier)

1. **Kahn's-algorithm frontier**: track each node's count of unfinished dependencies and
   decrement on completion, instead of copying and pruning the graph.  O(out-degree) per
   node; visit order unchanged.
2. **CSafeLoader** for the info store.  It is written with CSafeDumper, so safe loading is
   lossless.  On the real 9.2 MB store, `FullLoader` and `CSafeLoader` produced identical
   data; uncontended load time was 8.8 s vs 1.1 s.
3. **Skip the final save** when the store existed at the start of the run and its
   trace_store is unchanged at the end.

## Conclusion

A no-op run on a 7k-asset target went from ~14.5 min to ~8 s (profiled).  The remaining
time is ~4 s of per-asset verification (hashing + stat) and ~3 s of YAML load; neither
justifies a format change or an "assume these assets are up to date" escape hatch for now.
