# Unit-test speed analysis

Measurements behind the test-speed changes in PR #1359. All timings come from a 16-core WSL2
machine and use the full unit suite (test_mecfs_bio/unit), with no testmon and with typeguard on
unless stated otherwise. Another job was sometimes running in the background, so expect a few
seconds of noise in the wall times.

## Results at a glance

| Command | main | this branch |
|---|---|---|
| Serial (invoke test without testmon) | 142s | ~76s |
| xdist, 8 workers, 2 threads per worker (invoke test-xdist) | n/a | ~22–26s |

What changed, from largest to smallest effect:

1. **pytest-xdist** (invoke test-xdist, now used by invoke green). It runs min(8, usable cores / 2) workers.
2. **gc.freeze() after test collection** (test_mecfs_bio/unit/conftest.py).
3. **cached_importr**, which builds rpy2 package wrappers once per process.
4. **Thread caps for xdist workers**: BLAS, OpenMP, numexpr and polars get 2 threads each.
5. **Shared module-scoped fixtures** for the SuSiE and vectorized S-LDSC tests.
6. **Removed test_gwaslab_genetic_corr_by_ct_ldsc_task.py**, because genetic correlation is moving to the
   pure-Python implementation.

## 1. Baseline profile

- The first serial run took 152s. Three TwoSampleMR tests failed only because that R package was missing
  from the fresh worktree; after invoke install-r-packages they passed. The later main baseline
  measured 142s, with the same three failures.
- A few tests dominate. The 20 slowest tests account for about 55% of test time, and 81% of tests take
  under 0.1s.
- The slowest modules were SuSiE fine-mapping (~24s), the gwaslab ct-LDSC genetic correlation test,
  vectorized S-LDSC (~16s) and gwaslab sumstats creation (~6s).

## 2. pytest-xdist

| Workers | Full-suite wall time (before the other changes) |
|---|---|
| 4 | 79s |
| 8 | 50s |
| 16 | 86s |

- **Why 16 workers is slower:** every worker collects the whole suite, about 9s warm, so that cost is paid
  16 times over. The workers also fight over the CPU (see section 4).
- **Lazy imports would not shrink collection:** scipy, pandas and gwaslab are needed on every worker. Import
  self-time (6.2s total) is spread across hundreds of packages.
- **Worker count rule:** since 8 was best on 16 cores, test-xdist uses min(8, usable cores / 2), with a floor of 1. It counts
  usable cores with os.process_cpu_count().
- **testmon still selects the right tests under xdist:** after a probe edit to SusieRFinemapTask.execute, xdist
  and serial runs selected the identical 13 tests.

## 3. gc.freeze() after collection

- **The cost:** gwaslab, the runner's post-execute gc hook, and the SuSiE and vectorized S-LDSC tasks call
  gc.collect(). A full collection walks every tracked object. In the test process that includes every
  module the suite imports, so each call costs about 0.28s. A serial run makes 55 such calls, about
  16s in total.
- **The fix:** calling gc.freeze() in pytest_collection_finish moves the roughly 359,000 objects alive after collection
  into a permanent generation that later collections skip. Production gc behavior is unchanged.

Measured with a pytest plugin that times every gc.collect() call (two runs each):

| | Without freeze | With freeze |
|---|---|---|
| Serial wall time | 96.9s, 88.9s | 75.8s, 76.3s |
| Serial time in gc.collect() | 15.7s, 15.1s | 1.0s, 1.0s |
| xdist wall time | 32.0s, 29.7s | 22.1s, 23.1s |

The freeze ran in all 8 xdist workers, which I confirmed with a per-worker probe.

## 4. Thread oversubscription under xdist

Without caps, each of the 8 workers starts a full-size BLAS, OpenMP and polars thread pool and the workers
oversubscribe the CPU. Capping OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS,
NUMEXPR_MAX_THREADS and POLARS_MAX_THREADS:

| Threads per worker | Run 1 | Run 2 |
|---|---|---|
| Uncapped | ~39s | |
| 1 | 33.7s | 28.1s |
| 2 | 29.8s | 28.8s |
| 4 | 32.5s | 28.7s |

The capped settings are within noise of each other. I chose 2 because its runs were the most consistent.
invoke test (serial) does not set the caps.

## 5. rpy2 importr caching

- **The cost:** rpy2's importr rebuilds Python wrappers for every exported R function on each call, even when
  the R package is already loaded. That costs about 0.6s for ggplot2 and about 0.12s for susieR.
- **The fix:** cached_importr (mecfs_bio/util/r_packages/cached_importr.py) wraps it in functools.cache. This is safe
  because each process has exactly one embedded R session. It is used by the SuSiE, LAVA, TwoSampleMR
  and RData-extraction tasks.
- **Checks:** the LAVA vignette system test and the TwoSampleMR unit tests pass with it.

## 6. SuSiE tests

The module dropped from 23.9s to 15.0s (serial, median of 3):

- Making susie_prerequisite_file_tasks module-scoped (using tmp_path_factory) took it from 23.9s to 21.6s.
- cached_importr then took it from 21.6s to 15.0s.

**Shrinking the synthetic problem is not worth it.** I swept the problem size, n individuals by m variants,
with 25 seeds for each causal configuration (two causal variants, or none). Each run used the real
SusieRFinemapTask and the test's own pass criteria:

| n | m | Pass rate | Median task time |
|---|---|---|---|
| 2500 | 100 (current) | 100% | 0.81–0.82s |
| 500 | 100 | 100% | 0.83–0.84s |
| 1000 | 50 | 100% | 0.68–0.72s |
| 500 | 30 | 100% | 0.60–0.63s |

- **Timing:** every size passed every time, so the signal is far from the threshold. Even the smallest problem
  saves only about 0.2s per run, because fixed R and task overhead dominates.
- **Reliability:** the margin is thinner in the null case. The worst maximum PIP was 0.10 at the current size and
  0.82 at n=500, m=50, so smaller problems leave less room for error.
- **Plot DPI:** lowering the ggsave DPI from 300 to 100 saves only about 0.15–0.25s per plot.

## 7. Vectorized S-LDSC tests

A module-scoped fixture now writes the fixture files, runs the two Prepare tasks (baseline and
annotation set) and pickles the clean sumstats once per module. Each test still gets its own store
under its own tmp_path, and tests that alter the sumstats pickle their own frame. The module dropped
from 15.6s to 10.1s (serial median), and all 9 tests pass.

## 8. Where the remaining serial time goes

A cProfile of one full serial run (before gc.freeze), split by category of the code doing the work:

| Theme | Approx. time |
|---|---|
| polars LazyFrame.collect (real dataframe work, ~13,000 calls across many tests) | ~22s |
| gc.collect() (since fixed by gc.freeze) | ~16s |
| Test collection and imports | ~16s |
| R via rpy2, mostly SuSiE | ~13s |
| Typeguard instrumentation | ~2s |
| pytest internals | ~2s |
| File hashing in the build system | ~0s |

The profile also attributed about 71s to matplotlib text measurement. That is a profiler artifact:
category times summed to 179s against 106s of wall time, and no plotting test is slow in the per-test
durations.

## 9. Filesystem access is not a significant cost

Many tests run Tasks that read and write files. I reran the suite with all pytest temp dirs (basetemp
and TMPDIR) on tmpfs, which lives in RAM, instead of the ext4 disk:

| | Disk | tmpfs |
|---|---|---|
| Serial | 96.9s, 88.9s | 90.6s, 95.8s |
| xdist, 8 workers | 32.0s, 28.9s | 27.5s, 27.7s |

The serial difference is within noise, and the xdist difference is at most about 2s. The test files are
small and stay in the page cache. In the profile, raw file calls came to about 7s of a 106s run, mostly
Python reading and parsing text rather than waiting on the disk.

## Remaining opportunities

- The next-largest serial module is gwaslab sumstats creation, at about 5.6s. Every other slow test takes a
  few seconds or less.
- Collection costs about 9s per process and is mostly unavoidable imports.
- polars work is real computation spread over many tests, so there is no single target.
