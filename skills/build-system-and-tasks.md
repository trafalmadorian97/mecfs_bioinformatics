# Build system & Task conventions

The build system is described in `docs/Codebase_Concepts/Build_System.md`; the central class is the
`Task`. This page collects the conventions for writing Tasks, plus two important properties of the
cache and the asset store that affect how you reason about changes.

## Writing a Task

- **`.create()` derives identifying metadata from input deps — it does not accept it.** Do not
  take `trait`, `project`, or similar identifying fields as direct parameters. Instead read them
  from an input task's `.meta`, asserting the expected meta type first. This is standard across the
  codebase (`UpSetPlotTask.create`, `GWASLabManhattanAndQQPlotTask.create`,
  `GeneticCorrelationClustermapTask.create_std_with_clustering`). Letting callers pass `trait` /
  `project` directly invites drift between an output's labels and its source data's labels.
  - Inside `create()`, read `some_task.meta`, `assert isinstance(meta, ExpectedMetaType)` (e.g.
    `ProcessedGwasDataDirectoryMeta`, `ResultTableMeta`, `GWASLabSumStatsMeta`,
    `ResultDirectoryMeta`), then build the output meta from `meta.trait` / `meta.project`.
  - When the primary input is wrapped in a source object, expose `trait` / `project` as properties
    on the source class that assert-and-extract internally.
- **`SimpleFileMeta` / `SimpleDirectoryMeta` are test-only.** Production Tasks must pick a
  domain-appropriate `Meta` subtype — the Meta carries the asset's organizational identity
  (group / sub_group / sub_folder / extension / read_spec) the build system uses to place and read
  the asset; the "Simple" variants' bare id is a test convenience that loses that structure.
  - For reference data, use `ReferenceFileMeta(group, sub_group, sub_folder: PurePath, extension,
    id, read_spec?)`; model the Task on `GetUniProtReferenceDataTask` (a no-deps source task with a
    `meta: Meta` field and a `.create()` that constructs the meta internally). Consult the `Meta`
    union in `meta/meta.py` for the right subtype.
- **Read tabular inputs via `scan_dataframe_asset`, not `pl.read_csv` with a hardcoded separator.**
  Use `scan_dataframe_asset(asset, task.meta)` (from
  `mecfs_bio/build_system/meta/read_spec/read_dataframe.py`), which consults the asset's metadata
  `read_spec` (`DataFrameReadSpec` with a text / whitespace-sep / parquet format) to decide how to
  parse, and returns a narwhals LazyFrame — collect with `.collect().to_polars()`. This decouples
  the downstream reader from the upstream's on-disk format (tsv, space-delimited, or parquet, read
  the same way); hardcoding a separator bakes in a format assumption that silently breaks if the
  producer changes. If the upstream has no read_spec, add one to its meta. Exception: multi-file
  directory assets (e.g. per-chromosome LD scores) aren't a single dataframe — keep those
  path-based.
- **Acquire dependencies explicitly; minimize assumptions about the runtime machine.** Builds must
  be reproducible on a fresh machine, not just the developer's box.
  - Do not read files from ad-hoc local caches like `~/.gwaslab/`. Wrap acquisition as a task (e.g.
    wrap gwaslab `download_ref("1kg_eur_hg38")`; helper at
    `mecfs_bio/build_system/task/gwaslab/gwaslab_util.py`).
  - Invoke external CLI tools through pixi (`pixi run bcftools ...`), not a bare binary on PATH.
  - Exception: a file bundled inside a pinned pyproject.toml package dependency is guaranteed
    present, so resolving it from the installed package path (e.g. gwaslab's bundled hapmap3 snplist
    under `Path(gwaslab.__file__).parent / "data"`) is fine — no rehosting needed.

## The build cache does NOT see code changes

The cache decides an asset is up to date by asset id, comparing the hash of the on-disk file and
its dependencies' hashes against recorded values. **Task implementation code is not part of the
trace** (see `verify_trace` in
`mecfs_bio/build_system/rebuilder/verifying_trace_rebuilder/verifying_trace_rebuilder_core.py`).

Consequences:

- Refactoring how a task writes its output (swapping a pandas writer for polars, changing parquet
  encoding) causes **no** rebuild and no invalidation. **Do not argue against a refactor on the
  grounds that it will force expensive recomputation — that reasoning is wrong here.**
- The real consequence of an output-writing change is *drift*: stored assets keep the old bytes
  until something independently forces a rebuild, at which point the new writer produces different
  bytes and downstream rebuilds then. Assess such a change on whether new and old output are
  semantically equivalent, not on cache cost.
- Corollary gotcha: because traces are byte-content hashes, a **polars/pandas version bump can
  silently invalidate traces** of affected assets even when the data is unchanged (parquet row
  order, pandas Categorical code ordering). See [environment-gotchas.md](environment-gotchas.md).

## Asset store architecture

- Traces are location-independent content hashes, so relocating or copying assets never invalidates
  the cache — a moved asset trace-verifies at its new root with zero rebuilds.
- The asset store can be **split across disks via a `path_remap` section in the gitignored,
  machine-local `default_runner_config.yaml`** — so the current split is not derivable from the
  repo; check the file, or run `pixi r invoke migrate-asset-store` (dry run by default) to see the
  declared rules.
- **Selection rule for what to remap onto a slower/secondary mount: large, FEW-file, rarely-read
  subtrees only.** The cost of a secondary mount (e.g. a Windows drive under WSL) is per-file
  latency, not throughput — so a 51G tree spread over 186k files stays local, while a 13G tree in
  179 files is a fine candidate.

## The 16GB runner is an intentional constraint, not a CI problem

System tests run on the standard 16GB `ubuntu-latest` GitHub runner **deliberately**: the design
intent is that the main analysis be runnable on a 16GB laptop, and CI is the regression guard for
that property. **When a system test OOMs in CI, do not propose a larger runner (or re-enabling
swap) as the fix** — the OOM is a true regression for end users; find and fix the memory
regression. Suggest a larger runner only if explicitly told the constraint was relaxed for that
test.

See [python-style.md](python-style.md) for `PurePath`, `Literal`, and column-constant conventions
referenced here, and [testing.md](testing.md) for how to test Tasks.
