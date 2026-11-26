# Build system and Tasks

The build system is described in `docs/Codebase_Concepts/Build_System.md`; its central class is
the `Task`. This page covers how to write Tasks, where code belongs, and how the cache and asset
store behave.

## Layering: where code goes

`mecfs_bio` is layered `analysis` > `assets` / `asset_generator` > `build_system` > `util` >
`constants`. A layer imports only from layers below it; import-linter enforces this
(`.importlinter`).

- **Classes go in `build_system`; the specifics that instantiate them go in `assets`.** Tasks and
  configuration classes (e.g. a release descriptor) live in `build_system`. URLs, URL templates,
  per-release constants, and chosen column lists live in `assets`, next to where the class is
  instantiated. Never hardcode a URL in `build_system`.
- **A vocabulary that both layers type against** (e.g. a `Literal` of group codes) goes in
  `constants`.
- **No repo-relative paths in `build_system`.** A file name or path relative to the repo root (a
  config file, a generated YAML) is defined by the caller that knows the repo layout, usually
  `mecfs_bio/analysis/runner/default_runner.py` or `tasks.py`, and passed into `build_system` code
  as a parameter.

## Writing a Task

- **`.create()` derives identifying metadata from its inputs.** Don't accept `trait`, `project`,
  or similar fields as parameters. Read the input task's `.meta`, assert its type, and build the
  output meta from `meta.trait` and `meta.project`. This keeps an output's labels consistent with
  its source data.
- **Production Tasks use a domain `Meta` subtype.** `SimpleFileMeta` and `SimpleDirectoryMeta`
  are for tests only. Pick a subtype from the `Meta` union in `mecfs_bio/build_system/meta/meta.py`;
  for reference data use `ReferenceFileMeta`, modeled on `GetUniProtReferenceDataTask`.
- **Write `execute()` as a high-level specification of the algorithm.** Each line is one named
  step (fetch inputs, pass A, fit, write outputs) that delegates to a module-level helper. No raw
  loops in `execute()`: wrap a per-item loop in a helper. Details such as JSON indentation belong
  in helpers.
- **Read tabular inputs with `scan_dataframe_asset(asset, task.meta)`**
  (`mecfs_bio/build_system/meta/read_spec/read_dataframe.py`), not `pl.read_csv` with a hardcoded
  separator. It uses the asset's `read_spec`, so the reader doesn't depend on the producer's
  on-disk format. If the upstream meta has no `read_spec`, add one. Multi-file directory assets
  stay path-based.
- **Docstrings of a Task or asset module say only what is specific to it.** If a sentence would be
  equally true of any Task, leave it out.

## External data and reproducibility

- **Acquire every input through a Task.** Builds must work on a fresh machine. Don't read from
  ad-hoc local caches such as `~/.gwaslab/`; wrap the download as a Task. Exception: a file bundled
  inside a pinned Python dependency may be resolved from the installed package path.
- **Run external CLI tools through pixi** (`pixi run bcftools ...`), not a binary assumed to be on
  `PATH`.
- **Normalize external data once, as its own asset**, next to the download. This is usually a
  `PipeDataFrameTask` built from the generic pipes in `mecfs_bio/build_system/task/pipes/`; check
  for an existing pipe before writing one. Analysis Tasks consume the normalized table. They never
  split delimited lists, drop source-specific row kinds, or rename raw columns themselves.

## How the cache behaves

- **Task code is not part of the cache key.** An asset is up to date when its file hash and its
  dependencies' hashes match the recorded trace (`verify_trace` in
  `mecfs_bio/build_system/rebuilder/verifying_trace_rebuilder/verifying_trace_rebuilder_core.py`).
  Changing how a Task writes its output therefore triggers no rebuild. Judge such a change by
  whether old and new output are semantically equivalent, not by recomputation cost. The real
  risk is drift: stored assets keep the old bytes until something else forces a rebuild.
- **Building a target materializes all of its transitive dependencies.** A hand-deleted
  intermediate is rebuilt the next time anything downstream is built, so deleting intermediates
  does not save disk. To keep an intermediate off the main disk, give it its own meta `sub_folder`
  and move it with a `path_remap` rule.

## The asset store

- **Traces are location-independent**, so moving or copying assets never invalidates the cache.
- **The store can be split across disks** with a `path_remap` section in the gitignored,
  machine-local `default_runner_config.yaml`. `pixi r invoke migrate-asset-store` (a dry run by
  default) shows the declared rules. Remap only large subtrees with few files: a secondary mount's
  cost is per-file latency, not throughput.
- **To find an asset on disk, compute its path; don't `find` the store**, which is huge and slow.
  Take the defining Task's `.meta`, call `simple_meta_to_relative_path`
  (`mecfs_bio/build_system/rebuilder/metadata_to_path/simple_meta_to_path.py`), and join the result
  onto the root (and any `path_remap`) from `default_runner_config.yaml`.

## Memory budget: 16 GB

The main analysis must run on a 16 GB laptop, and CI's standard 16 GB GitHub runner enforces this
on purpose. **When a system test runs out of memory in CI, fix the memory regression.** Don't
propose a larger runner or swap.
