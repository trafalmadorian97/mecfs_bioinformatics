# Python & coding conventions

General style rules for new code in this repo. These favor type safety, self-documenting
call sites, and making invalid states unrepresentable.

## Paths

- **Represent filesystem paths as `pathlib.Path`, not `str`** — including dict keys, function
  parameters, return types, and in-memory data structures. Convert to `str` (usually via
  `.as_posix()`) only at the last moment, when serializing to JSON / text / another format that
  requires a string. Mixing `str` and `Path` for the same conceptual value invites silent bugs
  and forces every consumer to re-derive POSIX semantics.
  - For directory-coverage checks prefer `Path.parents` over string-prefix tricks:
    `child == parent or parent in child.parents`. This rejects false-prefix matches automatically
    (`Path("foo_extra.html")` is not under `Path("foo")`).
  - When loading a path from JSON, convert with `Path(s)` at the deserialization site, not at
    every downstream call site.
- **Use `PurePath` for a path that is relative to a currently-unspecified base directory** (e.g.
  a location relative to a package dir). `PurePath` gives type safety while correctly withholding
  filesystem operations (`write_text`, `unlink`, `is_file`) that only make sense once concretized.
  It signals "join me onto a base before use." Concretize via `Path(base) / the_purepath`.
  Example: `GWASLAB_HAPMAP3_HG38_SNPLIST_RELPATH = PurePath("data/...")`, resolved at use as
  `Path(os.path.dirname(gl.__file__)) / GWASLAB_HAPMAP3_HG38_SNPLIST_RELPATH`.

## Types & invariants

- **Type a string parameter with a small fixed value set as a `Literal`, not bare `str`.** Define
  a named alias once, near the values it describes, and reuse it on both the attrs field and the
  `create()` param. Makes invalid states unrepresentable at the type level, documents allowed
  values at the call site, and lets the typechecker catch typos. Example:
  `NanPolicy = Literal["raise", "ignore", "mean", "zero"]`.
- **Make invalid states unrepresentable via `__attrs_post_init__`.** When an attrs class has
  fields whose validity depends on each other, enforce the relationship in `__attrs_post_init__`
  rather than leaving an attribute that is sometimes silently unused. If a value can be ignored in
  some configurations, the type is too loose — widen it to admit the absent case (`str | None`)
  and assert the cross-field rule at construction.
- **Enforce numpy shape/dtype invariants in `__attrs_post_init__`.** When an attrs class holds
  numpy arrays that must satisfy shape/dtype invariants (all 1-D and the same length, a specific
  dtype kind), assert them — loop over `(name, arr, expected_kind)` tuples asserting `arr.ndim`,
  `arr.shape[0] == expected_length`, and `arr.dtype.kind == expected_kind` with clear messages.
  Comments describing dimensionality/dtype drift and are not enforced.
- **Shift left: assert and fail fast on invalid input.** When a code path depends on a
  precondition (e.g. "the input frame already carries an N column"), do not silently assume it —
  assert it (or explicitly raise) as early as possible, with a message naming what was expected vs.
  seen. Surfacing a problem at the boundary is cheaper than letting invalid input flow downstream
  and fail somewhere confusing, or produce wrong results silently.

## Function signatures

- **Don't return bare/anonymous tuples** (`-> tuple[float, float, float]`). Return an `@frozen`
  attrs class (or a `NamedTuple` for many fields) whose attributes have meaningful names, so the
  call site is self-documenting and positional-swap bugs are impossible. Define the small class
  next to the function. The one acceptable exception is a nested/local helper whose tuple is
  consumed a line or two later, where the reader sees both ends at once.
- **Pass same-typed parameters as named kwargs.** When a call passes two or more args of the same
  type (several `str` or `Path` params), name them: `f(source_url=..., uri=...)`, not `f(a, b)`.
  Positional args of identical type can be transposed silently — the typechecker won't catch it.
- **Prefer module-level helper free functions over helper methods.** When a class (e.g. a Task)
  accumulates several `def _helper(self, ...)` methods, hoist them to module-level free functions
  that accept only the fields they actually use (pass `config`, `index_task`, `fetch`, not `self`).
  A helper method implicitly receives all of `self`; a free function's signature reveals its real,
  narrow dependencies. Keep on the class only what genuinely needs the whole object (`execute`,
  `create`).

## Dataframes & I/O

- **Prefer polars over pandas for new dataframe code** (`pl.DataFrame` / `pl.LazyFrame`, `scan_*`).
  This is the approved default for new work (performance, lazy scans, cleaner expression API) — not
  a mandate to rewrite existing pandas (e.g. the genomic_sem port can stay pandas).
- **When converting a narwhals frame to polars, use `.to_polars()`, not `.to_native()`.**
  `to_native()` returns whatever backend narwhals wraps (pandas, pyarrow, polars), so polars-only
  downstream code breaks whenever the backend isn't polars. `to_polars()` always yields a genuine
  polars DataFrame.
- **Use shared string constants for common dataframe column names, not repeated literals.** Put
  each constant in the file matching the source/tool that *defines* that column. Known homes:
  `mecfs_bio/constants/gwaslab_constants.py` for gwaslab output columns (strictly gwaslab's own),
  `ldsc_constants.py` for LDSC-internal names, `regenie_constants.py` for raw regenie columns,
  `ppp_index_constants.py` for PPP-variant-index columns. Create a new file when a source has no
  home. Only genuinely one-off local id columns may stay as literals.

## Subprocess

- **Prefer `execute_command` (`mecfs_bio/util/subproc/run_command.py`) over raw
  `subprocess.run`/`Popen`.** It's the repo's standard wrapper — streams+logs output via structlog,
  asserts args are strings, returns combined stdout+stderr, and raises `CalledProcessError` on
  non-zero exit. `execute_command_with_retries` adds injectable backoff for flaky remote services.
- For a **new flaky remote call**, use `call_with_retries`
  (`mecfs_bio/util/retry/call_with_retries.py`) — the repo's single retry-with-backoff loop. It
  takes a no-arg callable (bind arguments with `functools.partial`) and an explicit `retry_on`
  tuple so exceptions outside it propagate immediately and genuine bugs are never retried. Do not
  write another backoff loop.

See also [testing.md](testing.md) and [build-system-and-tasks.md](build-system-and-tasks.md).
