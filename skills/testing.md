# Testing conventions

The guiding principle: a test should fail only when observable behavior changes, never when a
behavior-preserving refactor moves internals around. Tests that break under refactoring impose
cost without catching bugs.

## Test at the Task level (the public API)

- **Unit-test the Task, not its private helpers.** The `Task` is the unit a user of the build
  system operates on, so it is the public API (in the "test via public APIs" sense from *Software
  Engineering at Google*). Refactoring a Task's internals while preserving behavior should not
  break its tests.
- Default to constructing the Task and calling `execute` with `FakeTask` dependencies, synthetic
  input dataframes written to `tmp_path`, and a `fetch` that maps asset id → Asset. See
  `test_ppp_protein_rg_task.py` and `test_consolidate_ld_scores_task.py`.
- Reserve helper-level tests for logic genuinely hard to reach through the Task, or a numerical
  kernel needing an exact reference check (e.g. `test_batched_ldsc_rg.py`).
- **If a Task cannot be tested at its own level, that is a signal about the Task's interface, not
  a reason to drop down a level.** Refactoring its dependencies to make it testable (e.g. taking
  plain dataframes instead of a bespoke directory layout or an .xlsx sheet) is usually worthwhile
  on its own — the simplification is the point, and testability follows.

## Every assertion needs a cost/benefit case

Each assertion has a cost (friction when refactoring) and a benefit (bugs caught); it earns its
place only when the benefit exceeds the cost. Before adding an assertion, name the bug it would
catch — if the honest answer is "someone reworded a string," drop it.

- **Don't assert on error-message text or log wording.** When a function raises, assert *that* it
  raises, not what the message says:
  ```python
  # no
  with pytest.raises(RemapRootUnavailableError) as raised:
      check_remap_roots_available((rule,))
  assert str(missing_root) in str(raised.value)

  # yes
  with pytest.raises(RemapRootUnavailableError):
      check_remap_roots_available((rule,))
  ```
  Same for logs: assert that a warning was emitted, not the words in it. If a test's *only* content
  is message composition, delete it. Prefer asserting on the return value, the raised type, the
  state left on disk, or a call count.
- **Don't duplicate an error-message substring in `pytest.raises(match=...)`.** It copies the
  wording across two files. Prefer omitting `match=` entirely; only when a test must tell two
  same-type errors apart, define a module-level string constant in the *source* module, use it in
  the f-string assert message, and import it into the test for `match=`.

## Avoid brittle / cargo-cult tests

Don't add tests that mechanically follow TDD form but catch no real bug, and don't write
assertions that break under behavior-preserving refactors.

- Don't unit-test plain constants, or construct a struct and immediately assert a field you just
  set. Trust the type system / fail-fast invariants.
- Avoid fragile count assertions (`sum("gctb" in c for c in cmds) == 3`) that a non-bug change (an
  extra call) would break; assert the essential invariant instead (`any("--gwfm RC" in c ...)`).
- When a literal appears in both the code-under-test and the assertion, define ONE shared constant
  used in both places. Define local constants at the top of a test rather than retyping the same
  number in setup and assertion.
- Prefer `isinstance(asset, DirectoryAsset)` over re-deriving the same fact
  (`Path(asset.path).is_dir()`).
- **Asset-generator / wiring builders don't warrant dedicated structural tests.** Functions in
  `mecfs_bio/asset_generator/` (`generate_assets_*`, `build_*_groups`) that just construct and wire
  Tasks together should not get tests that assert facts the builder set directly (group counts,
  distinct asset ids, deps membership) — those break on refactors and catch nothing. A
  demonstrator/asset module constructs at import time, so a wiring error surfaces on import; that
  is enough coverage. Reserve tests for Tasks with real `execute()` logic.

## Inject dependencies — never monkeypatch or mock

Never use pytest `monkeypatch`, `unittest.mock`, or any mocking library. Instead make the
dependency an injectable parameter with a sensible production default, so a test passes its own
implementation:

```python
def execute_command_with_retries(
    cmd, ...,
    executor: Callable[[list[str]], str] = execute_command,
    sleep: Callable[[float], None] = time.sleep,
): ...
```

The test passes its own `flaky_executor` and a list-appending `sleep` — no monkeypatch. Mirror the
existing `FakeDownloader` / `robust_download(downloader=...)` pattern in `mecfs_bio/util/download/`.
Functional injection keeps the seams explicit in the signature and avoids brittle patching of
module internals.

## Don't gate tests on library presence

Avoid `@pytest.mark.skipif(not _have_rpy2(...), ...)` or try/except import probes that skip. The
environment is fully managed by pixi, so declared dependencies are always present. Such guards are
unnecessary and let a *real* failure (a genuinely broken library, an erroring test) masquerade as a
harmless "skip" nobody notices. If some tests genuinely need a heavy/optional dependency (e.g.
R-comparison tests via rpy2), separate them into their own module (`test_*_r_comparison.py`) —
still no skipif.

## System tests for Docker-based tools

For Tasks that invoke external tools via Docker (MiXeR, LAVA, SBayesRC):

1. **Use hello-world/vignette data from the tool's repo** — download a small tarball via
   `DownloadFileTask` + `ExtractTarGzipTask`, then prepare reference data via a preprocessing Task.
2. **Docker mount gotcha** — the asset store (under `tmp_path`) is not under `$PWD`, so
   `-v $PWD:/home` can't see it. Mount the reference dir explicitly via `extra_mounts` at a
   container path like `/ref_data`; use relative paths (under `$PWD`) for output files.
3. **Add `PreformattedDataSource` variants** when test data is already in the tool's native format,
   to skip column conversion.
4. **Define the Task inside the test function** (not at module level) if its `__attrs_post_init__`
   invokes Docker — otherwise collection fails on machines without Docker.
5. **Use fast/demo args** (`--fit-sequence diffevo-fast`, `--chr2use 21-22`) to keep tests under
   ~30 seconds. Pattern follows `test_lava.py`.

See also [python-style.md](python-style.md) for the injection and typing conventions these build on.
