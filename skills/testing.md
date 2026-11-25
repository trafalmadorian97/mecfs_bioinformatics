# Testing

Guiding principle: a test should fail when observable behavior changes, and only then. A test
that breaks under a behavior-preserving refactor costs effort without catching bugs.

## What to test

- **Test the Task, not its private helpers.** The Task is the public API of the build system.
  Construct it, call `execute` with `FakeTask` dependencies (`mecfs_bio/build_system/task/fake_task.py`),
  synthetic inputs written to `tmp_path`, and a `fetch` that maps asset ids to assets. Example:
  `test_mecfs_bio/unit/build_system/task/ppp_ldsc/test_ppp_protein_rg_task.py`.
- Test a helper directly only when its logic is hard to reach through the Task, or when it is a
  numerical kernel that needs an exact reference check.
- **If a Task is hard to test, simplify its interface.** For example, have it take plain
  dataframes instead of a bespoke directory layout. Don't drop down to testing internals.
- **Don't write structural tests for asset generators.** Code in `mecfs_bio/asset_generator/`
  that only constructs and wires Tasks runs at import time, so an import already catches wiring
  errors.

## What to assert

Every assertion costs refactoring friction, so it must catch a real bug. Before adding one, name
the bug it would catch.

- **Assert that an exception is raised, not what its message says.** Likewise, assert that a
  warning was logged, not its wording. Don't copy message text into `pytest.raises(match=...)`. If
  two errors of the same type must be told apart, define the message as a constant in the source
  module and import it into the test.
- **Don't assert what you just set**, and don't test plain constants.
- **Avoid fragile counts** (`assert len(calls) == 3`) that a harmless extra call would break.
  Assert the essential invariant instead (`assert any("--flag" in c for c in calls)`).
- When the code and the test both need the same literal, share one constant.

## How to write tests

- **Never use `monkeypatch`, `unittest.mock`, or any mocking library.** Make the dependency a
  parameter with a production default, and pass a fake from the test. Example:
  `robust_download` (`mecfs_bio/util/download/`) takes a `Downloader`, and its test passes a
  `FakeDownloader`.
- **Don't skip tests based on whether a library is installed** (`skipif` or try-import probes).
  pixi guarantees declared dependencies, and a skip guard hides real breakage. Put tests that need
  heavy dependencies in their own module instead.
- **Keep unit tests fast**, about a second each. Size synthetic data by the structure the test
  needs (blocks, row groups, seek points), not by realism; shrink tunable knobs such as chunk sizes
  instead of growing data. Build expensive fixtures once per module. If you're unsure how fast a
  test is, time it.

## Running tests

- `pixi r invoke green` is the full check (see [workflow.md](workflow.md)).
- **For a targeted run, list explicit paths under `test_mecfs_bio/unit/`.** Never pass the
  `test_mecfs_bio` root with a `-k` filter: it also matches slow system tests (Docker, R, real
  data). Add `--durations=5` so slow tests are visible.

## System tests for Docker-based tools

For Tasks that wrap a tool run through Docker (e.g. MiXeR, LAVA, SBayesRC); `test_mecfs_bio/system/test_lava.py`
is the model:

1. Use the tool's own hello-world or vignette data, downloaded with `DownloadFileTask` and
   `ExtractTarGzipTask`.
2. The asset store under `tmp_path` is outside `$PWD`, so a `-v $PWD:...` mount can't see it.
   Mount reference data explicitly with `extra_mounts`.
3. If the Task's `__attrs_post_init__` calls Docker, define it inside the test function so test
   collection doesn't need Docker.
4. Use the tool's fast or demo settings (e.g. a couple of chromosomes) to stay under about 30 s.
