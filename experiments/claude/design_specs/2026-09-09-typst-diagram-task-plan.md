# Typst Diagram Task Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a build-system Task that compiles Typst source into an SVG diagram asset and plugs it into the existing figure system for use in MkDocs docs.

**Architecture:** A new leaf `GeneratingTask` (`TypstDiagramTask`) takes Typst source (either a committed `.typ` file path or an inline string, behind one self-validating `TypstSource` type), shells out to the `typst` CLI via the injected `execute_command` runner, and returns a `FileAsset` pointing at the produced `.svg`. A new lean `DiagramFileMeta` describes the asset; it is registered in the `Meta` union, the asset-store path map, and the figure exporter so diagrams flow through the existing generate/publish machinery unchanged.

**Tech Stack:** Python 3.13, attrs (frozen), pixi (conda-forge `typst`), pytest, structlog. Typst packages (CeTZ/Fletcher) are fetched from the Typst registry at compile time and are not needed for any test in this plan.

**Spec:** `experiments/claude/design_specs/2026-09-09-typst-diagram-task-design.md`

## Global Constraints

- Run all commands via pixi: `pixi r <cmd>`; Python scripts via `pixi r python <script>`.
- After each task, run `pixi r invoke green` (lintfix, format, spellcheck, link, import, typecheck, test). `pytest-testmon` skips unaffected tests.
- Docstrings: no backticks around inline code, no RST.
- Paths: use `Path`/`PurePath` in memory; strings only at the CLI/text boundary. Relative store paths are `PurePath`.
- Inject dependencies; never monkeypatch/mock. The subprocess runner is a param with a production default of `execute_command`.
- Do not assert on error-message/log wording in tests; do not pass `match=` with a duplicated string. Avoid brittle "assert a field I just set" tests.
- Do not add skipif library-presence guards; pixi guarantees `typst` is present.
- Commits: end each commit message with the repo's standard attribution trailer (Co-Authored-By + Claude-Session lines).
- Single source of truth for the diagram file extension: the module-level `DIAGRAM_EXTENSION = ".svg"` constant defined in Task 1; reference it everywhere, never a repeated `".svg"` literal.

---

### Task 1: `DiagramFileMeta` + registration in the Meta union and asset-store path map

**Files:**
- Create: `mecfs_bio/build_system/meta/diagram_file_meta.py`
- Modify: `mecfs_bio/build_system/meta/meta.py` (add to the `Meta` union)
- Modify: `mecfs_bio/build_system/rebuilder/metadata_to_path/simple_meta_to_path.py` (import + `isinstance` branch + `_DIAGRAMS` constant)
- Test: `test_mecfs_bio/unit/build_system/meta/test_diagram_file_meta.py`

**Interfaces:**
- Produces:
  - `DIAGRAM_EXTENSION: str = ".svg"` (module constant in `diagram_file_meta.py`)
  - `DiagramFileMeta(FileMeta)` — frozen attrs, field `id: AssetId` (converter `AssetId`); property `asset_id -> AssetId` returning `id`.
  - Asset-store relative path for a `DiagramFileMeta` is `PurePath("diagrams") / (asset_id + DIAGRAM_EXTENSION)`.

Note on design: `DiagramFileMeta` carries only `id` (no `trait`/`project`, no `sub_dir`). Asset ids are globally unique, and the figure destination is flat, so no per-diagram sub-directory is needed. Extension is fixed to `.svg` via the constant rather than a caller field.

- [ ] **Step 1: Write the failing test**

trivial test removed

- [ ] **Step 2: Run test to verify it fails**

skipped, because trivial test was removed

- [ ] **Step 3: Create `DiagramFileMeta`**

Create `mecfs_bio/build_system/meta/diagram_file_meta.py`:

```python
"""
Metadata describing a single SVG diagram asset generated from Typst source.
"""

from attrs import field, frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta

DIAGRAM_EXTENSION = ".svg"


@frozen
class DiagramFileMeta(FileMeta):
    """
    Metadata describing a single SVG diagram.

    A diagram is a conceptual illustration for the documentation, not tied to
    any GWAS, so unlike GWASPlotFileMeta it carries no trait or project. The
    file extension is always .svg.
    """

    id: AssetId = field(converter=AssetId)

    @property
    def asset_id(self) -> AssetId:
        return self.id
```

- [ ] **Step 4: Register in the `Meta` union**

Modify `mecfs_bio/build_system/meta/meta.py`: add the import
```python
from mecfs_bio.build_system.meta.diagram_file_meta import DiagramFileMeta
```
and add `| DiagramFileMeta` to the `Meta = ( ... )` union.

- [ ] **Step 5: Add the asset-store path branch**

Modify `mecfs_bio/build_system/rebuilder/metadata_to_path/simple_meta_to_path.py`:

Add the import near the other meta imports:
```python
from mecfs_bio.build_system.meta.diagram_file_meta import (
    DIAGRAM_EXTENSION,
    DiagramFileMeta,
)
```
Add a module-level constant alongside `_GWAS` etc.:
```python
_DIAGRAMS = PurePath("diagrams")
```
Add this branch inside `simple_meta_to_relative_path` (before the final `raise`):
```python
    if isinstance(m, DiagramFileMeta):
        return _DIAGRAMS / (m.asset_id + DIAGRAM_EXTENSION)
```

- [ ] **Step 6: Run test to verify it passes**

skipped because trivial test removed

- [ ] **Step 7: Commit**

```bash
git add mecfs_bio/build_system/meta/diagram_file_meta.py \
        mecfs_bio/build_system/meta/meta.py \
        mecfs_bio/build_system/rebuilder/metadata_to_path/simple_meta_to_path.py \
        test_mecfs_bio/unit/build_system/meta/test_diagram_file_meta.py
git commit -m "feat: add DiagramFileMeta and register it in Meta union and path map"
```

---

### Task 2: `TypstSource` — self-validating path-or-string source

**Files:**
- Create: `mecfs_bio/build_system/task/typst_source.py`
- Test: `test_mecfs_bio/unit/build_system/task/test_typst_source.py`

**Interfaces:**
- Produces: `TypstSource` — frozen attrs with `path: Path | None = None` and `code: str | None = None`.
  - `__attrs_post_init__` asserts exactly one of `path`/`code` is set, and (path variant) that the file exists.
  - `resolve_to_path(self, scratch_dir: Path) -> Path` returns the existing `.typ` path, or writes `code` to `scratch_dir / "diagram.typ"` and returns that.

- [ ] **Step 1: Write the failing test**

Create `test_mecfs_bio/unit/build_system/task/test_typst_source.py`:

```python
from pathlib import Path

import pytest

from mecfs_bio.build_system.task.typst_source import TypstSource


def test_rejects_both_path_and_code(tmp_path: Path):
    f = tmp_path / "d.typ"
    f.write_text("x")
    with pytest.raises(AssertionError):
        TypstSource(path=f, code="x")


def test_rejects_neither():
    with pytest.raises(AssertionError):
        TypstSource()


def test_rejects_missing_file(tmp_path: Path):
    with pytest.raises(AssertionError):
        TypstSource(path=tmp_path / "does_not_exist.typ")


def test_resolves_code_to_scratch_file(tmp_path: Path):
    src = TypstSource(code="hello typst")
    resolved = src.resolve_to_path(tmp_path)
    assert resolved.read_text() == "hello typst"


def test_returns_existing_path_unchanged(tmp_path: Path):
    f = tmp_path / "d.typ"
    f.write_text("content")
    src = TypstSource(path=f)
    assert src.resolve_to_path(tmp_path) == f
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/test_typst_source.py -v`
Expected: FAIL (ModuleNotFoundError: `typst_source`).

- [ ] **Step 3: Implement `TypstSource`**

Create `mecfs_bio/build_system/task/typst_source.py`:

```python
"""
A Typst diagram source: exactly one of a committed .typ file path or an inline
string of Typst code.
"""

from pathlib import Path

from attrs import frozen

_SOURCE_FILENAME = "diagram.typ"


@frozen
class TypstSource:
    """
    Holds exactly one of a path to a .typ file or a literal string of Typst
    code. Enforcing exactly-one at construction makes the invalid states (both
    set, neither set, missing file) unrepresentable.
    """

    path: Path | None = None
    code: str | None = None

    def __attrs_post_init__(self):
        assert (self.path is None) != (
            self.code is None
        ), "TypstSource requires exactly one of path or code"
        if self.path is not None:
            assert (
                self.path.is_file()
            ), f"Typst source file does not exist: {self.path}"

    def resolve_to_path(self, scratch_dir: Path) -> Path:
        """
        Return a concrete .typ path for the source. For the path variant this is
        the existing file; for the code variant the string is written into
        scratch_dir and that new path is returned.
        """
        if self.path is not None:
            return self.path
        assert self.code is not None  # guaranteed by __attrs_post_init__
        target = scratch_dir / _SOURCE_FILENAME
        target.write_text(self.code)
        return target
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/test_typst_source.py -v`
Expected: PASS (5 passed).

- [ ] **Step 5: Commit**

```bash
git add mecfs_bio/build_system/task/typst_source.py \
        test_mecfs_bio/unit/build_system/task/test_typst_source.py
git commit -m "feat: add self-validating TypstSource"
```

---

### Task 3: `typst` dependency + `TypstDiagramTask`

**Files:**
- Modify: `pyproject.toml` (add `typst` under `[tool.pixi.feature.analysis.dependencies]`)
- Create: `mecfs_bio/build_system/task/typst_diagram_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/test_typst_diagram_task.py`

**Interfaces:**
- Consumes: `DiagramFileMeta`, `DIAGRAM_EXTENSION` (Task 1); `TypstSource` (Task 2); `execute_command` (`mecfs_bio/util/subproc/run_command.py`, signature `(cmd: list[str]) -> str`); `FileAsset`, `Fetch`, `WF`, `GeneratingTask`/`Task`.
- Produces: `TypstDiagramTask(GeneratingTask)` — frozen attrs with `meta: DiagramFileMeta`, `source: TypstSource`, `executor: Callable[[list[str]], str] = execute_command`. `deps -> []`. `execute(scratch_dir, fetch, wf) -> FileAsset` compiling `typst compile <src.typ> <out.svg>`.

- [ ] **Step 1: Add the `typst` dependency**

Modify `pyproject.toml`, in the `[tool.pixi.feature.analysis.dependencies]` table, add:
```toml
typst = ">=0.15,<0.16"
```

- [ ] **Step 2: Install so the binary is available**

Run: `pixi install`
Then confirm the CLI is on PATH: `pixi r typst --version`
Expected: prints a `typst 0.15.x` version line.

- [ ] **Step 3: Write the failing test**

Create `test_mecfs_bio/unit/build_system/task/test_typst_diagram_task.py`:

```python
from pathlib import Path

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.diagram_file_meta import DiagramFileMeta
from mecfs_bio.build_system.task.typst_diagram_task import TypstDiagramTask
from mecfs_bio.build_system.task.typst_source import TypstSource
from mecfs_bio.build_system.wf.base_wf import make_wf

# A single-page Typst document with visible content, no @preview imports, so the
# render is hermetic (no package-registry fetch) and fast.
_SMOKE_TYPST = (
    "#set page(width: auto, height: auto, margin: 4pt)\n"
    "#circle(radius: 8pt)\n"
)


def _no_fetch(asset_id: AssetId) -> Asset:
    raise ValueError("TypstDiagramTask has no dependencies")


def test_renders_svg_from_inline_code(tmp_path: Path):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    task = TypstDiagramTask(
        meta=DiagramFileMeta(AssetId("smoke")),
        source=TypstSource(code=_SMOKE_TYPST),
    )
    result = task.execute(scratch_dir=scratch, fetch=_no_fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    head = result.path.read_bytes()[:256].lstrip()
    assert head.startswith(b"<svg") or head.startswith(b"<?xml")


def test_renders_svg_from_file_source(tmp_path: Path):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    typ_file = tmp_path / "diagram.typ"
    typ_file.write_text(_SMOKE_TYPST)
    task = TypstDiagramTask(
        meta=DiagramFileMeta(AssetId("smoke_file")),
        source=TypstSource(path=typ_file),
    )
    result = task.execute(scratch_dir=scratch, fetch=_no_fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    head = result.path.read_bytes()[:256].lstrip()
    assert head.startswith(b"<svg") or head.startswith(b"<?xml")
```

- [ ] **Step 4: Run test to verify it fails**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/test_typst_diagram_task.py -v`
Expected: FAIL (ModuleNotFoundError: `typst_diagram_task`).

- [ ] **Step 5: Implement `TypstDiagramTask`**

Create `mecfs_bio/build_system/task/typst_diagram_task.py`:

```python
"""
Task that compiles Typst source into an SVG diagram asset.
"""

from collections.abc import Callable
from pathlib import Path

from attrs import field, frozen

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.diagram_file_meta import (
    DIAGRAM_EXTENSION,
    DiagramFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.task.typst_source import TypstSource
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.subproc.run_command import execute_command


@frozen
class TypstDiagramTask(GeneratingTask):
    """
    Compiles a Typst source (a .typ file or an inline string) into an SVG
    diagram by shelling out to the typst CLI.
    """

    meta: DiagramFileMeta
    source: TypstSource
    executor: Callable[[list[str]], str] = field(default=execute_command)

    @property
    def deps(self) -> list[Task]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> FileAsset:
        src_path = self.source.resolve_to_path(scratch_dir)
        out_path = scratch_dir / (self.meta.asset_id + DIAGRAM_EXTENSION)
        self.executor(["typst", "compile", str(src_path), str(out_path)])
        return FileAsset(out_path)
```

- [ ] **Step 6: Run test to verify it passes**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/test_typst_diagram_task.py -v`
Expected: PASS (2 passed). If the SVG check fails, print the first bytes to see what typst emitted (some versions emit a leading XML declaration, both of which the assertion already accepts).

- [ ] **Step 7: Commit**

```bash
git add pyproject.toml pixi.lock \
        mecfs_bio/build_system/task/typst_diagram_task.py \
        test_mecfs_bio/unit/build_system/task/test_typst_diagram_task.py
git commit -m "feat: add TypstDiagramTask and typst dependency"
```

---

### Task 4: figure-system integration

**Files:**
- Modify: `mecfs_bio/figures/figure_exporter.py` (import; add to `ValidFigureMeta`; add branch + helper in `get_figure_destination`)
- Test: `test_mecfs_bio/unit/figures/test_diagram_figure_destination.py`

**Interfaces:**
- Consumes: `DiagramFileMeta`, `DIAGRAM_EXTENSION` (Task 1); existing `get_figure_destination(meta, fig_dir) -> Path`.
- Produces: `get_figure_destination` returns `fig_dir / (asset_id + DIAGRAM_EXTENSION)` for a `DiagramFileMeta`; `DiagramFileMeta` is a member of `ValidFigureMeta` so the exporter accepts diagram tasks and its `FileAsset` copy branch handles them.

- [ ] **Step 1: Write the failing test**

Create `test_mecfs_bio/unit/figures/test_diagram_figure_destination.py`:

```python
from pathlib import Path

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.diagram_file_meta import (
    DIAGRAM_EXTENSION,
    DiagramFileMeta,
)
from mecfs_bio.figures.figure_exporter import get_figure_destination


def test_diagram_lands_flat_in_fig_dir(tmp_path: Path):
    meta = DiagramFileMeta(AssetId("concept"))
    assert get_figure_destination(meta=meta, fig_dir=tmp_path) == tmp_path / (
        "concept" + DIAGRAM_EXTENSION
    )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi r pytest test_mecfs_bio/unit/figures/test_diagram_figure_destination.py -v`
Expected: FAIL — `get_figure_destination` raises `ValueError("Unknown meta type ...")` for `DiagramFileMeta`.

- [ ] **Step 3: Wire `DiagramFileMeta` into the exporter**

Modify `mecfs_bio/figures/figure_exporter.py`:

Add the import near the other meta imports:
```python
from mecfs_bio.build_system.meta.diagram_file_meta import (
    DIAGRAM_EXTENSION,
    DiagramFileMeta,
)
```
Add `DiagramFileMeta` to the `ValidFigureMeta` union:
```python
ValidFigureMeta = (
    GWASPlotFileMeta
    | DirectoryFigureMeta
    | MarkdownFileMeta
    | GWASLabManhattanQQPlotMeta
    | ResultTableMeta
    | DiagramFileMeta
)
```
Add a helper next to the other `get_*_fig_file_path` helpers:
```python
def get_diagram_fig_file_path(meta: DiagramFileMeta, fig_dir: Path) -> Path:
    return fig_dir / (meta.asset_id + DIAGRAM_EXTENSION)
```
Add a branch in `get_figure_destination` (before the final `raise`):
```python
    if isinstance(meta, DiagramFileMeta):
        return get_diagram_fig_file_path(meta=meta, fig_dir=fig_dir)
```
No change to `export` is needed: `DiagramFileMeta` is not a `DirectoryFigureMeta`, so the existing `FileAsset` copy branch handles it.

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi r pytest test_mecfs_bio/unit/figures/test_diagram_figure_destination.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add mecfs_bio/figures/figure_exporter.py \
        test_mecfs_bio/unit/figures/test_diagram_figure_destination.py
git commit -m "feat: accept DiagramFileMeta in the figure exporter"
```

---

### Task 5: `mecfs_bio/diagrams/` package + one example diagram wired into the figure list

**Files:**
- Create: `mecfs_bio/diagrams/__init__.py`
- Create: `mecfs_bio/diagrams/example_genotype_phenotype.typ`
- Create: `mecfs_bio/diagrams/example_genotype_phenotype.py`
- Modify: `mecfs_bio/figures/figure_task_list.py` (import the singleton; append to `ALL_FIGURE_TASKS`)

**Interfaces:**
- Consumes: `TypstDiagramTask` (Task 3), `TypstSource` (Task 2), `DiagramFileMeta` (Task 1).
- Produces: module-level singleton `EXAMPLE_GENOTYPE_PHENOTYPE_DIAGRAM: TypstDiagramTask`, referenced from `ALL_FIGURE_TASKS`.

This task establishes the source-package convention the spec calls for and proves the end-to-end path (source -> task -> figure list). Per project convention there is no structural unit test for a wiring-only builder; import-time construction plus `invoke green` is the check.

- [ ] **Step 1: Create the diagrams package**

Create `mecfs_bio/diagrams/__init__.py`:
```python
"""
Typst sources for conceptual documentation diagrams and the tasks that build
them into SVG figures.
"""
```

- [ ] **Step 2: Create the example Typst source**

Create `mecfs_bio/diagrams/example_genotype_phenotype.typ`:
```typst
#set page(width: auto, height: auto, margin: 6pt)
#set text(size: 11pt)

#let r = 14pt
#circle(radius: r, name: "g")[Genotype]
#place(dx: 90pt, dy: 0pt, circle(radius: r, name: "p")[Phenotype])
#line((2 * r, r), (90pt, r), mark: (end: ">"))
```

Note: this uses only bare Typst (no `@preview` imports), so the build does not touch the package registry. If you later add a CeTZ/Fletcher diagram, pin the exact version in the import line, e.g. `#import "@preview/cetz:0.5.2"`.

- [ ] **Step 3: Create the task singleton module**

Create `mecfs_bio/diagrams/example_genotype_phenotype.py`:
```python
"""
Example diagram: a genotype node pointing to a phenotype node. Establishes the
diagram-source convention; safe to replace with a real documentation diagram.
"""

from pathlib import Path

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.diagram_file_meta import DiagramFileMeta
from mecfs_bio.build_system.task.typst_diagram_task import TypstDiagramTask
from mecfs_bio.build_system.task.typst_source import TypstSource

_SOURCE_PATH = Path(__file__).parent / "example_genotype_phenotype.typ"

EXAMPLE_GENOTYPE_PHENOTYPE_DIAGRAM = TypstDiagramTask(
    meta=DiagramFileMeta(AssetId("example_genotype_phenotype")),
    source=TypstSource(path=_SOURCE_PATH),
)
```

- [ ] **Step 4: Wire it into the figure list**

Modify `mecfs_bio/figures/figure_task_list.py`:

Add the import with the other imports at the top:
```python
from mecfs_bio.diagrams.example_genotype_phenotype import (
    EXAMPLE_GENOTYPE_PHENOTYPE_DIAGRAM,
)
```
Add an entry to the `ALL_FIGURE_TASKS: list[Task]` list (a `# Diagrams` comment plus the singleton):
```python
    # Diagrams
    EXAMPLE_GENOTYPE_PHENOTYPE_DIAGRAM,
```

- [ ] **Step 5: Verify the wiring imports and constructs**

Run: `pixi r python -c "from mecfs_bio.figures.figure_task_list import ALL_FIGURE_TASKS; from mecfs_bio.diagrams.example_genotype_phenotype import EXAMPLE_GENOTYPE_PHENOTYPE_DIAGRAM as d; assert d in ALL_FIGURE_TASKS; print('wired', d.asset_id)"`
Expected: prints `wired example_genotype_phenotype`.

- [ ] **Step 6: Optionally verify a real end-to-end render**

Run: `pixi r python -c "from pathlib import Path; import tempfile; from mecfs_bio.diagrams.example_genotype_phenotype import EXAMPLE_GENOTYPE_PHENOTYPE_DIAGRAM as d; from mecfs_bio.build_system.wf.base_wf import make_wf; t=Path(tempfile.mkdtemp()); a=d.execute(scratch_dir=t, fetch=lambda i: (_ for _ in ()).throw(ValueError()), wf=make_wf()); print(a.path, a.path.stat().st_size)"`
Expected: prints a path ending in `example_genotype_phenotype.svg` and a nonzero size.

- [ ] **Step 7: Commit**

```bash
git add mecfs_bio/diagrams/ mecfs_bio/figures/figure_task_list.py
git commit -m "feat: add diagrams package with an example Typst diagram"
```

---

### Final verification

- [ ] **Step 1: Run the full green suite**

Run: `pixi r invoke green 2>&1 | tee /tmp/green_typst.log`
Expected: exits 0. Because testmon skips unaffected tests, confirm the new tests actually ran by grepping the log or re-running them explicitly:
`pixi r pytest test_mecfs_bio/unit/build_system/task/test_typst_diagram_task.py test_mecfs_bio/unit/build_system/task/test_typst_source.py test_mecfs_bio/unit/build_system/meta/test_diagram_file_meta.py test_mecfs_bio/unit/figures/test_diagram_figure_destination.py -v`
Expected: all pass.

- [ ] **Step 2: Confirm the branch is clean and push if desired**

Run: `git status` (expect clean) and `git log --oneline main..HEAD` (expect the five feature commits).
```

## Self-Review

**Spec coverage** — every spec section maps to a task:
- DiagramFileMeta → Task 1; Meta union + path map registration → Task 1 (the spec's Edit 2 and the meta-union edit).
- TypstDiagramTask + TypstSource + executor injection → Tasks 2 and 3.
- typst pixi dependency, version pin → Task 3.
- figure_exporter ValidFigureMeta + get_figure_destination (spec Edit 3) → Task 4.
- mecfs_bio/diagrams/ source layout + figure_task_list wiring → Task 5.
- Error handling (invalid/missing source fails at construction; compile failure surfaces via execute_command's CalledProcessError) → Task 2 (construction) and inherent to Task 3's `execute` (no swallow).
- Testing conventions (real binary, no skipif, inject runner, no message-text asserts, hermetic smoke snippet) → Tasks 2–4.
- Out-of-scope items (Graphviz/drawsvg/SVGO/TikZ/vendoring; manual cache invalidation) → correctly produce no tasks.

**Placeholder scan** — no TBD/TODO/"add error handling"/"similar to Task N"; all code is spelled out.

**Type consistency** — `DIAGRAM_EXTENSION` and `DiagramFileMeta.asset_id` used identically across Tasks 1, 3, 4; `TypstSource(path=, code=)` and `resolve_to_path(scratch_dir)` consistent between Tasks 2, 3, 5; `execute(scratch_dir, fetch, wf)` signature matches `GeneratingTask.execute` and the test call sites; `execute_command` signature `(list[str]) -> str` matches the injected default.
