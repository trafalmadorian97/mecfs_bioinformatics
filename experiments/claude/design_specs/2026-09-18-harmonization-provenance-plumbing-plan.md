# Harmonization Provenance Plumbing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Thread genome-build harmonization provenance through the metadata layer so that a FASTA-harmonized GWAS and the reference tables it is fine-mapped against both declare, in their metadata, which genome build they are oriented to.

**Architecture:** Add a `HarmonizationInfo` value type and carry it on `FilteredGWASDataMeta`; introduce two sibling reference-meta leaves (`FASTAMeta`, `HarmonizableReferenceTableMeta`) under the existing abstract `FileMeta`/`DirMeta`; add a capability-restricted `RenameColsTask` that can update that provenance accurately across a column rename. Tag the hg19/hg38 FASTA, the 1kg EUR panel, and the Broad UKBB LD labels; have `GenomeReferenceHarmonizationTask.create` stamp its output; propagate through the one orientation-preserving intermediate task (`JoinDataFramesTask`). This plan changes no build results and enforces nothing at runtime — it is pure plumbing that Plan 2 consumes.

**Tech Stack:** Python 3.13, attrs (`@frozen`), polars, the repo's `Task`/`Asset`/`Meta` build system, pixi, pytest (+ pytest-testmon).

**Spec:** `experiments/claude/design_specs/susie_polyfun_unordered_allele_key_cleanup.md`

## Global Constraints

- All commands run via pixi: `pixi r <cmd>`, `pixi r python <script>`, `pixi r invoke green`.
- After any significant change run `pixi r invoke green`. testmon skips unaffected tests, so exit 0 alone does not prove a test ran — name the test in `pixi r pytest ... -v` when you need to see it execute.
- `GenomeBuild = Literal["19", "38"]` from `mecfs_bio.constants.genomic_coordinate_constants`.
- Column-name constants from `mecfs_bio.constants.gwaslab_constants` (`GWASLAB_CHROM_COL="CHR"`, `GWASLAB_POS_COL="POS"`, `GWASLAB_EFFECT_ALLELE_COL="EA"`, `GWASLAB_NON_EFFECT_ALLELE_COL="NEA"`, `GWASLAB_RSID_COL`). Never repeat the string literals in code.
- Prefer helper free functions over methods; Path over str; assert-narrow over cast; no monkeypatching (inject dependencies). The test task stub is `FakeTask(meta=...)` from `mecfs_bio.build_system.task.fake_task` (its `execute` raises, so use it only for `create()`/meta assertions).
- **All imports at the top of the file** — never inside a function or test body (the sole exception in this repo is `tasks.py`).
- Docstrings: no backticks around inline code, no RST.
- Meta classes are `@frozen` attrs. A required field after defaulted parent fields needs `kw_only=True` to avoid the "non-default follows default" error.

### Testing philosophy for this plan
Write a test only where it asserts real logic with a real failure mode. Do **not** test:
- that an asset generator wired A into B (tautological — asset generators are wiring);
- that a `create()` stamped a constant we just passed in (field-you-just-set);
- thin wrappers over library calls (e.g. a polars `.rename`).
Only two tasks below carry tests: the path-resolution branches (Task 2 — a wrong path silently
orphans assets) and `RenameColsTask`'s provenance remap + collision guard (Task 3 — the one piece
of novel computation). Everything else is verified by `pixi r invoke green` and, at runtime, by
Plan 2's SUSIE build-match assertion (which fails loudly at graph construction if provenance does
not reach SUSIE).

---

### Task 1: `HarmonizationInfo` value type + `FilteredGWASDataMeta` field

**Files:**
- Create: `mecfs_bio/build_system/meta/harmonization_info.py`
- Modify: `mecfs_bio/build_system/meta/filtered_gwas_data_meta.py`
- Test: none (a data holder; a field-you-just-set test is cargo-cult).

**Interfaces:**
- Produces: `HarmonizationInfo(build: GenomeBuild, ref_allele_col: str, pos_col: str)` (frozen); `FilteredGWASDataMeta.harmonization_info: HarmonizationInfo | None` (default `None`).

- [ ] **Step 1: Create the value type**

```python
# mecfs_bio/build_system/meta/harmonization_info.py
"""Provenance of genome-reference harmonization: which build an asset is oriented to,
and which columns carry that orientation."""

from attrs import frozen

from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild


@frozen
class HarmonizationInfo:
    """Records that, at the positions in pos_col, ref_allele_col equals the reference
    base of the given genome build. build is the load-bearing field for construction-time
    build-match checks; the column names document the claim and support an optional
    data-level verification."""

    build: GenomeBuild
    ref_allele_col: str
    pos_col: str
```

- [ ] **Step 2: Add the field to `FilteredGWASDataMeta`**

Add the import at the top and a trailing optional field (keep existing fields/order):

```python
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
# ...
@frozen
class FilteredGWASDataMeta(FileMeta):
    # (existing fields unchanged: id, trait, project, sub_dir, read_spec, extension)
    harmonization_info: HarmonizationInfo | None = None
```

- [ ] **Step 3: Verify it imports**

Run: `pixi r python -c "from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta; print(FilteredGWASDataMeta.__attrs_attrs__[-1].name)"`
Expected: prints `harmonization_info`.

- [ ] **Step 4: Green + commit**

```bash
pixi r invoke green
git add mecfs_bio/build_system/meta/harmonization_info.py mecfs_bio/build_system/meta/filtered_gwas_data_meta.py
git commit -m "feat: add HarmonizationInfo and FilteredGWASDataMeta.harmonization_info"
```

---

### Task 2: sibling reference-meta types + path-resolver branches

**Files:**
- Create: `mecfs_bio/build_system/meta/reference_meta/fasta_meta.py`
- Create: `mecfs_bio/build_system/meta/reference_meta/harmonizable_reference_table_meta.py`
- Modify: `mecfs_bio/build_system/rebuilder/metadata_to_path/simple_meta_to_path.py`
- Test: `test_mecfs_bio/unit/build_system/metadata_to_path/test_simple_meta_to_path_reference_provenance.py`

**Interfaces:**
- Produces: `FASTAMeta(group, sub_group, sub_folder, build: GenomeBuild, id, dirname=None)` (DirMeta sibling of `ReferenceDataDirectoryMeta`); `HarmonizableReferenceTableMeta(group, sub_group, sub_folder, extension, id, filename=None, read_spec=None, harmonization_info: HarmonizationInfo | None = None)` (FileMeta sibling of `ReferenceFileMeta`). Both resolve to the same store-relative path as their concrete sibling. This test guards a real failure mode: a wrong path silently orphans built assets and forces a rebuild.

- [ ] **Step 1: Write the failing path test**

```python
# test_mecfs_bio/unit/build_system/metadata_to_path/test_simple_meta_to_path_reference_provenance.py
from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.reference_meta.fasta_meta import FASTAMeta
from mecfs_bio.build_system.meta.reference_meta.harmonizable_reference_table_meta import (
    HarmonizableReferenceTableMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import ReferenceFileMeta
from mecfs_bio.build_system.rebuilder.metadata_to_path.simple_meta_to_path import (
    simple_meta_to_relative_path,
)


def test_fasta_meta_resolves_like_reference_directory():
    fasta = FASTAMeta(
        group="genome_sequence", sub_group="ucsc_hg19", sub_folder=PurePath("processed"),
        id=AssetId("idx"), build="19",
    )
    twin = ReferenceDataDirectoryMeta(
        group="genome_sequence", sub_group="ucsc_hg19", sub_folder=PurePath("processed"),
        id=AssetId("idx"),
    )
    assert simple_meta_to_relative_path(fasta) == simple_meta_to_relative_path(twin)


def test_harmonizable_table_resolves_like_reference_file():
    table = HarmonizableReferenceTableMeta(
        group="g", sub_group="s", sub_folder=PurePath("processed"), extension=".parquet",
        id=AssetId("t"), filename="f",
        harmonization_info=HarmonizationInfo(build="19", ref_allele_col="REF", pos_col="POS"),
    )
    twin = ReferenceFileMeta(
        group="g", sub_group="s", sub_folder=PurePath("processed"), extension=".parquet",
        id=AssetId("t"), filename="f",
    )
    assert simple_meta_to_relative_path(table) == simple_meta_to_relative_path(twin)
```

- [ ] **Step 2: Run it to confirm it fails**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/metadata_to_path/test_simple_meta_to_path_reference_provenance.py -v`
Expected: FAIL with ImportError (FASTAMeta not defined).

- [ ] **Step 3: Create `FASTAMeta`**

```python
# mecfs_bio/build_system/meta/reference_meta/fasta_meta.py
"""Directory metadata for an indexed genome FASTA, tagged with its genome build. A sibling
of ReferenceDataDirectoryMeta (not a subclass), so the concrete reference metas stay leaves."""

from pathlib import PurePath

from attrs import field, frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import DirMeta
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild


@frozen
class FASTAMeta(DirMeta):
    group: str
    sub_group: str
    sub_folder: PurePath
    build: GenomeBuild
    id: AssetId = field(converter=AssetId)
    dirname: str | None = None

    @property
    def asset_id(self) -> AssetId:
        return self.id
```

- [ ] **Step 4: Create `HarmonizableReferenceTableMeta`**

```python
# mecfs_bio/build_system/meta/reference_meta/harmonizable_reference_table_meta.py
"""File metadata for a tabular reference oriented to a genome build (LD panel labels, the
1kg allele-frequency panel). A sibling of ReferenceFileMeta (not a subclass): generic
transformer tasks that isinstance-check ReferenceFileMeta will not capture it, so piping one
through a generic transformer fails closed. Use RenameColsTask for renames."""

from pathlib import PurePath

from attrs import field, frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.read_spec.read_spec import ReadSpec


@frozen
class HarmonizableReferenceTableMeta(FileMeta):
    group: str
    sub_group: str
    sub_folder: PurePath
    extension: str
    id: AssetId = field(converter=AssetId)
    filename: str | None = None
    read_spec: ReadSpec | None = None
    harmonization_info: HarmonizationInfo | None = None

    def __attrs_post_init__(self):
        assert self.extension.startswith(".") or self.extension == ""

    @property
    def asset_id(self) -> AssetId:
        return self.id
```

- [ ] **Step 5: Add two branches to `simple_meta_to_relative_path`**

Import both new types at the top of `simple_meta_to_path.py` and add branches reusing the existing
file/dir path logic (place each near its family; siblings, so exact order is not required):

```python
    if isinstance(m, HarmonizableReferenceTableMeta):
        pth = _REFERENCE_DATA / m.group / m.sub_group / m.sub_folder
        if m.filename is not None:
            pth = pth / (m.filename + m.extension)
        else:
            pth = pth / str(m.id + m.extension)
        return pth

    if isinstance(m, FASTAMeta):
        dirname = m.dirname if m.dirname is not None else m.id
        return _REFERENCE_DATA / m.group / m.sub_group / m.sub_folder / dirname
```

- [ ] **Step 6: Run the test to verify it passes**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/metadata_to_path/test_simple_meta_to_path_reference_provenance.py -v`
Expected: PASS (both tests).

- [ ] **Step 7: Green + commit**

```bash
pixi r invoke green
git add mecfs_bio/build_system/meta/reference_meta/fasta_meta.py mecfs_bio/build_system/meta/reference_meta/harmonizable_reference_table_meta.py mecfs_bio/build_system/rebuilder/metadata_to_path/simple_meta_to_path.py test_mecfs_bio/unit/build_system/metadata_to_path/test_simple_meta_to_path_reference_provenance.py
git commit -m "feat: add FASTAMeta and HarmonizableReferenceTableMeta sibling reference metas"
```

---

### Task 3: `RenameColsTask` (rename that preserves provenance)

**Files:**
- Create: `mecfs_bio/build_system/task/rename_cols_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/test_rename_cols_task.py`

**Interfaces:**
- Consumes: `HarmonizableReferenceTableMeta`, `FilteredGWASDataMeta`, `HarmonizationInfo`, `ReferenceFileMeta` (fallback), `ParquetOutFormat`, `scan_dataframe_asset`, `FakeTask`.
- Produces: `RenameColsTask.create(source_task: Task, asset_id: str, renames: Mapping[str, str], out_format: OutFormat = ParquetOutFormat(), backend: ValidBackend = "polars") -> RenameColsTask`. When the source meta carries `harmonization_info`, the output meta carries it with `ref_allele_col`/`pos_col` remapped through `renames` (build unchanged), and asserts the two key columns do not collide. The tested logic is `_remap_harmonization_info` and the collision guard (real computation, not wiring).

- [ ] **Step 1: Write the failing tests**

```python
# test_mecfs_bio/unit/build_system/task/test_rename_cols_task.py
from pathlib import PurePath

import pytest

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.harmonizable_reference_table_meta import (
    HarmonizableReferenceTableMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.rename_cols_task import RenameColsTask


def _harmonizable_source(task_id: str) -> FakeTask:
    return FakeTask(
        meta=HarmonizableReferenceTableMeta(
            group="ukbb_reference_ld", sub_group="chr1_1_2", sub_folder=PurePath("raw"),
            extension=".gz", id=AssetId(task_id), filename="chr1_1_2",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            harmonization_info=HarmonizationInfo(
                build="19", ref_allele_col="allele1", pos_col="position"
            ),
        )
    )


def test_create_remaps_harmonization_info():
    task = RenameColsTask.create(
        source_task=_harmonizable_source("src"),
        asset_id="renamed",
        renames={"allele1": "NEA", "allele2": "EA", "position": "POS"},
    )
    assert isinstance(task.meta, HarmonizableReferenceTableMeta)
    info = task.meta.harmonization_info
    assert info is not None
    assert (info.build, info.ref_allele_col, info.pos_col) == ("19", "NEA", "POS")


def test_create_rejects_colliding_provenance_columns():
    with pytest.raises(AssertionError):
        RenameColsTask.create(
            source_task=_harmonizable_source("src"),
            asset_id="renamed",
            renames={"allele1": "position"},  # ref-allele col collides onto pos col name
        )
```

- [ ] **Step 2: Run to confirm failure**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/test_rename_cols_task.py -v`
Expected: FAIL with ImportError (RenameColsTask not defined).

- [ ] **Step 3: Implement `RenameColsTask`**

```python
# mecfs_bio/build_system/task/rename_cols_task.py
"""Rename dataframe columns. A deliberately narrow alternative to PipeDataFrameTask: because the
rename map is the whole operation, create() can update HarmonizationInfo (ref/pos columns)
accurately across the rename, which a general pipe cannot."""

from collections.abc import Mapping
from pathlib import Path, PurePath

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.reference_meta.harmonizable_reference_table_meta import (
    HarmonizableReferenceTableMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import ReferenceFileMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import OutFormat, ParquetOutFormat
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ValidBackend,
    get_extension_and_read_spec_from_format,
)
from mecfs_bio.build_system.wf.base_wf import WF

_OUT_FILENAME = "renamed.parquet"


def apply_renames(frame: pl.LazyFrame, renames: Mapping[str, str]) -> pl.LazyFrame:
    return frame.rename(dict(renames))


def remap_harmonization_info(
    info: HarmonizationInfo, renames: Mapping[str, str]
) -> HarmonizationInfo:
    new_ref = renames.get(info.ref_allele_col, info.ref_allele_col)
    new_pos = renames.get(info.pos_col, info.pos_col)
    assert new_ref != new_pos, (
        f"rename collides the ref-allele column and position column onto {new_ref!r}"
    )
    return HarmonizationInfo(build=info.build, ref_allele_col=new_ref, pos_col=new_pos)


@frozen
class RenameColsTask(Task):
    meta: Meta
    source_task: Task
    renames: Mapping[str, str]
    out_format: OutFormat = ParquetOutFormat()
    backend: ValidBackend = "polars"

    @property
    def deps(self) -> list[Task]:
        return [self.source_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source = fetch(self.source_task.asset_id)
        frame = apply_renames(
            scan_dataframe_asset(source, self.source_task.meta).to_polars(), self.renames
        )
        out_path = scratch_dir / _OUT_FILENAME
        frame.sink_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        source_task: Task,
        asset_id: str,
        renames: Mapping[str, str],
        out_format: OutFormat = ParquetOutFormat(),
        backend: ValidBackend = "polars",
    ) -> "RenameColsTask":
        source_meta = source_task.meta
        extension, read_spec = get_extension_and_read_spec_from_format(out_format=out_format)
        meta: Meta
        if isinstance(source_meta, HarmonizableReferenceTableMeta):
            meta = HarmonizableReferenceTableMeta(
                group=source_meta.group, sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"), extension=extension,
                id=AssetId(asset_id), filename=None, read_spec=read_spec,
                harmonization_info=(
                    remap_harmonization_info(source_meta.harmonization_info, renames)
                    if source_meta.harmonization_info is not None else None
                ),
            )
        elif isinstance(source_meta, FilteredGWASDataMeta):
            meta = FilteredGWASDataMeta(
                id=AssetId(asset_id), trait=source_meta.trait, project=source_meta.project,
                sub_dir=source_meta.sub_dir, read_spec=read_spec,
                harmonization_info=(
                    remap_harmonization_info(source_meta.harmonization_info, renames)
                    if source_meta.harmonization_info is not None else None
                ),
            )
        elif isinstance(source_meta, ReferenceFileMeta):
            meta = ReferenceFileMeta(
                group=source_meta.group, sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"), extension=extension,
                id=AssetId(asset_id), read_spec=read_spec,
            )
        else:
            raise ValueError(f"RenameColsTask: unsupported source meta {type(source_meta).__name__}")
        return cls(
            meta=meta, source_task=source_task, renames=dict(renames),
            out_format=out_format, backend=backend,
        )
```

Note: confirm `ValidBackend` and `get_extension_and_read_spec_from_format` are importable from `pipe_dataframe_task` (they are used there); if they live in a sibling module, import from there. If `scan_dataframe_asset(...).to_polars()` is not the local idiom for an eager polars frame, match how `harmonize_gwas_with_reference_table_via_chrom_pos_alleles.py` reads a df asset.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/test_rename_cols_task.py -v`
Expected: PASS (both).

- [ ] **Step 5: Green + commit**

```bash
pixi r invoke green
git add mecfs_bio/build_system/task/rename_cols_task.py test_mecfs_bio/unit/build_system/task/test_rename_cols_task.py
git commit -m "feat: add RenameColsTask that preserves harmonization provenance"
```

---

### Task 4: tag the FASTA assets with `FASTAMeta` + build

**Files:**
- Modify: `mecfs_bio/build_system/task/genome_reference_harmonization/indexed_fasta_task.py`
- Modify: `mecfs_bio/assets/reference_data/genome_sequence/ucsc_hg19_fasta.py`
- Modify: `mecfs_bio/assets/reference_data/genome_sequence/ucsc_hg38_fasta.py`
- Test: none (asset wiring; verified by green).

**Interfaces:**
- Produces: `IndexedFastaTask.create(fasta_gz_task, asset_id, build: GenomeBuild)` returns a task whose `.meta` is `FASTAMeta` carrying `build`. `DiscardDepsWrapper.meta` is `attrs.evolve(inner.meta, id=...)`, so the wrapped asset's meta stays `FASTAMeta`.

- [ ] **Step 1: Update `IndexedFastaTask.create` to emit `FASTAMeta`**

Change the `meta` field annotation to `FASTAMeta`, import `FASTAMeta` and `GenomeBuild` at the top, drop the now-unused `ReferenceDataDirectoryMeta` import if nothing else uses it, and update `create`:

```python
    meta: FASTAMeta
    # ...
    @classmethod
    def create(cls, fasta_gz_task: Task, asset_id: str, build: GenomeBuild) -> "IndexedFastaTask":
        source_meta = fasta_gz_task.meta
        assert isinstance(source_meta, ReferenceFileMeta), (
            f"expected a ReferenceFileMeta source for {asset_id}, got {type(source_meta).__name__}"
        )
        return cls(
            meta=FASTAMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                build=build,
            ),
            fasta_gz_task=fasta_gz_task,
        )
```

- [ ] **Step 2: Pass `build` at the two asset definitions**

In `ucsc_hg19_fasta.py`: `IndexedFastaTask.create(fasta_gz_task=UCSC_HG19_FASTA_GZ, asset_id="ucsc_hg19_indexed_fasta", build="19")`.
In `ucsc_hg38_fasta.py`: read the file first for its exact variable names, then add `build="38"` to the analogous `IndexedFastaTask.create(...)` call.

- [ ] **Step 3: Green + commit**

Run: `pixi r invoke green` (catches any other consumer of the old FASTA directory meta type).

```bash
git add mecfs_bio/build_system/task/genome_reference_harmonization/indexed_fasta_task.py mecfs_bio/assets/reference_data/genome_sequence/ucsc_hg19_fasta.py mecfs_bio/assets/reference_data/genome_sequence/ucsc_hg38_fasta.py
git commit -m "feat: tag indexed FASTA assets with FASTAMeta build"
```

---

### Task 5: tag the 1kg EUR panel with `HarmonizableReferenceTableMeta` + build

**Files:**
- Modify: `mecfs_bio/build_system/task/genome_reference_harmonization/reference_panel_task.py`
- Modify: `mecfs_bio/assets/reference_data/thousand_genomes/eur_panel_allele_frequencies.py`
- Test: none (asset wiring; verified by green).

**Interfaces:**
- Produces: `ReferencePanelAlleleFrequencyTask.create(vcf_task, asset_id, build: GenomeBuild)` returns a task whose `.meta` is `HarmonizableReferenceTableMeta` with `harmonization_info=HarmonizationInfo(build, ref_allele_col=PANEL_REF_COL, pos_col=GWASLAB_POS_COL)`. `PANEL_REF_COL="REF"` is already defined in `reference_panel_task.py`.

- [ ] **Step 1: Update `ReferencePanelAlleleFrequencyTask.create`**

Change the task's `meta:` field annotation to `HarmonizableReferenceTableMeta`, add a `build: GenomeBuild` param, import `HarmonizableReferenceTableMeta`, `HarmonizationInfo`, and `GenomeBuild` at the top (`PANEL_REF_COL` and `GWASLAB_POS_COL` are already in the module), and emit:

```python
    @classmethod
    def create(cls, vcf_task: Task, asset_id: str, build: GenomeBuild) -> "ReferencePanelAlleleFrequencyTask":
        source_meta = vcf_task.meta
        assert isinstance(source_meta, ReferenceFileMeta), (
            f"expected a ReferenceFileMeta source for {asset_id}, got {type(source_meta).__name__}"
        )
        return cls(
            meta=HarmonizableReferenceTableMeta(
                group="reference_panel_allele_frequencies",
                sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                filename="panel_allele_frequencies",
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                harmonization_info=HarmonizationInfo(
                    build=build, ref_allele_col=PANEL_REF_COL, pos_col=GWASLAB_POS_COL
                ),
            ),
            vcf_task=vcf_task,
        )
```

- [ ] **Step 2: Pass `build` at the asset definitions**

In `eur_panel_allele_frequencies.py`: add `build="19"` to the HG19 call and `build="38"` to the HG38 call.

- [ ] **Step 3: Green + commit**

Run: `pixi r invoke green`.

```bash
git add mecfs_bio/build_system/task/genome_reference_harmonization/reference_panel_task.py mecfs_bio/assets/reference_data/thousand_genomes/eur_panel_allele_frequencies.py
git commit -m "feat: tag 1kg EUR panel with HarmonizableReferenceTableMeta build"
```

---

### Task 6: `GenomeReferenceHarmonizationTask.create` checks build match + stamps output

**Files:**
- Modify: `mecfs_bio/build_system/task/genome_reference_harmonization/genome_reference_harmonization_task.py:364-391`
- Test: none as a standalone (the stamp is a field-you-just-set; the cross-build guard is exercised by `invoke green` — existing harmonization unit tests construct the task, and any that pass an untagged fasta/panel will now fail and must be updated to tagged assets or `FASTAMeta`/`HarmonizableReferenceTableMeta` values, which is the guard doing its job).

**Interfaces:**
- Consumes: `FASTAMeta` (`fasta_task.meta`), `HarmonizableReferenceTableMeta` (`panel_task.meta`), `HarmonizationInfo`, `GWASLAB_NON_EFFECT_ALLELE_COL`, `GWASLAB_POS_COL`, `GenomeBuild`.
- Produces: output `FilteredGWASDataMeta.harmonization_info = HarmonizationInfo(build, ref_allele_col="NEA", pos_col="POS")` where `build == fasta.build == panel.harmonization_info.build`; a mismatch raises at construction.

- [ ] **Step 1: Add the build-match helper + stamping**

Add a free helper (import `FASTAMeta`, `HarmonizableReferenceTableMeta`, `HarmonizationInfo`, `GWASLAB_NON_EFFECT_ALLELE_COL`, `GWASLAB_POS_COL`, `GenomeBuild` at the top):

```python
def resolve_harmonized_build(fasta_task: Task, panel_task: Task) -> GenomeBuild:
    fasta_meta = fasta_task.meta
    panel_meta = panel_task.meta
    assert isinstance(fasta_meta, FASTAMeta), (
        f"fasta_task must carry FASTAMeta, got {type(fasta_meta).__name__}"
    )
    assert isinstance(panel_meta, HarmonizableReferenceTableMeta) and panel_meta.harmonization_info is not None, (
        "panel_task must carry HarmonizableReferenceTableMeta with harmonization_info"
    )
    assert fasta_meta.build == panel_meta.harmonization_info.build, (
        f"fasta build {fasta_meta.build} != panel build {panel_meta.harmonization_info.build}"
    )
    return fasta_meta.build
```

Then in `create`, compute `build = resolve_harmonized_build(fasta_task, panel_task)` and add to the returned `FilteredGWASDataMeta(...)`:

```python
                harmonization_info=HarmonizationInfo(
                    build=build,
                    ref_allele_col=GWASLAB_NON_EFFECT_ALLELE_COL,
                    pos_col=GWASLAB_POS_COL,
                ),
```

- [ ] **Step 2: Green + commit**

Run: `pixi r invoke green`. If existing genome-reference-harmonization unit tests construct the task with a bare/untagged fasta or panel meta, update those call sites to use the tagged assets or `FASTAMeta`/`HarmonizableReferenceTableMeta` values — the new assertion is meant to reject untagged inputs.

```bash
git add mecfs_bio/build_system/task/genome_reference_harmonization/genome_reference_harmonization_task.py
git commit -m "feat: GenomeReferenceHarmonizationTask stamps harmonization_info and checks build match"
```

---

### Task 7: propagate `harmonization_info` through `JoinDataFramesTask`

**Files:**
- Modify: `mecfs_bio/build_system/task/join_dataframes_task.py` (the `FilteredGWASDataMeta` branch of `create_from_result_df`, ~lines 139-146)
- Test: none (a single field copy on an existing branch; the real check is Plan 2's runtime SUSIE build-match assertion, which fails at construction if this drops).

**Interfaces:**
- Produces: the `FilteredGWASDataMeta` branch copies `harmonization_info=source_meta.harmonization_info` (a join adds columns and subsets rows; it never flips alleles, so the orientation claim is preserved).

- [ ] **Step 1: Copy the field in the `FilteredGWASDataMeta` branch**

```python
        elif isinstance(source_meta, FilteredGWASDataMeta):
            meta = FilteredGWASDataMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                extension=extension,
                read_spec=read_spec,
                sub_dir=source_meta.sub_dir,
                harmonization_info=source_meta.harmonization_info,
            )
```

- [ ] **Step 2: Green + commit**

Run: `pixi r invoke green`.

```bash
git add mecfs_bio/build_system/task/join_dataframes_task.py
git commit -m "feat: JoinDataFramesTask preserves harmonization_info (orientation-preserving join)"
```

---

### Task 8: tag the Broad UKBB LD labels and rename them via `RenameColsTask`

**Files:**
- Modify: `mecfs_bio/asset_generator/ukbb_broad_ld_matrix_generator.py` (`get_genomic_interval_ld_labels_task`)
- Modify: `mecfs_bio/asset_generator/fine_mapping_asset_generator.py` (`ld_labels_task_renamed`)
- Modify: `mecfs_bio/asset_generator/polyfun_explain_fine_mapping_asset_generator.py` (`ld_labels_task_renamed` in `_build_shared_locus_inputs`)
- Test: none (asset-generator wiring). A manual end-to-end check is in Step 4.

**Interfaces:**
- Produces: the raw label download task's `.meta` is `HarmonizableReferenceTableMeta` with `harmonization_info=HarmonizationInfo(build="19", ref_allele_col="allele1", pos_col="position")` (verified: allele1 == REF for the Broad panel); the renamed labels task is a `RenameColsTask` carrying `HarmonizationInfo(build="19", ref_allele_col="NEA", pos_col="POS")`.

- [ ] **Step 1: Tag the raw label download meta**

In `get_genomic_interval_ld_labels_task`, replace the `ReferenceFileMeta(...)` on the `DownloadFileTask` with `HarmonizableReferenceTableMeta(...)` — same `group`/`sub_group`/`sub_folder`/`id`/`filename`/`extension`/`read_spec` — plus (import `HarmonizableReferenceTableMeta` and `HarmonizationInfo` at the top):

```python
            harmonization_info=HarmonizationInfo(
                build="19", ref_allele_col="allele1", pos_col="position"
            ),
```

- [ ] **Step 2: Swap both generators' label rename to `RenameColsTask`**

In `fine_mapping_asset_generator.py`, replace the `PipeDataFrameTask.create(... RenameColPipe ...)` that builds `ld_labels_task_renamed` with (import `RenameColsTask` at the top):

```python
    ld_labels_task_renamed = RenameColsTask.create(
        source_task=ld_labels_task,
        asset_id=ld_labels_task.asset_id + "_renamed",
        renames={
            "rsid": GWASLAB_RSID_COL,
            "chromosome": GWASLAB_CHROM_COL,
            "position": GWASLAB_POS_COL,
            "allele1": GWASLAB_NON_EFFECT_ALLELE_COL,
            "allele2": GWASLAB_EFFECT_ALLELE_COL,
        },
    )
```

Apply the identical replacement in `polyfun_explain_fine_mapping_asset_generator.py`. Remove `PipeDataFrameTask`/`RenameColPipe`/`ParquetOutFormat` imports from each file if no longer referenced.

- [ ] **Step 3: Green**

Run: `pixi r invoke green`. If a generic transformer now raises on the label task (a sibling meta is not a `ReferenceFileMeta`), the only intended consumer is the rename just converted to `RenameColsTask`; any OTHER breaking consumer must be listed and handled. The LD matrix `.npz` download task is separate and untouched.

- [ ] **Step 4: Manual end-to-end check (not a committed test)**

Run:
```bash
pixi r python -c "
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED as G
info = G.join_task.meta.harmonization_info
print(info)
assert info is not None and info.build == '19' and info.ref_allele_col == 'NEA' and info.pos_col == 'POS'
print('OK: provenance reaches the finemap generator input')
"
```
Expected: prints the `HarmonizationInfo` and `OK: ...`. This confirms Tasks 6+7 carry provenance harmonize -> join -> `.join_task`. If it fails, a task on the chain rebuilt `FilteredGWASDataMeta` without copying the field; trace `harmonize_task -> join_task`.

- [ ] **Step 5: Commit**

```bash
git add mecfs_bio/asset_generator/ukbb_broad_ld_matrix_generator.py mecfs_bio/asset_generator/fine_mapping_asset_generator.py mecfs_bio/asset_generator/polyfun_explain_fine_mapping_asset_generator.py
git commit -m "feat: tag Broad LD labels and rename them via RenameColsTask (provenance-preserving)"
```

---

## Self-Review

**Spec coverage (Plan 1 scope = spec's D minus the SUSIE runtime check):**
- `HarmonizationInfo` + `FilteredGWASDataMeta` field → Task 1. ✓
- `FASTAMeta`, `HarmonizableReferenceTableMeta` siblings + path branches → Task 2. ✓
- `RenameColsTask` with provenance remap → Task 3. ✓
- Tag FASTA (hg19/hg38) → Task 4; tag 1kg panel → Task 5; tag Broad LD labels → Task 8. ✓
- `GenomeReferenceHarmonizationTask` build-match + stamping → Task 6. ✓
- Spot-2 propagation (`JoinDataFramesTask`) → Task 7. ✓
- LD-label rename via `RenameColsTask` in both generators → Task 8. ✓
- SUSIE build-match assertion → deferred to Plan 2 (documented in the spec's decomposition). ✓ (out of scope)

**Placeholder scan:** no "TBD"/"handle edge cases"; each code step shows real code. Two steps say "read the file first for exact names" (hg38 fasta vars, and confirming `ValidBackend`/`get_extension_and_read_spec_from_format` import location) — verification instructions, not placeholders.

**Type consistency:** `HarmonizationInfo(build, ref_allele_col, pos_col)` identical in Tasks 1/2/3/5/6/8. `FASTAMeta.build` and `HarmonizableReferenceTableMeta.harmonization_info` consistent across Tasks 2/4/5/6/8. `RenameColsTask.create(source_task, asset_id, renames, out_format, backend)` consistent in Tasks 3/8. Ref-allele column is `"NEA"` (gwas / renamed-LD), `"REF"` (panel), `"allele1"` (raw LD) — intentional per the verified orientation.

**Trivial-test scan (per the testing philosophy):** only Task 2 (path resolution — real logic, silent-orphan failure mode) and Task 3 (`RenameColsTask` remap + collision — the novel computation) carry tests. No asset-generator/wiring tests, no field-you-just-set tests, no thin-wrapper tests.

## Risks / notes for the executor
- Task 6's `invoke green` may surface existing harmonization unit tests that construct the task with an untagged fasta/panel; update those to tagged assets or `FASTAMeta`/`HarmonizableReferenceTableMeta` values — the new assertion is meant to reject untagged inputs.
- Metadata-only changes do not alter asset content, and Task 2 proves the tagged assets' store paths are unchanged, so existing built assets stay findable — no rebuilds.
