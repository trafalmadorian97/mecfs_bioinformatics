# Vectorized Cell-Type-Specific S-LDSC Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A source-agnostic, vectorized cell-type-specific S-LDSC Task that reproduces gwaslab's `cell_type_specific` exactly (default) or runs a faster shared-weights approximation, plus cached "prepared" reference assets and generators wiring it to the six standard LDSC-SEG sets.

**Architecture:** Adapter Tasks convert LDSC-format reference directories into a format-neutral parquet contract (prepared annotation set; prepared baseline + weights), built once and reused by every trait. A per-trait Task aligns sumstats to the baseline SNP order once and calls a pure-numpy kernel that runs all cell types' weighted regressions + 200-block jackknife as batched GEMMs. Generators reuse the existing downstream wiring (multiple testing, markdown, labels, plots) with `_vec`-suffixed asset ids.

**Tech Stack:** Python 3.13, numpy, polars 2.0, scipy.stats, attrs, gwaslab (oracle only), pytest. All commands via `pixi r`.

**Spec:** `experiments/claude/design_specs/2026-10-16-vectorized-sldsc-design.md` (read it first; §3 is the exact computation).

## Global Constraints

- Run everything through pixi: `pixi r pytest ...`, `pixi r python ...`, `pixi r invoke green`.
- After each task: `pixi r invoke green 2>&1 | tee /tmp/green_vec_sldsc.log; echo EXIT=$?` and read the log (a bare exit code with testmon "no tests ran" proves nothing).
- Docstrings: no backticks around inline code, no RST.
- Imports at the top of files (never inside functions/tests).
- No mocking/monkeypatching; inject dependencies.
- Return `@frozen` attrs objects, not bare tuples; numpy-array attrs assert ndim/shape/dtype.kind in `__attrs_post_init__`.
- Use polars (not pandas) for new dataframe code; narwhals → polars via `.to_polars()`.
- Use column-name constants (`GWASLAB_*` from `mecfs_bio/constants/gwaslab_constants.py`, new constants defined below), not repeated literals.
- Use `Path`/`PurePath` objects; str only at text boundaries.
- Do not use pytest `match=` with a duplicated error string (omit match).
- Exact mode must reproduce gwaslab to: |Δcoef|/se and |Δz| < 1e-4 on real data; unit-test oracle rtol 1e-5.
- Missing N on any regression SNP raises (intentional divergence from gwaslab's median fill). N is cast to int64 with truncation toward zero, as gwaslab does.
- New generated asset ids carry a `_vec` suffix; the old `CellAnalysisByLDSCTask` and `SLDSCTaskGenerator` ids/behaviour are unchanged.
- Never `find` over the asset store; locate assets via Task.meta + `simple_meta_to_relative_path` + `default_runner_config.yaml`.
- Never commit logs or PDFs; stage explicit paths.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. Duplicate rsIDs in the sumstats among baseline SNPs: gwaslab's merge silently duplicates rows; the new Task must raise instead (test in Task 4).
2. Degenerate Z rows (SE = 0 giving NaN or ±inf Z): must be dropped (NaN) or removed by the χ² filter (inf), with results identical to running without those rows (test in Task 4).
3. Fewer than 200 regression SNPs: n_blocks must become min(n_snp, 200) exactly as ldsc does (test in Task 3).
4. Models spread over multiple control groups and over several chunks: results must not depend on grouping or chunk_size (test in Task 3).
5. Sumstats without an N column: must fail fast with a clear assertion, not a polars ColumnNotFoundError deep in alignment (test in Task 4).

---

## File Structure

Create:

- `mecfs_bio/constants/sldsc_constants.py` — `SLDSCWeighting` Literal; result-table column names.
- `mecfs_bio/build_system/task/sldsc/__init__.py`
- `mecfs_bio/build_system/task/sldsc/prepared_sldsc_reference.py` — the parquet contract: `LdScoreTable`, `SLDSCAnnotationModel`, `PreparedAnnotationSet`, `PreparedBaseline`, readers/writers.
- `mecfs_bio/build_system/task/sldsc/ld_score_files.py` — reading ldsc chromosome-split `.l2.ldscore.gz` / `.l2.M_5_50` files into an `LdScoreTable`.
- `mecfs_bio/build_system/task/sldsc/prepare_ldcts_annotation_set_task.py` — `PrepareLdctsAnnotationSetTask`.
- `mecfs_bio/build_system/task/sldsc/prepare_sldsc_baseline_task.py` — `PrepareSLDSCBaselineTask`.
- `mecfs_bio/build_system/task/sldsc/batched_sldsc_cts.py` — the kernel.
- `mecfs_bio/build_system/task/sldsc/sldsc_regression_inputs.py` — sumstats loading + alignment to prepared references.
- `mecfs_bio/build_system/task/sldsc/vectorized_cell_analysis_by_sldsc_task.py` — `VectorizedCellAnalysisBySLDSCTask`.
- `mecfs_bio/build_system/task_generator/vectorized_sldsc_task_generator.py` — `VectorizedSLDSCTaskGenerator`, `PreparedAnnotationSetRecord`.
- `mecfs_bio/assets/reference_data/linkage_disequilibrium_score_reference_data/prepared/__init__.py`
- `mecfs_bio/assets/reference_data/linkage_disequilibrium_score_reference_data/prepared/prepared_sldsc_reference_sets.py` — seven prepared-asset Task instances.
- `mecfs_bio/asset_generator/concrete_vectorized_sldsc_generator.py` — `standard_vectorized_sldsc_task_generator`.
- `test_mecfs_bio/unit/build_system/task/sldsc/__init__.py`
- `test_mecfs_bio/unit/build_system/task/sldsc/sldsc_fixture.py` — synthetic LDSC-SEG reference + sumstats writer (test helper, not a test module).
- `test_mecfs_bio/unit/build_system/task/sldsc/test_prepared_sldsc_reference.py`
- `test_mecfs_bio/unit/build_system/task/sldsc/test_prepare_tasks.py`
- `test_mecfs_bio/unit/build_system/task/sldsc/test_batched_sldsc_cts.py`
- `test_mecfs_bio/unit/build_system/task/sldsc/test_vectorized_cell_analysis_by_sldsc_task.py`
- `experiments/claude/vectorized_sldsc/vectorized_sldsc_experiment.py`

Modify:

- `mecfs_bio/build_system/task_generator/sldsc_task_generator.py` — extract downstream wiring into a free function; widen `CellAnalysisTaskGroup.cell_analysis_task` to `Task`.
- `mecfs_bio/asset_generator/concrete_sldsc_generator.py` — hoist the three `CellOrTissueLabelRecord`s to module constants for reuse.

---

### Task 1: Constants and the prepared-reference contract

**Files:**
- Create: `mecfs_bio/constants/sldsc_constants.py`
- Create: `mecfs_bio/build_system/task/sldsc/__init__.py` (empty)
- Create: `mecfs_bio/build_system/task/sldsc/prepared_sldsc_reference.py`
- Create: `test_mecfs_bio/unit/build_system/task/sldsc/__init__.py` (empty)
- Test: `test_mecfs_bio/unit/build_system/task/sldsc/test_prepared_sldsc_reference.py`

**Interfaces:**
- Produces (constants): `SLDSCWeighting = Literal["exact", "shared_controls"]`; `SLDSC_CTS_NAME_COL = "Name"`, `SLDSC_CTS_COEF_COL = "Coefficient"`, `SLDSC_CTS_COEF_SE_COL = "Coefficient_std_error"`, `SLDSC_CTS_P_COL = "Coefficient_P_value"`.
- Produces (contract): column constants `SNP_COL="SNP"`, `CHR_COL="CHR"`, `BP_COL="BP"`, `W_LD_COL="w_ld"`; `LdScoreTable(snps: pl.DataFrame, column_ids: tuple[str, ...], sources: tuple[str, ...], matrix: np.ndarray, m_5_50: np.ndarray)`; `SLDSCAnnotationModel(label: str, target_column: str, control_columns: tuple[str, ...])`; `PreparedAnnotationSet(ld: LdScoreTable, models: tuple[SLDSCAnnotationModel, ...])`; `PreparedBaseline(ld: LdScoreTable, w_ld: pl.DataFrame)`; `write_prepared_annotation_set(directory: Path, annotation_set: PreparedAnnotationSet) -> None`; `read_prepared_annotation_set(directory: Path) -> PreparedAnnotationSet`; `write_prepared_baseline(directory: Path, baseline: PreparedBaseline) -> None`; `read_prepared_baseline(directory: Path) -> PreparedBaseline`.

- [ ] **Step 1: Write the constants module**

```python
"""
Constants for cell-type-specific stratified LD score regression (S-LDSC).
"""

from typing import Literal

# exact: per-cell-type regression weights, reproducing ldsc's cell_type_specific.
# shared_controls: one weight vector per control group (target column excluded from the
# weight formula), enabling a Frisch-Waugh shortcut; faster but approximate.
SLDSCWeighting = Literal["exact", "shared_controls"]

# Result-table columns written by ldsc's cell_type_specific and by the vectorized Task.
SLDSC_CTS_NAME_COL = "Name"
SLDSC_CTS_COEF_COL = "Coefficient"
SLDSC_CTS_COEF_SE_COL = "Coefficient_std_error"
SLDSC_CTS_P_COL = "Coefficient_P_value"
```

- [ ] **Step 2: Write the failing tests**

`test_prepared_sldsc_reference.py`:

```python
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from mecfs_bio.build_system.task.sldsc.prepared_sldsc_reference import (
    BP_COL,
    CHR_COL,
    SNP_COL,
    W_LD_COL,
    LdScoreTable,
    PreparedAnnotationSet,
    PreparedBaseline,
    SLDSCAnnotationModel,
    read_prepared_annotation_set,
    read_prepared_baseline,
    write_prepared_annotation_set,
    write_prepared_baseline,
)


def _table(n_cols: int) -> LdScoreTable:
    snps = pl.DataFrame(
        {SNP_COL: ["rs1", "rs2", "rs3"], CHR_COL: [1, 1, 2], BP_COL: [10, 20, 5]}
    )
    return LdScoreTable(
        snps=snps,
        column_ids=tuple(f"ld_{i}" for i in range(n_cols)),
        sources=tuple(f"src:{i}" for i in range(n_cols)),
        matrix=np.arange(3 * n_cols, dtype=np.float32).reshape(3, n_cols),
        m_5_50=np.arange(1, n_cols + 1, dtype=np.float64),
    )


def test_annotation_set_round_trips(tmp_path: Path) -> None:
    original = PreparedAnnotationSet(
        ld=_table(3),
        models=(
            SLDSCAnnotationModel("a", "ld_0", ("ld_2",)),
            SLDSCAnnotationModel("b", "ld_1", ("ld_2",)),
        ),
    )
    write_prepared_annotation_set(tmp_path, original)
    loaded = read_prepared_annotation_set(tmp_path)
    assert loaded.models == original.models
    assert loaded.ld.column_ids == original.ld.column_ids
    np.testing.assert_array_equal(loaded.ld.matrix, original.ld.matrix)
    np.testing.assert_array_equal(loaded.ld.m_5_50, original.ld.m_5_50)
    assert loaded.ld.snps.equals(original.ld.snps)


def test_baseline_round_trips(tmp_path: Path) -> None:
    w_ld = pl.DataFrame({SNP_COL: ["rs1", "rs2", "rs3"], W_LD_COL: [1.0, 2.0, 3.0]})
    original = PreparedBaseline(ld=_table(2), w_ld=w_ld)
    write_prepared_baseline(tmp_path, original)
    loaded = read_prepared_baseline(tmp_path)
    assert loaded.w_ld.equals(original.w_ld)
    np.testing.assert_array_equal(loaded.ld.matrix, original.ld.matrix)


def test_model_referencing_missing_column_is_rejected() -> None:
    with pytest.raises(AssertionError):
        PreparedAnnotationSet(
            ld=_table(2), models=(SLDSCAnnotationModel("a", "ld_0", ("ld_9",)),)
        )


def test_target_among_its_controls_is_rejected() -> None:
    with pytest.raises(AssertionError):
        SLDSCAnnotationModel("a", "ld_0", ("ld_0",))
```

(The round-trip tests guard the one lossy boundary — parquet list columns and float32 — that every consumer depends on.)

- [ ] **Step 3: Run tests to verify they fail**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/sldsc/test_prepared_sldsc_reference.py -v`
Expected: FAIL with ModuleNotFoundError for prepared_sldsc_reference.

- [ ] **Step 4: Implement the contract**

`prepared_sldsc_reference.py`:

```python
"""
Format-neutral on-disk contract for S-LDSC reference data.

A prepared annotation set is a directory holding:
    ld_scores.parquet: SNP, CHR, BP, then one float32 column per unique LD-score source,
        rows sorted by CHR then BP with unique SNPs.
    m_5_50.parquet: one row per LD column: column id, genome-wide M_5_50, source.
    models.parquet: one row per analysis label: label, target column id, control column ids.

A prepared baseline holds ld_scores.parquet and m_5_50.parquet for the baseline annotations
plus w_ld.parquet (SNP and the regression-weight LD score).

Adapter Tasks (one per source format) write this contract; the vectorized S-LDSC Task reads it
and knows nothing about the source format.
"""

from pathlib import Path

import numpy as np
import polars as pl
from attrs import frozen

SNP_COL = "SNP"
CHR_COL = "CHR"
BP_COL = "BP"
W_LD_COL = "w_ld"

LD_SCORES_FILENAME = "ld_scores.parquet"
M_5_50_FILENAME = "m_5_50.parquet"
MODELS_FILENAME = "models.parquet"
W_LD_FILENAME = "w_ld.parquet"

_COLUMN_ID_COL = "column_id"
_M_5_50_COL = "m_5_50"
_SOURCE_COL = "source"
_LABEL_COL = "label"
_TARGET_COL = "target_column"
_CONTROLS_COL = "control_columns"
_KEY_COLS = (SNP_COL, CHR_COL, BP_COL)


@frozen(slots=True)
class LdScoreTable:
    """
    LD scores for a fixed SNP universe.

    snps: SNP, CHR, BP in row order (sorted by CHR then BP, SNP unique).
    column_ids: internal id of each LD column (matrix column order).
    sources: human-readable provenance of each column (e.g. file prefix and column name).
    matrix: (n_snps, n_columns) float32 LD scores.
    m_5_50: (n_columns,) genome-wide M_5_50 per column.
    """

    snps: pl.DataFrame
    column_ids: tuple[str, ...]
    sources: tuple[str, ...]
    matrix: np.ndarray
    m_5_50: np.ndarray

    def __attrs_post_init__(self) -> None:
        assert self.snps.columns == list(_KEY_COLS), (
            f"snps must have columns {_KEY_COLS}, got {self.snps.columns}"
        )
        assert self.snps[SNP_COL].is_unique().all(), "SNP ids must be unique"
        n_cols = len(self.column_ids)
        assert len(set(self.column_ids)) == n_cols, "column ids must be unique"
        assert len(self.sources) == n_cols, "one source per column"
        assert self.matrix.ndim == 2, f"matrix must be 2-D, got {self.matrix.shape}"
        assert self.matrix.shape == (self.snps.height, n_cols), (
            f"matrix shape {self.matrix.shape} != ({self.snps.height}, {n_cols})"
        )
        assert self.matrix.dtype == np.float32, f"matrix must be float32, got {self.matrix.dtype}"
        assert self.m_5_50.ndim == 1 and self.m_5_50.shape[0] == n_cols, (
            f"m_5_50 must have shape ({n_cols},), got {self.m_5_50.shape}"
        )
        assert self.m_5_50.dtype.kind == "f", f"m_5_50 must be float, got {self.m_5_50.dtype}"

    def column_index(self, column_id: str) -> int:
        return self.column_ids.index(column_id)


@frozen(slots=True)
class SLDSCAnnotationModel:
    """One S-LDSC analysis: the reported target column plus the control columns it is fit with."""

    label: str
    target_column: str
    control_columns: tuple[str, ...]

    def __attrs_post_init__(self) -> None:
        assert self.target_column not in self.control_columns, (
            f"model {self.label}: target {self.target_column} is among its own controls"
        )
        assert len(set(self.control_columns)) == len(self.control_columns), (
            f"model {self.label}: duplicate control columns"
        )


@frozen(slots=True)
class PreparedAnnotationSet:
    ld: LdScoreTable
    models: tuple[SLDSCAnnotationModel, ...]

    def __attrs_post_init__(self) -> None:
        labels = [m.label for m in self.models]
        assert len(set(labels)) == len(labels), "model labels must be unique"
        known = set(self.ld.column_ids)
        for model in self.models:
            missing = {model.target_column, *model.control_columns} - known
            assert not missing, f"model {model.label} references unknown columns {missing}"


@frozen(slots=True)
class PreparedBaseline:
    """Baseline LD scores plus the regression-weight LD score (w_ld: SNP, w_ld; SNP unique)."""

    ld: LdScoreTable
    w_ld: pl.DataFrame

    def __attrs_post_init__(self) -> None:
        assert self.w_ld.columns == [SNP_COL, W_LD_COL], (
            f"w_ld must have columns {[SNP_COL, W_LD_COL]}, got {self.w_ld.columns}"
        )
        assert self.w_ld[SNP_COL].is_unique().all(), "w_ld SNP ids must be unique"


def write_prepared_annotation_set(
    directory: Path, annotation_set: PreparedAnnotationSet
) -> None:
    _write_ld_score_table(directory, annotation_set.ld)
    pl.DataFrame(
        {
            _LABEL_COL: [m.label for m in annotation_set.models],
            _TARGET_COL: [m.target_column for m in annotation_set.models],
            _CONTROLS_COL: [list(m.control_columns) for m in annotation_set.models],
        },
        schema={
            _LABEL_COL: pl.String,
            _TARGET_COL: pl.String,
            _CONTROLS_COL: pl.List(pl.String),
        },
    ).write_parquet(directory / MODELS_FILENAME)


def read_prepared_annotation_set(directory: Path) -> PreparedAnnotationSet:
    models = pl.read_parquet(directory / MODELS_FILENAME)
    return PreparedAnnotationSet(
        ld=_read_ld_score_table(directory),
        models=tuple(
            SLDSCAnnotationModel(
                label=row[_LABEL_COL],
                target_column=row[_TARGET_COL],
                control_columns=tuple(row[_CONTROLS_COL]),
            )
            for row in models.iter_rows(named=True)
        ),
    )


def write_prepared_baseline(directory: Path, baseline: PreparedBaseline) -> None:
    _write_ld_score_table(directory, baseline.ld)
    baseline.w_ld.write_parquet(directory / W_LD_FILENAME)


def read_prepared_baseline(directory: Path) -> PreparedBaseline:
    return PreparedBaseline(
        ld=_read_ld_score_table(directory),
        w_ld=pl.read_parquet(directory / W_LD_FILENAME),
    )


def _write_ld_score_table(directory: Path, table: LdScoreTable) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    values = pl.from_numpy(table.matrix, schema=list(table.column_ids), orient="row")
    pl.concat([table.snps, values], how="horizontal").write_parquet(
        directory / LD_SCORES_FILENAME, compression="zstd"
    )
    pl.DataFrame(
        {
            _COLUMN_ID_COL: list(table.column_ids),
            _M_5_50_COL: table.m_5_50,
            _SOURCE_COL: list(table.sources),
        }
    ).write_parquet(directory / M_5_50_FILENAME)


def _read_ld_score_table(directory: Path) -> LdScoreTable:
    m_frame = pl.read_parquet(directory / M_5_50_FILENAME)
    column_ids = tuple(m_frame[_COLUMN_ID_COL].to_list())
    frame = pl.read_parquet(directory / LD_SCORES_FILENAME)
    return LdScoreTable(
        snps=frame.select(_KEY_COLS),
        column_ids=column_ids,
        sources=tuple(m_frame[_SOURCE_COL].to_list()),
        matrix=frame.select(list(column_ids)).to_numpy().astype(np.float32, copy=False),
        m_5_50=m_frame[_M_5_50_COL].to_numpy().astype(np.float64),
    )
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/sldsc/test_prepared_sldsc_reference.py -v`
Expected: 4 passed.

- [ ] **Step 6: Green and commit**

```bash
pixi r invoke green 2>&1 | tee /tmp/green_vec_sldsc.log; echo EXIT=$?
git add mecfs_bio/constants/sldsc_constants.py mecfs_bio/build_system/task/sldsc/__init__.py mecfs_bio/build_system/task/sldsc/prepared_sldsc_reference.py test_mecfs_bio/unit/build_system/task/sldsc/__init__.py test_mecfs_bio/unit/build_system/task/sldsc/test_prepared_sldsc_reference.py
git commit -m "Add prepared S-LDSC reference contract and constants

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: LD-score file reader, synthetic fixture, and the two prepare Tasks

**Files:**
- Create: `mecfs_bio/build_system/task/sldsc/ld_score_files.py`
- Create: `mecfs_bio/build_system/task/sldsc/prepare_ldcts_annotation_set_task.py`
- Create: `mecfs_bio/build_system/task/sldsc/prepare_sldsc_baseline_task.py`
- Create: `test_mecfs_bio/unit/build_system/task/sldsc/sldsc_fixture.py`
- Test: `test_mecfs_bio/unit/build_system/task/sldsc/test_prepare_tasks.py`

**Interfaces:**
- Consumes: Task 1 contract.
- Produces: `read_ld_score_sources(root: Path, prefixes: Sequence[str], *, drop_zero_variance: bool) -> LdScoreTable`; `PrepareLdctsAnnotationSetTask.create(asset_id: str, source_directory_task: Task, ldcts_filename: str) -> PrepareLdctsAnnotationSetTask`; `PrepareSLDSCBaselineTask.create(asset_id: str, baseline_directory_task: Task, baseline_pattern: str, weights_directory_task: Task, weights_pattern: str) -> PrepareSLDSCBaselineTask`. Both return a `DirectoryAsset` holding the Task 1 contract.
- Produces (test helper): `write_sldsc_fixture(root: Path, *, seed: int = 0) -> SLDSCFixture` with fields `baseline_dir: Path`, `weights_dir: Path`, `cts_dir: Path`, `sumstats_frame: pd.DataFrame`; constants `FIXTURE_BASELINE_PATTERN = "baseline/baseline.@"`, `FIXTURE_WEIGHTS_PATTERN = "weights/weights.@"`, `FIXTURE_LDCTS_FILENAME = "Set.ldcts"`, `FIXTURE_LABELS = ("Cell_A", "Cell_B", "Cell_C")`; `sumstats_from_frame(frame: pd.DataFrame) -> gl.Sumstats`.

Reproduction notes (from gwaslab/extension/ldsc/ldsc_parse.py):
- `sub_chr(s, chrom)`: if "@" not in s, append "@"; then replace "@" with the chromosome number. Present chromosomes are those 1..22 whose file exists.
- `ldscore()`: concatenate chromosomes, `sort_values(["CHR", "BP"])`, drop CHR/BP, `drop_duplicates(subset="SNP")` (keep first). Use polars `sort(..., maintain_order=True)` + `unique(subset=SNP, keep="first", maintain_order=True)`.
- `ldscore_fromlist`: every file must have the identical SNP column after that step.
- `M()`: first line of `<prefix><chr>.l2.M_5_50`, whitespace-split floats, summed over present chromosomes.
- `_check_variance` drops baseline columns with zero (pandas, ddof=1) variance over the full reference, and their M. Only the baseline, never cts columns.
- Files are tab-separated with header `CHR SNP BP <cols...>`; parse value columns as Float64.

- [ ] **Step 1: Write the fixture helper**

`sldsc_fixture.py`:

```python
"""
Synthetic LDSC-SEG-format reference data and sumstats for S-LDSC unit tests.

Three chromosomes, 1,000 reference SNPs each. The baseline has a total-LD column, two
annotations, and an all-zero annotation (exercises ldsc's zero-variance drop). One shared
control plus three cell-type targets; Cell_A carries a planted enrichment. Values are printed
with 3 decimals, like the real files. The sumstats omit some reference SNPs, add SNPs absent
from the reference (so gwaslab's smart_merge takes its pd.merge path, as in production), and
include a few chi-square outliers removed by the max(0.001 N, 80) filter.
"""

import gzip
from pathlib import Path

import gwaslab as gl
import numpy as np
import pandas as pd
from attrs import frozen

FIXTURE_BASELINE_PATTERN = "baseline/baseline.@"
FIXTURE_WEIGHTS_PATTERN = "weights/weights.@"
FIXTURE_LDCTS_FILENAME = "Set.ldcts"
FIXTURE_LABELS = ("Cell_A", "Cell_B", "Cell_C")
FIXTURE_N = 50_000

_CHROMS = (1, 2, 3)
_SNPS_PER_CHROM = 1_000
_BASELINE_COLS = ("baseL2", "AnnotAL2", "AnnotBL2", "ZeroL2")
_BASELINE_M = (5000.0, 1000.0, 500.0, 0.0)
_CONTROL_M = 2000.0
_TARGET_M = (300.0, 400.0, 500.0)


@frozen(slots=True)
class SLDSCFixture:
    baseline_dir: Path
    weights_dir: Path
    cts_dir: Path
    sumstats_frame: pd.DataFrame


def _write_ld_file(path: Path, chrom: int, rsids: list[str], bps: np.ndarray, cols: dict[str, np.ndarray]) -> None:
    with gzip.open(path, "wt") as f:
        f.write("\t".join(["CHR", "SNP", "BP", *cols]) + "\n")
        for i, rsid in enumerate(rsids):
            values = "\t".join(f"{cols[c][i]:.3f}" for c in cols)
            f.write(f"{chrom}\t{rsid}\t{bps[i]}\t{values}\n")


def write_sldsc_fixture(root: Path, *, seed: int = 0) -> SLDSCFixture:
    rng = np.random.default_rng(seed)
    baseline_dir = root / "baseline_root"
    weights_dir = root / "weights_root"
    cts_dir = root / "cts_root"
    for d in (baseline_dir / "baseline", weights_dir / "weights", cts_dir / "cts"):
        d.mkdir(parents=True)
    rows = []
    for chrom in _CHROMS:
        rsids = [f"rs{chrom}{i:05d}" for i in range(_SNPS_PER_CHROM)]
        bps = 10_000 + 100 * np.arange(_SNPS_PER_CHROM)
        base = {
            "baseL2": rng.uniform(1, 50, _SNPS_PER_CHROM),
            "AnnotAL2": rng.uniform(0, 10, _SNPS_PER_CHROM),
            "AnnotBL2": rng.uniform(0, 5, _SNPS_PER_CHROM),
            "ZeroL2": np.zeros(_SNPS_PER_CHROM),
        }
        _write_ld_file(baseline_dir / f"baseline/baseline.{chrom}.l2.ldscore.gz", chrom, rsids, bps, base)
        (baseline_dir / f"baseline/baseline.{chrom}.l2.M_5_50").write_text(
            "\t".join(str(m / len(_CHROMS)) for m in _BASELINE_M) + "\n"
        )
        _write_ld_file(
            weights_dir / f"weights/weights.{chrom}.l2.ldscore.gz", chrom, rsids, bps,
            {"L2": rng.uniform(1, 40, _SNPS_PER_CHROM)},
        )
        control = rng.uniform(0, 20, _SNPS_PER_CHROM)
        _write_ld_file(cts_dir / f"cts/Set.control.{chrom}.l2.ldscore.gz", chrom, rsids, bps, {"L2": control})
        (cts_dir / f"cts/Set.control.{chrom}.l2.M_5_50").write_text(f"{_CONTROL_M / len(_CHROMS)}\n")
        targets = [rng.uniform(0, 10, _SNPS_PER_CHROM) for _ in FIXTURE_LABELS]
        for j, target in enumerate(targets, start=1):
            _write_ld_file(cts_dir / f"cts/Set.{j}.{chrom}.l2.ldscore.gz", chrom, rsids, bps, {"L2": target})
            (cts_dir / f"cts/Set.{j}.{chrom}.l2.M_5_50").write_text(f"{_TARGET_M[j - 1] / len(_CHROMS)}\n")
        expected_chi2 = 1.0 + 0.02 * base["baseL2"] + 0.08 * targets[0]
        z = rng.normal(0.0, 1.0, _SNPS_PER_CHROM) * np.sqrt(expected_chi2)
        rows.extend(
            {"rsID": r, "CHR": chrom, "POS": int(bp), "Z": float(zz)}
            for r, bp, zz in zip(rsids, bps, z)
        )
    (cts_dir / FIXTURE_LDCTS_FILENAME).write_text(
        "\n".join(f"{label}\tcts/Set.{j}.,cts/Set.control." for j, label in enumerate(FIXTURE_LABELS, start=1))
        + "\n"
    )
    frame = pd.DataFrame(rows)
    frame = frame.iloc[50:].reset_index(drop=True)  # some reference SNPs absent from the sumstats
    extra = pd.DataFrame(
        {"rsID": [f"rs_extra{i}" for i in range(20)], "CHR": 1, "POS": 9_000_000 + np.arange(20), "Z": rng.normal(0, 1, 20)}
    )
    frame = pd.concat([frame, extra], ignore_index=True)
    frame.loc[[100, 200, 300], "Z"] = 12.0  # chi2 = 144 > max(0.001 * N, 80)
    se = np.full(len(frame), 0.01)
    frame = frame.assign(
        BETA=frame["Z"] * se, SE=se, EA="A", NEA="G", P=0.5, N=FIXTURE_N
    ).drop(columns="Z")
    return SLDSCFixture(baseline_dir=baseline_dir, weights_dir=weights_dir, cts_dir=cts_dir, sumstats_frame=frame)


def sumstats_from_frame(frame: pd.DataFrame) -> gl.Sumstats:
    return gl.Sumstats(
        frame, rsid="rsID", chrom="CHR", pos="POS", ea="EA", nea="NEA",
        beta="BETA", se="SE", p="P", n="N", build="19", verbose=False,
    )
```

(Format with `pixi r invoke green` later; the long lines above will be reformatted.)

- [ ] **Step 2: Write the failing tests**

`test_prepare_tasks.py`:

```python
from pathlib import Path

import numpy as np

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.sldsc.prepare_ldcts_annotation_set_task import (
    PrepareLdctsAnnotationSetTask,
)
from mecfs_bio.build_system.task.sldsc.prepare_sldsc_baseline_task import (
    PrepareSLDSCBaselineTask,
)
from mecfs_bio.build_system.task.sldsc.prepared_sldsc_reference import (
    read_prepared_annotation_set,
    read_prepared_baseline,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from test_mecfs_bio.unit.build_system.task.sldsc.sldsc_fixture import (
    FIXTURE_BASELINE_PATTERN,
    FIXTURE_LABELS,
    FIXTURE_LDCTS_FILENAME,
    FIXTURE_WEIGHTS_PATTERN,
    write_sldsc_fixture,
)


def _dir_task(asset_id: str) -> FakeTask:
    return FakeTask(
        ReferenceDataDirectoryMeta(
            group="g", sub_group="s", sub_folder=Path("extracted"), id=AssetId(asset_id)
        )
    )


def test_prepare_annotation_set_stores_shared_control_once(tmp_path: Path) -> None:
    fixture = write_sldsc_fixture(tmp_path)
    task = PrepareLdctsAnnotationSetTask.create(
        asset_id="prepared_set",
        source_directory_task=_dir_task("cts"),
        ldcts_filename=FIXTURE_LDCTS_FILENAME,
    )

    def fetch(asset_id: AssetId) -> Asset:
        assert asset_id == "cts"
        return DirectoryAsset(fixture.cts_dir)

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)
    prepared = read_prepared_annotation_set(result.path)
    assert [m.label for m in prepared.models] == list(FIXTURE_LABELS)
    controls = {m.control_columns for m in prepared.models}
    assert len(controls) == 1  # one shared control column, referenced by every model
    assert prepared.ld.matrix.shape == (3_000, len(FIXTURE_LABELS) + 1)
    np.testing.assert_allclose(prepared.ld.m_5_50.sum(), 2000.0 + 300.0 + 400.0 + 500.0)


def test_prepare_baseline_drops_zero_variance_column(tmp_path: Path) -> None:
    fixture = write_sldsc_fixture(tmp_path)
    task = PrepareSLDSCBaselineTask.create(
        asset_id="prepared_baseline",
        baseline_directory_task=_dir_task("baseline"),
        baseline_pattern=FIXTURE_BASELINE_PATTERN,
        weights_directory_task=_dir_task("weights"),
        weights_pattern=FIXTURE_WEIGHTS_PATTERN,
    )

    def fetch(asset_id: AssetId) -> Asset:
        return DirectoryAsset({"baseline": fixture.baseline_dir, "weights": fixture.weights_dir}[asset_id])

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)
    prepared = read_prepared_baseline(result.path)
    assert not any("ZeroL2" in s for s in prepared.ld.sources)
    assert prepared.ld.matrix.shape == (3_000, 3)
    assert prepared.w_ld.height == 3_000
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/sldsc/test_prepare_tasks.py -v`
Expected: FAIL with ModuleNotFoundError for prepare_ldcts_annotation_set_task.

- [ ] **Step 4: Implement the file reader**

`ld_score_files.py`:

```python
"""
Readers for ldsc's chromosome-split LD-score files (<prefix><chr>.l2.ldscore.gz and
<prefix><chr>.l2.M_5_50), reproducing the parsing in ldsc's ldsc_parse module: chromosome
substitution, CHR/BP sort, first-occurrence SNP de-duplication, identical SNP columns across
files, and M_5_50 summed over the chromosomes present.
"""

from pathlib import Path
from typing import Sequence

import numpy as np
import polars as pl

from mecfs_bio.build_system.task.sldsc.prepared_sldsc_reference import (
    BP_COL,
    CHR_COL,
    SNP_COL,
    LdScoreTable,
)

_MAX_CHROM = 22
_LD_SUFFIX = ".l2.ldscore.gz"
_M_SUFFIX = ".l2.M_5_50"
_KEY_COLS = (SNP_COL, CHR_COL, BP_COL)


def chromosome_file_prefix(prefix: str, chrom: int) -> str:
    """ldsc's sub_chr: append @ if absent, then substitute the chromosome number."""
    with_marker = prefix if "@" in prefix else prefix + "@"
    return with_marker.replace("@", str(chrom))


def present_chromosomes(root: Path, prefix: str) -> list[int]:
    return [
        c
        for c in range(1, _MAX_CHROM + 1)
        if (root / (chromosome_file_prefix(prefix, c) + _LD_SUFFIX)).is_file()
    ]


def read_ld_score_prefix(root: Path, prefix: str) -> pl.DataFrame:
    chroms = present_chromosomes(root, prefix)
    assert chroms, f"no LD-score files for prefix {prefix} under {root}"
    frames = [
        pl.read_csv(root / (chromosome_file_prefix(prefix, c) + _LD_SUFFIX), separator="\t")
        for c in chroms
    ]
    frame = pl.concat(frames)
    value_cols = [c for c in frame.columns if c not in _KEY_COLS]
    return (
        frame.select(
            pl.col(SNP_COL).cast(pl.String),
            pl.col(CHR_COL).cast(pl.Int64),
            pl.col(BP_COL).cast(pl.Int64),
            *[pl.col(c).cast(pl.Float64) for c in value_cols],
        )
        .sort([CHR_COL, BP_COL], maintain_order=True)
        .unique(subset=SNP_COL, keep="first", maintain_order=True)
    )


def read_m_5_50_prefix(root: Path, prefix: str) -> np.ndarray:
    per_chrom = [
        [
            float(v)
            for v in (root / (chromosome_file_prefix(prefix, c) + _M_SUFFIX))
            .read_text()
            .splitlines()[0]
            .split()
        ]
        for c in present_chromosomes(root, prefix)
    ]
    return np.sum(np.array(per_chrom, dtype=np.float64), axis=0)


def read_ld_score_sources(
    root: Path, prefixes: Sequence[str], *, drop_zero_variance: bool
) -> LdScoreTable:
    """
    Read and side-by-side concatenate the LD scores of several file prefixes (ldsc's
    ldscore_fromlist). Column ids are ld_0, ld_1, ... in prefix-then-column order; sources
    record prefix:column. Streams one prefix at a time into a preallocated float32 matrix so
    that hundreds of prefixes do not all sit in memory as dataframes.

    drop_zero_variance: drop columns whose values have zero variance (ldsc's _check_variance,
    applied to the baseline only).
    """
    assert prefixes, "need at least one prefix"
    snps: pl.DataFrame | None = None
    blocks: list[np.ndarray] = []
    sources: list[str] = []
    m_values: list[np.ndarray] = []
    for prefix in prefixes:
        frame = read_ld_score_prefix(root, prefix)
        if snps is None:
            snps = frame.select(_KEY_COLS)
        assert frame[SNP_COL].equals(snps[SNP_COL]), (
            f"LD scores for {prefix} do not have the same SNPs as {prefixes[0]}"
        )
        value_cols = [c for c in frame.columns if c not in _KEY_COLS]
        m = read_m_5_50_prefix(root, prefix)
        assert m.shape[0] == len(value_cols), (
            f"{prefix}: {m.shape[0]} M_5_50 values for {len(value_cols)} LD columns"
        )
        values = frame.select(value_cols).to_numpy()
        keep = np.var(values, axis=0, ddof=1) != 0 if drop_zero_variance else np.ones(len(value_cols), dtype=bool)
        blocks.append(values[:, keep].astype(np.float32))
        sources.extend(f"{prefix}:{c}" for c, k in zip(value_cols, keep) if k)
        m_values.append(m[keep])
    assert snps is not None
    matrix = np.concatenate(blocks, axis=1) if len(blocks) > 1 else blocks[0]
    return LdScoreTable(
        snps=snps,
        column_ids=tuple(f"ld_{i}" for i in range(matrix.shape[1])),
        sources=tuple(sources),
        matrix=matrix,
        m_5_50=np.concatenate(m_values),
    )
```

Memory note: for the 490-prefix chromosome set each block is a (1.2M, 1) float32 column (~5 MB), so the block list stays at ~2.4 GB and the concatenate peaks at ~4.8 GB — within the 16 GB budget.

- [ ] **Step 5: Implement the two prepare Tasks**

`prepare_ldcts_annotation_set_task.py`:

```python
"""
Adapter Task: convert an LDSC-SEG style directory (a .ldcts index plus chromosome-split LD-score
files) into a prepared S-LDSC annotation set (see prepared_sldsc_reference).

Each .ldcts row is "label<TAB>target_prefix,control_prefix". Every unique prefix is read once,
so a control shared by every row is stored once and referenced by every model.
"""

from pathlib import Path, PurePath

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_cell_analysis_by_sldsc import (
    LDCTSFileEntry,
    ldcts_raw_to_processed,
)
from mecfs_bio.build_system.task.sldsc.ld_score_files import read_ld_score_sources
from mecfs_bio.build_system.task.sldsc.prepared_sldsc_reference import (
    PreparedAnnotationSet,
    SLDSCAnnotationModel,
    write_prepared_annotation_set,
)
from mecfs_bio.build_system.wf.base_wf import WF

PREPARED_SLDSC_SUB_FOLDER = PurePath("prepared_sldsc")


@frozen(slots=True)
class PrepareLdctsAnnotationSetTask(Task):
    """
    source_directory_task: directory containing ldcts_filename and the LD-score files its
        rows reference (paths in the index are relative to this directory).
    ldcts_filename: name of the .ldcts index inside that directory.
    """

    meta: Meta
    source_directory_task: Task
    ldcts_filename: str

    @property
    def deps(self) -> list[Task]:
        return [self.source_directory_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source = fetch(self.source_directory_task.asset_id)
        assert isinstance(source, DirectoryAsset)
        entries = _read_ldcts_entries(source.path / self.ldcts_filename)
        prefixes = _unique_prefixes(entries)
        ld = read_ld_score_sources(source.path, prefixes, drop_zero_variance=False)
        assert len(ld.column_ids) == len(prefixes), "each cts prefix must hold exactly one LD column"
        column_of = dict(zip(prefixes, ld.column_ids))
        annotation_set = PreparedAnnotationSet(
            ld=ld, models=tuple(_model_from_entry(e, column_of) for e in entries)
        )
        out_dir = scratch_dir / "prepared_annotation_set"
        write_prepared_annotation_set(out_dir, annotation_set)
        return DirectoryAsset(out_dir)

    @classmethod
    def create(
        cls, asset_id: str, source_directory_task: Task, ldcts_filename: str
    ) -> "PrepareLdctsAnnotationSetTask":
        source_meta = source_directory_task.meta
        assert isinstance(source_meta, ReferenceDataDirectoryMeta)
        return cls(
            meta=ReferenceDataDirectoryMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=PREPARED_SLDSC_SUB_FOLDER,
                id=AssetId(asset_id),
            ),
            source_directory_task=source_directory_task,
            ldcts_filename=ldcts_filename,
        )


def _read_ldcts_entries(path: Path) -> list[LDCTSFileEntry]:
    assert path.is_file(), f"missing .ldcts index {path}"
    entries = ldcts_raw_to_processed(path.read_text().strip())
    assert entries, f"empty .ldcts index {path}"
    return entries


def _unique_prefixes(entries: list[LDCTSFileEntry]) -> list[str]:
    seen: dict[str, None] = {}
    for entry in entries:
        seen.setdefault(str(entry.target_gene_data_path), None)
        seen.setdefault(str(entry.control_gene_data_path), None)
    return list(seen)


def _model_from_entry(entry: LDCTSFileEntry, column_of: dict[str, str]) -> SLDSCAnnotationModel:
    return SLDSCAnnotationModel(
        label=entry.label,
        target_column=column_of[str(entry.target_gene_data_path)],
        control_columns=(column_of[str(entry.control_gene_data_path)],),
    )
```

Check while implementing: `str(Path("cts/Set.1."))` must keep the trailing dot (it does on POSIX); if not, keep the raw string from the index instead.

`prepare_sldsc_baseline_task.py`:

```python
"""
Adapter Task: convert ldsc baseline LD scores (chromosome-split, @ pattern) and regression-weight
LD scores into a prepared S-LDSC baseline (see prepared_sldsc_reference). Zero-variance baseline
columns are dropped here, as ldsc's _check_variance does on the full reference.
"""

from pathlib import Path

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.sldsc.ld_score_files import (
    read_ld_score_prefix,
    read_ld_score_sources,
)
from mecfs_bio.build_system.task.sldsc.prepare_ldcts_annotation_set_task import (
    PREPARED_SLDSC_SUB_FOLDER,
)
from mecfs_bio.build_system.task.sldsc.prepared_sldsc_reference import (
    BP_COL,
    CHR_COL,
    SNP_COL,
    W_LD_COL,
    PreparedBaseline,
    write_prepared_baseline,
)
from mecfs_bio.build_system.wf.base_wf import WF


@frozen(slots=True)
class PrepareSLDSCBaselineTask(Task):
    """
    baseline_pattern: path inside the baseline directory with @ standing for the chromosome,
        e.g. baseline_v1.2/baseline.@
    weights_pattern: same, for the regression-weight LD scores (exactly one LD column).
    """

    meta: Meta
    baseline_directory_task: Task
    baseline_pattern: str
    weights_directory_task: Task
    weights_pattern: str

    @property
    def deps(self) -> list[Task]:
        return [self.baseline_directory_task, self.weights_directory_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        baseline_dir = _fetch_directory(fetch, self.baseline_directory_task.asset_id)
        weights_dir = _fetch_directory(fetch, self.weights_directory_task.asset_id)
        baseline = PreparedBaseline(
            ld=read_ld_score_sources(baseline_dir, [self.baseline_pattern], drop_zero_variance=True),
            w_ld=_read_weight_ld(weights_dir, self.weights_pattern),
        )
        out_dir = scratch_dir / "prepared_baseline"
        write_prepared_baseline(out_dir, baseline)
        return DirectoryAsset(out_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        baseline_directory_task: Task,
        baseline_pattern: str,
        weights_directory_task: Task,
        weights_pattern: str,
    ) -> "PrepareSLDSCBaselineTask":
        source_meta = baseline_directory_task.meta
        assert isinstance(source_meta, ReferenceDataDirectoryMeta)
        return cls(
            meta=ReferenceDataDirectoryMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=PREPARED_SLDSC_SUB_FOLDER,
                id=AssetId(asset_id),
            ),
            baseline_directory_task=baseline_directory_task,
            baseline_pattern=baseline_pattern,
            weights_directory_task=weights_directory_task,
            weights_pattern=weights_pattern,
        )


def _fetch_directory(fetch: Fetch, asset_id: AssetId) -> Path:
    asset = fetch(asset_id)
    assert isinstance(asset, DirectoryAsset)
    return asset.path


def _read_weight_ld(root: Path, pattern: str) -> pl.DataFrame:
    frame = read_ld_score_prefix(root, pattern)
    value_cols = [c for c in frame.columns if c not in (SNP_COL, CHR_COL, BP_COL)]
    assert len(value_cols) == 1, f"weight LD scores must have one column, got {value_cols}"
    return frame.select(pl.col(SNP_COL), pl.col(value_cols[0]).alias(W_LD_COL))
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/sldsc/test_prepare_tasks.py -v`
Expected: 2 passed.

- [ ] **Step 7: Green and commit**

```bash
pixi r invoke green 2>&1 | tee /tmp/green_vec_sldsc.log; echo EXIT=$?
git add mecfs_bio/build_system/task/sldsc/ld_score_files.py mecfs_bio/build_system/task/sldsc/prepare_ldcts_annotation_set_task.py mecfs_bio/build_system/task/sldsc/prepare_sldsc_baseline_task.py test_mecfs_bio/unit/build_system/task/sldsc/sldsc_fixture.py test_mecfs_bio/unit/build_system/task/sldsc/test_prepare_tasks.py
git commit -m "Add LDSC-SEG and baseline adapter Tasks for prepared S-LDSC references

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Kernel, exact mode (oracle: gwaslab's Hsq)

**Files:**
- Create: `mecfs_bio/build_system/task/sldsc/batched_sldsc_cts.py`
- Test: `test_mecfs_bio/unit/build_system/task/sldsc/test_batched_sldsc_cts.py`

**Interfaces:**
- Consumes: `SLDSCWeighting` (Task 1).
- Produces: `SLDSCRegressionInputs(chi2, n, w_ld, baseline_ld, baseline_m, annotation_ld, annotation_m)`; `SLDSCModel(label: str, target_index: int, control_indices: tuple[int, ...])`; `BatchedSLDSCResult(labels: tuple[str, ...], coef, coef_se, p)`; `batched_sldsc_cts(inputs: SLDSCRegressionInputs, models: Sequence[SLDSCModel], *, weighting: SLDSCWeighting, n_blocks: int = LDSC_DEFAULT_N_BLOCKS, chunk_size: int = DEFAULT_CHUNK_SIZE) -> BatchedSLDSCResult` (results in input model order, unsorted). In this task `weighting="shared_controls"` raises `NotImplementedError`; Task 5 fills it in.

Math (spec §3): per model c, unscaled LD columns L_c = [t_c, controls, baseline]; x_tot,c = row-sum of L_c; M_tot,c = sum of their M_5_50; h²_agg,c = clip(M_tot,c (mean χ² − 1) / mean(x_tot,c N), 0, 1); w_c = 1 / (2 (1 + h²_agg,c N max(x_tot,c, 1) / M_tot,c)²) / max(w_ld, 1). Design [t_c, controls, baseline] · N / N̄ plus an intercept column; WLS with weight w_c; ldsc separators floor(linspace(0, S, n_blocks + 1)) with n_blocks = min(S, 200); delete values from (totals − block); pseudo = nb·est − (nb−1)·delete; se = sqrt(var(pseudo, ddof=1)/nb)/N̄; coef = est/N̄; p = norm.sf(coef/se).

Batched exact Gram per block (A = shared scaled design [controls, baseline, 1], q columns; T = scaled targets for a chunk of k models; W = (s_b, k) weights): the shared q×q part for all k models is one GEMM P_bᵀ W_b where P_b holds the upper-triangle products A_i ⊗ A_i (s_b × q(q+1)/2); cross terms A_bᵀ(W_b ∘ T_b); diagonal Σ W T²; right-hand sides A_bᵀ(W_b ∘ y) and Σ W T y.

- [ ] **Step 1: Write the failing tests**

`test_batched_sldsc_cts.py`:

```python
import numpy as np
import pytest
from gwaslab.extension.ldsc.ldsc_regressions import Hsq

from mecfs_bio.build_system.task.sldsc.batched_sldsc_cts import (
    SLDSCModel,
    SLDSCRegressionInputs,
    batched_sldsc_cts,
)


def _random_inputs(n_snps: int, seed: int = 0) -> SLDSCRegressionInputs:
    rng = np.random.default_rng(seed)
    baseline_ld = np.column_stack(
        [rng.uniform(1, 50, n_snps), rng.uniform(0, 10, n_snps), rng.uniform(0, 5, n_snps)]
    )
    annotation_ld = rng.uniform(0, 10, (n_snps, 5)).astype(np.float32)
    n = np.full(n_snps, 40_000.0) + rng.integers(0, 2_000, n_snps)
    chi2 = rng.chisquare(1, n_snps) * (1.0 + 0.02 * baseline_ld[:, 0] + 0.1 * annotation_ld[:, 0])
    return SLDSCRegressionInputs(
        chi2=chi2,
        n=n,
        w_ld=rng.uniform(1, 40, n_snps),
        baseline_ld=baseline_ld,
        baseline_m=np.array([5000.0, 1000.0, 500.0]),
        annotation_ld=annotation_ld,
        annotation_m=np.array([300.0, 400.0, 500.0, 2000.0, 1500.0]),
    )


# Two control groups: columns 3 and 4 are controls; 0, 1, 2 are targets.
_MODELS = (
    SLDSCModel("a", 0, (3,)),
    SLDSCModel("b", 1, (3,)),
    SLDSCModel("c", 2, (4,)),
    SLDSCModel("d", 0, (3, 4)),
)


def _gwaslab_oracle(inputs: SLDSCRegressionInputs, model: SLDSCModel, n_blocks: int) -> tuple[float, float]:
    cols = [model.target_index, *model.control_indices]
    ref_ld = np.column_stack([inputs.annotation_ld[:, cols].astype(np.float64), inputs.baseline_ld])
    m = np.concatenate([inputs.annotation_m[cols], inputs.baseline_m]).reshape(1, -1)
    s = inputs.chi2.shape[0]
    hsq = Hsq(
        inputs.chi2.reshape(s, 1), ref_ld, inputs.w_ld.reshape(s, 1), inputs.n.reshape(s, 1), m,
        n_blocks=min(s, n_blocks), intercept=None, twostep=None, old_weights=True,
    )
    return float(hsq.coef[0]), float(hsq.coef_se[0])


@pytest.mark.parametrize("n_snps", [3_000, 150])  # 150 < 200 exercises n_blocks = min(S, 200)
def test_exact_matches_gwaslab_hsq(n_snps: int) -> None:
    inputs = _random_inputs(n_snps)
    result = batched_sldsc_cts(inputs, _MODELS, weighting="exact", chunk_size=2)
    for i, model in enumerate(_MODELS):
        coef, se = _gwaslab_oracle(inputs, model, n_blocks=200)
        np.testing.assert_allclose(result.coef[i], coef, rtol=1e-8)
        np.testing.assert_allclose(result.coef_se[i], se, rtol=1e-8)
    assert result.labels == tuple(m.label for m in _MODELS)


def test_exact_is_independent_of_chunk_size() -> None:
    inputs = _random_inputs(3_000)
    small = batched_sldsc_cts(inputs, _MODELS, weighting="exact", chunk_size=1)
    large = batched_sldsc_cts(inputs, _MODELS, weighting="exact", chunk_size=64)
    np.testing.assert_allclose(small.coef, large.coef, rtol=1e-12)
    np.testing.assert_allclose(small.coef_se, large.coef_se, rtol=1e-12)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/sldsc/test_batched_sldsc_cts.py -v`
Expected: FAIL with ModuleNotFoundError for batched_sldsc_cts.

- [ ] **Step 3: Implement the kernel (exact mode)**

`batched_sldsc_cts.py`:

```python
"""
Batched cell-type-specific stratified LD score regression.

Reproduces ldsc's cell_type_specific (Hsq with old_weights=True: a single weighted least
squares, no IRWLS, no two-step) for many models at once. Each model regresses chi-square on
[target, controls, baseline] LD scores scaled by N / mean(N), plus an intercept, and reports the
target coefficient with a 200-block jackknife standard error.

Models sharing a control set share the design apart from the target column, so they are
processed together. In exact mode each model keeps its own regression weights (the weights
depend on the model's total LD and total M, which include the target); the per-block normal
equations of a chunk of models are formed with a few large matrix products rather than one
regression at a time. ldsc normalizes the square-root weights to sum to one before weighting
the design and response; that scalar cancels in every solve, so the effective weight here is
the raw ldsc weight.

In shared_controls mode every model in a control group uses the weight of the model without a
target column, so the shared Gram matrix is formed once and only the target's border changes;
the target coefficient follows from the Frisch-Waugh-Lovell identity.
"""

from collections import defaultdict
from typing import Sequence

import numpy as np
from attrs import frozen
from scipy.stats import norm

from mecfs_bio.constants.sldsc_constants import SLDSCWeighting

LDSC_DEFAULT_N_BLOCKS = 200
DEFAULT_CHUNK_SIZE = 64
_NULL_INTERCEPT = 1.0


@frozen(slots=True)
class SLDSCRegressionInputs:
    """
    Regression SNPs in ldsc row order (the jackknife blocks depend on it).

    chi2, n, w_ld: (S,) chi-square, integer-valued sample size, weight LD score.
    baseline_ld: (S, B) baseline LD scores; baseline_m: (B,) M_5_50.
    annotation_ld: (S, U) target/control LD scores (float32 allowed); annotation_m: (U,).
    """

    chi2: np.ndarray
    n: np.ndarray
    w_ld: np.ndarray
    baseline_ld: np.ndarray
    baseline_m: np.ndarray
    annotation_ld: np.ndarray
    annotation_m: np.ndarray

    def __attrs_post_init__(self) -> None:
        s = self.chi2.shape[0]
        for name, arr in (("chi2", self.chi2), ("n", self.n), ("w_ld", self.w_ld)):
            assert arr.ndim == 1 and arr.shape[0] == s, f"{name} must have shape ({s},), got {arr.shape}"
            assert arr.dtype.kind == "f", f"{name} must be float, got {arr.dtype}"
        for name, ld, m in (
            ("baseline", self.baseline_ld, self.baseline_m),
            ("annotation", self.annotation_ld, self.annotation_m),
        ):
            assert ld.ndim == 2 and ld.shape[0] == s, f"{name}_ld must be (S={s}, k), got {ld.shape}"
            assert ld.dtype.kind == "f", f"{name}_ld must be float, got {ld.dtype}"
            assert m.ndim == 1 and m.shape[0] == ld.shape[1], (
                f"{name}_m must have shape ({ld.shape[1]},), got {m.shape}"
            )
        assert np.all(np.isfinite(self.chi2)), "chi2 must be finite (filter before the kernel)"


@frozen(slots=True)
class SLDSCModel:
    """target_index and control_indices are columns of SLDSCRegressionInputs.annotation_ld."""

    label: str
    target_index: int
    control_indices: tuple[int, ...]

    def __attrs_post_init__(self) -> None:
        assert self.target_index not in self.control_indices, f"{self.label}: target among controls"


@frozen(slots=True)
class BatchedSLDSCResult:
    """Per-model target coefficient, jackknife SE, and one-sided p-value, in input model order."""

    labels: tuple[str, ...]
    coef: np.ndarray
    coef_se: np.ndarray
    p: np.ndarray

    def __attrs_post_init__(self) -> None:
        k = len(self.labels)
        for name, arr in (("coef", self.coef), ("coef_se", self.coef_se), ("p", self.p)):
            assert arr.ndim == 1 and arr.shape[0] == k, f"{name} must have shape ({k},), got {arr.shape}"
            assert arr.dtype.kind == "f", f"{name} must be float, got {arr.dtype}"


@frozen(slots=True)
class _SharedDesign:
    """Columns common to a control group: unscaled LD row-sum, M total, scaled design with intercept."""

    ld_sum: np.ndarray  # (S,)
    m_sum: float
    design: np.ndarray  # (S, q), last column is the intercept


@frozen(slots=True)
class _BlockNormalEquations:
    """Per-block weighted normal equations: gram (n_blocks, k, p, p), rhs (n_blocks, k, p), target first."""

    gram: np.ndarray
    rhs: np.ndarray


@frozen(slots=True)
class _TargetEstimates:
    """Target coefficient on the full data (k,) and with each block deleted (n_blocks, k), unscaled."""

    est: np.ndarray
    delete: np.ndarray


def ldsc_block_separators(n_snps: int, n_blocks: int) -> np.ndarray:
    """ldsc's Jackknife.get_separators: floor(linspace(0, n_snps, n_blocks + 1))."""
    return np.floor(np.linspace(0, n_snps, n_blocks + 1)).astype(int)


def ldsc_hsq_weights(
    ld_tot: np.ndarray, w_ld: np.ndarray, n: np.ndarray, m_tot: np.ndarray, chi2: np.ndarray
) -> np.ndarray:
    """
    ldsc's Hsq regression weights (old_weights=True path), for k models at once.

    ld_tot: (S, k) unscaled total LD per model; m_tot: (k,) total M per model.
    Returns (S, k) weights het * oc (the square-root weighting in ldsc squares back to this).
    """
    h2_agg = np.clip(
        m_tot * (chi2.mean() - _NULL_INTERCEPT) / (ld_tot * n[:, None]).mean(axis=0), 0.0, 1.0
    )
    c = h2_agg[None, :] * n[:, None] / m_tot[None, :]
    het = 1.0 / (2.0 * (_NULL_INTERCEPT + c * np.maximum(ld_tot, 1.0)) ** 2)
    return het / np.maximum(w_ld, 1.0)[:, None]


def batched_sldsc_cts(
    inputs: SLDSCRegressionInputs,
    models: Sequence[SLDSCModel],
    *,
    weighting: SLDSCWeighting,
    n_blocks: int = LDSC_DEFAULT_N_BLOCKS,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> BatchedSLDSCResult:
    assert models, "need at least one model"
    s = inputs.chi2.shape[0]
    n_blocks = min(s, n_blocks)
    separators = ldsc_block_separators(s, n_blocks)
    n_bar = float(inputs.n.mean())
    scale = inputs.n / n_bar
    est = np.empty(len(models))
    delete = np.empty((n_blocks, len(models)))
    for controls, positions in _group_by_controls(models).items():
        shared = _shared_design(inputs, controls, scale)
        targets = [models[i].target_index for i in positions]
        group = (
            _exact_group(inputs, shared, targets, scale, separators, chunk_size)
            if weighting == "exact"
            else _shared_controls_group(inputs, shared, targets, scale, separators, chunk_size)
        )
        est[positions] = group.est
        delete[:, positions] = group.delete
    pseudo = n_blocks * est[None, :] - (n_blocks - 1) * delete
    coef = est / n_bar
    coef_se = np.sqrt(np.var(pseudo, axis=0, ddof=1) / n_blocks) / n_bar
    return BatchedSLDSCResult(
        labels=tuple(m.label for m in models), coef=coef, coef_se=coef_se, p=norm.sf(coef / coef_se)
    )


def _group_by_controls(models: Sequence[SLDSCModel]) -> dict[tuple[int, ...], list[int]]:
    groups: dict[tuple[int, ...], list[int]] = defaultdict(list)
    for i, model in enumerate(models):
        groups[model.control_indices].append(i)
    return dict(groups)


def _shared_design(inputs: SLDSCRegressionInputs, controls: tuple[int, ...], scale: np.ndarray) -> _SharedDesign:
    control_ld = inputs.annotation_ld[:, list(controls)].astype(np.float64)
    unscaled = np.hstack([control_ld, inputs.baseline_ld])
    return _SharedDesign(
        ld_sum=unscaled.sum(axis=1),
        m_sum=float(inputs.annotation_m[list(controls)].sum() + inputs.baseline_m.sum()),
        design=np.hstack([unscaled * scale[:, None], np.ones((unscaled.shape[0], 1))]),
    )


def _exact_group(
    inputs: SLDSCRegressionInputs,
    shared: _SharedDesign,
    targets: list[int],
    scale: np.ndarray,
    separators: np.ndarray,
    chunk_size: int,
) -> _TargetEstimates:
    q = shared.design.shape[1]
    tri = np.triu_indices(q)
    ests, deletes = [], []
    for start in range(0, len(targets), chunk_size):
        chunk = targets[start : start + chunk_size]
        t = inputs.annotation_ld[:, chunk].astype(np.float64)
        w = ldsc_hsq_weights(
            ld_tot=shared.ld_sum[:, None] + t,
            w_ld=inputs.w_ld,
            n=inputs.n,
            m_tot=shared.m_sum + inputs.annotation_m[chunk],
            chi2=inputs.chi2,
        )
        blocks = _exact_block_normal_equations(shared.design, w, t * scale[:, None], inputs.chi2, separators, tri)
        gram_tot, rhs_tot = blocks.gram.sum(axis=0), blocks.rhs.sum(axis=0)
        est = _solve_first(gram_tot, rhs_tot)
        delete = _solve_first(gram_tot[None] - blocks.gram, rhs_tot[None] - blocks.rhs)
        ests.append(est)
        deletes.append(delete)
    return _TargetEstimates(est=np.concatenate(ests), delete=np.concatenate(deletes, axis=1))


def _exact_block_normal_equations(
    a: np.ndarray,
    w: np.ndarray,
    t: np.ndarray,
    y: np.ndarray,
    separators: np.ndarray,
    tri: tuple[np.ndarray, np.ndarray],
) -> _BlockNormalEquations:
    """
    Per-block weighted normal equations for designs [t_c, a] with per-model weights.

    a: (S, q) shared scaled design; w, t: (S, k); y: (S,).
    Returns gram (n_blocks, k, p, p) and rhs (n_blocks, k, p) with p = q + 1, target first.
    """
    q = a.shape[1]
    k = w.shape[1]
    n_blocks = len(separators) - 1
    gram = np.empty((n_blocks, k, q + 1, q + 1))
    rhs = np.empty((n_blocks, k, q + 1))
    diag = np.arange(q)
    for b in range(n_blocks):
        lo, hi = separators[b], separators[b + 1]
        a_b, w_b, t_b, y_b = a[lo:hi], w[lo:hi], t[lo:hi], y[lo:hi]
        wt = w_b * t_b
        upper = np.zeros((k, q, q))
        upper[:, tri[0], tri[1]] = ((a_b[:, tri[0]] * a_b[:, tri[1]]).T @ w_b).T
        shared = upper + np.transpose(upper, (0, 2, 1))
        shared[:, diag, diag] /= 2.0
        cross = (a_b.T @ wt).T  # (k, q)
        gram[b, :, 1:, 1:] = shared
        gram[b, :, 0, 1:] = cross
        gram[b, :, 1:, 0] = cross
        gram[b, :, 0, 0] = (wt * t_b).sum(axis=0)
        rhs[b, :, 0] = (wt * y_b[:, None]).sum(axis=0)
        rhs[b, :, 1:] = (a_b.T @ (w_b * y_b[:, None])).T
    return _BlockNormalEquations(gram=gram, rhs=rhs)


def _solve_first(gram: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """First (target) coefficient of batched solves; gram (..., p, p), rhs (..., p)."""
    return np.linalg.solve(gram, rhs[..., None])[..., 0, 0]


def _shared_controls_group(
    inputs: SLDSCRegressionInputs,
    shared: _SharedDesign,
    targets: list[int],
    scale: np.ndarray,
    separators: np.ndarray,
    chunk_size: int,
) -> _TargetEstimates:
    raise NotImplementedError("shared_controls weighting is implemented in a later step")
```


- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/sldsc/test_batched_sldsc_cts.py -v`
Expected: 3 passed. If coefficients match but SEs differ, check that ldsc's jackknife covariance is cov(pseudo, ddof=1)/n_blocks (ldsc_jackknife.py Jackknife.jknife) and that separators use linspace(0, S, ...) — NOT genomic_sem_ldsc.block_bounds, which uses a different scheme.

- [ ] **Step 5: Green and commit**

```bash
pixi r invoke green 2>&1 | tee /tmp/green_vec_sldsc.log; echo EXIT=$?
git add mecfs_bio/build_system/task/sldsc/batched_sldsc_cts.py test_mecfs_bio/unit/build_system/task/sldsc/test_batched_sldsc_cts.py
git commit -m "Add batched exact cell-type-specific S-LDSC kernel

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Alignment and VectorizedCellAnalysisBySLDSCTask (oracle: CellAnalysisByLDSCTask)

**Files:**
- Create: `mecfs_bio/build_system/task/sldsc/sldsc_regression_inputs.py`
- Create: `mecfs_bio/build_system/task/sldsc/vectorized_cell_analysis_by_sldsc_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/sldsc/test_vectorized_cell_analysis_by_sldsc_task.py`

**Interfaces:**
- Consumes: Task 1 contract, Task 2 prepare Tasks + fixture, Task 3 kernel.
- Produces: `load_sldsc_sumstats(asset: Asset, pre_pipe: DataProcessingPipe) -> pl.DataFrame` (columns SNP, Z, N); `build_regression_inputs(sumstats: pl.DataFrame, baseline: PreparedBaseline, annotation_set: PreparedAnnotationSet) -> AlignedSLDSCRegression` with fields `inputs: SLDSCRegressionInputs`, `models: tuple[SLDSCModel, ...]`; `VectorizedCellAnalysisBySLDSCTask.create(asset_id: str, source_sumstats_task: Task, baseline_task: Task, annotation_set_task: Task, weighting: SLDSCWeighting = "exact", pre_pipe: DataProcessingPipe = IdentityPipe()) -> VectorizedCellAnalysisBySLDSCTask`. Output: CSV FileAsset with columns Name, Coefficient, Coefficient_std_error, Coefficient_P_value sorted by p ascending; meta is `ResultTableMeta` identical in shape to `CellAnalysisByLDSCTask.create`.

Reproduction notes (gwaslab util_ex_ldsc._estimate_h2_cts_by_ldsc, ldsc_sumstats._read_ld_sumstats / cell_type_specific):
- Z = BETA/SE when no Z column; rsID → SNP; drop rows with null/NaN SNP, Z, EA, NEA.
- Inner join baseline (left, baseline row order) with sumstats, then inner join w_ld.
- N: gwaslab median-fills then casts to int64; here missing N on any joined row raises, then cast Float64 → Int64 (truncation).
- Duplicate SNP among joined rows raises (gwaslab would duplicate rows).
- χ² filter: keep Z² < max(0.001 · max(N), 80), N max over the joined rows.
- Annotation set must contain every regression SNP (left join, assert no nulls).

- [ ] **Step 1: Write the failing tests**

`test_vectorized_cell_analysis_by_sldsc_task.py`:

```python
from pathlib import Path

import gwaslab as gl
import numpy as np
import pandas as pd
import polars as pl
import pytest
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.gwaslab.gwaslab_cell_analysis_by_sldsc import (
    CellAnalysisByLDSCTask,
)
from mecfs_bio.build_system.task.sldsc.prepare_ldcts_annotation_set_task import (
    PrepareLdctsAnnotationSetTask,
)
from mecfs_bio.build_system.task.sldsc.prepare_sldsc_baseline_task import (
    PrepareSLDSCBaselineTask,
)
from mecfs_bio.build_system.task.sldsc.vectorized_cell_analysis_by_sldsc_task import (
    VectorizedCellAnalysisBySLDSCTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.sldsc_constants import (
    SLDSC_CTS_COEF_COL,
    SLDSC_CTS_COEF_SE_COL,
    SLDSC_CTS_NAME_COL,
    SLDSC_CTS_P_COL,
)
from test_mecfs_bio.unit.build_system.task.sldsc.sldsc_fixture import (
    FIXTURE_BASELINE_PATTERN,
    FIXTURE_LABELS,
    FIXTURE_LDCTS_FILENAME,
    FIXTURE_WEIGHTS_PATTERN,
    SLDSCFixture,
    sumstats_from_frame,
    write_sldsc_fixture,
)


def _dir_task(asset_id: str) -> FakeTask:
    return FakeTask(
        ReferenceDataDirectoryMeta(group="g", sub_group="s", sub_folder=Path("extracted"), id=AssetId(asset_id))
    )


_SUMSTATS = FakeTask(GWASLabSumStatsMeta(id=AssetId("sumstats"), trait="t", project="p"))
_BASELINE_DIR = _dir_task("baseline_dir")
_WEIGHTS_DIR = _dir_task("weights_dir")
_CTS_DIR = _dir_task("cts_dir")
_PREP_BASELINE = PrepareSLDSCBaselineTask.create(
    asset_id="prep_baseline",
    baseline_directory_task=_BASELINE_DIR,
    baseline_pattern=FIXTURE_BASELINE_PATTERN,
    weights_directory_task=_WEIGHTS_DIR,
    weights_pattern=FIXTURE_WEIGHTS_PATTERN,
)
_PREP_SET = PrepareLdctsAnnotationSetTask.create(
    asset_id="prep_set", source_directory_task=_CTS_DIR, ldcts_filename=FIXTURE_LDCTS_FILENAME
)


@frozen(slots=True)
class _Store:
    """Materializes fixture assets and runs Tasks against them, like the build system's fetch."""

    root: Path
    assets: dict[str, Asset]

    def fetch(self, asset_id: AssetId) -> Asset:
        return self.assets[str(asset_id)]

    def run(self, task: FakeTask | PrepareSLDSCBaselineTask | PrepareLdctsAnnotationSetTask | VectorizedCellAnalysisBySLDSCTask | CellAnalysisByLDSCTask) -> Asset:
        scratch = self.root / f"scratch_{task.asset_id}"
        scratch.mkdir()
        asset = task.execute(scratch_dir=scratch, fetch=self.fetch, wf=make_wf())
        self.assets[str(task.asset_id)] = asset
        return asset


def _store(tmp_path: Path, fixture: SLDSCFixture, sumstats_frame: pd.DataFrame) -> _Store:
    tmp_path.mkdir(parents=True, exist_ok=True)
    pickle_path = tmp_path / "sumstats.pickle"
    gl.dump_pickle(sumstats_from_frame(sumstats_frame), path=str(pickle_path))
    store = _Store(
        root=tmp_path,
        assets={
            "sumstats": FileAsset(pickle_path),
            "baseline_dir": DirectoryAsset(fixture.baseline_dir),
            "weights_dir": DirectoryAsset(fixture.weights_dir),
            "cts_dir": DirectoryAsset(fixture.cts_dir),
        },
    )
    store.run(_PREP_BASELINE)
    store.run(_PREP_SET)
    return store


def _vectorized(asset_id: str = "vec") -> VectorizedCellAnalysisBySLDSCTask:
    return VectorizedCellAnalysisBySLDSCTask.create(
        asset_id=asset_id, source_sumstats_task=_SUMSTATS, baseline_task=_PREP_BASELINE, annotation_set_task=_PREP_SET
    )


def _read_result(asset: Asset) -> pl.DataFrame:
    assert isinstance(asset, FileAsset)
    return pl.read_csv(asset.path).sort(SLDSC_CTS_NAME_COL)


def test_exact_mode_matches_gwaslab_cell_type_specific(tmp_path: Path) -> None:
    fixture = write_sldsc_fixture(tmp_path)
    store = _store(tmp_path, fixture, fixture.sumstats_frame)
    old = CellAnalysisByLDSCTask.create(
        asset_id="old",
        source_sumstats_task=_SUMSTATS,
        ref_ld_chr_cts_task=_CTS_DIR,
        ref_ld_chr_cts_filename=FIXTURE_LDCTS_FILENAME,
        ref_ld_chr_task=_BASELINE_DIR,
        ref_ld_chr_inner_dirname=FIXTURE_BASELINE_PATTERN,
        w_ld_chr_task=_WEIGHTS_DIR,
        w_ld_chr_inner_dirname=FIXTURE_WEIGHTS_PATTERN,
    )
    expected = _read_result(store.run(old))
    actual = _read_result(store.run(_vectorized()))
    assert actual[SLDSC_CTS_NAME_COL].to_list() == sorted(FIXTURE_LABELS)
    for col in (SLDSC_CTS_COEF_COL, SLDSC_CTS_COEF_SE_COL, SLDSC_CTS_P_COL):
        np.testing.assert_allclose(actual[col].to_numpy(), expected[col].to_numpy(), rtol=1e-5)


def test_degenerate_z_rows_do_not_change_results(tmp_path: Path) -> None:
    fixture = write_sldsc_fixture(tmp_path)
    clean = _read_result(_store(tmp_path / "clean", fixture, fixture.sumstats_frame).run(_vectorized()))
    degenerate = fixture.sumstats_frame.copy()
    degenerate.loc[[10, 11], "SE"] = 0.0  # +-inf Z: removed by the chi-square filter
    degenerate.loc[[12], ["BETA", "SE"]] = 0.0  # NaN Z: dropped
    with_bad = _read_result(_store(tmp_path / "degenerate", fixture, degenerate).run(_vectorized()))
    without_rows = fixture.sumstats_frame.drop(index=[10, 11, 12])
    dropped = _read_result(_store(tmp_path / "dropped", fixture, without_rows).run(_vectorized()))
    np.testing.assert_allclose(with_bad[SLDSC_CTS_COEF_COL].to_numpy(), dropped[SLDSC_CTS_COEF_COL].to_numpy(), rtol=1e-12)
    assert clean.height == with_bad.height


@pytest.mark.parametrize(
    "corrupt",
    [
        lambda f: f.assign(N=f["N"].astype(float).where(f.index != 500)),  # missing N on a regression SNP
        lambda f: pd.concat([f, f.iloc[[400]]], ignore_index=True),  # duplicate rsID
        lambda f: f.drop(columns="N"),  # no N column
    ],
)
def test_invalid_sumstats_fail_fast(tmp_path: Path, corrupt) -> None:  # noqa: ANN001
    fixture = write_sldsc_fixture(tmp_path)
    store = _store(tmp_path, fixture, corrupt(fixture.sumstats_frame))
    with pytest.raises(AssertionError):
        store.run(_vectorized())


def test_annotation_set_missing_regression_snps_fails(tmp_path: Path) -> None:
    fixture = write_sldsc_fixture(tmp_path)
    store = _store(tmp_path, fixture, fixture.sumstats_frame)
    small_cts = tmp_path / "small_cts"
    write_sldsc_fixture(small_cts)  # same SNP ids, then drop chromosome 3 from the cts files
    for f in (small_cts / "cts_root" / "cts").glob("*.3.l2.*"):
        f.unlink()
    store.assets["small_cts_dir"] = DirectoryAsset(small_cts / "cts_root")
    small_set = PrepareLdctsAnnotationSetTask.create(
        asset_id="small_set", source_directory_task=_dir_task("small_cts_dir"), ldcts_filename=FIXTURE_LDCTS_FILENAME
    )
    store.run(small_set)
    task = VectorizedCellAnalysisBySLDSCTask.create(
        asset_id="vec_small", source_sumstats_task=_SUMSTATS, baseline_task=_PREP_BASELINE, annotation_set_task=small_set
    )
    with pytest.raises(AssertionError):
        store.run(task)
```

Notes for the implementer:
- Type the `corrupt` parameter properly instead of `# noqa: ANN001` (ruff ANN001 is enabled): define `class _FrameCorruption(Protocol): def __call__(self, frame: pd.DataFrame) -> pd.DataFrame: ...` and annotate `corrupt: _FrameCorruption`; replace the lambdas with three small named functions.
- Type `_Store.run`'s parameter as `Task` (from base_task) rather than the union shown.
- If `gl.Sumstats` rejects a frame without N (third corruption case), construct that case with `n=None` in `sumstats_from_frame` via an optional parameter rather than skipping it.

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/sldsc/test_vectorized_cell_analysis_by_sldsc_task.py -v`
Expected: FAIL with ModuleNotFoundError for vectorized_cell_analysis_by_sldsc_task.

- [ ] **Step 3: Implement alignment**

`sldsc_regression_inputs.py`:

```python
"""
Turn gwaslab sumstats plus prepared S-LDSC references into kernel inputs, reproducing the
preprocessing of ldsc's cell_type_specific as wrapped by gwaslab: Z from BETA/SE, rows with
missing SNP/Z/alleles dropped, inner joins onto the baseline SNP order and the weight LD score,
N cast to int64 by truncation, and the chi-square filter Z^2 < max(0.001 * max(N), 80).

Deliberate divergences from gwaslab, both fail-fast: a regression SNP with missing N raises
(gwaslab fills with the median N), and a SNP appearing more than once raises (gwaslab's merge
would duplicate the row).
"""

import narwhals
import numpy as np
import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.read_spec.read_sumstats import read_sumstats
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.sldsc.batched_sldsc_cts import (
    SLDSCModel,
    SLDSCRegressionInputs,
)
from mecfs_bio.build_system.task.sldsc.prepared_sldsc_reference import (
    SNP_COL,
    W_LD_COL,
    LdScoreTable,
    PreparedAnnotationSet,
    PreparedBaseline,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
    GWASLAB_Z_COL,
)

_CHISQ_FILTER_N_COEF = 0.001
_CHISQ_FILTER_FLOOR = 80.0
_ROW_COL = "_row"


@frozen(slots=True)
class AlignedSLDSCRegression:
    inputs: SLDSCRegressionInputs
    models: tuple[SLDSCModel, ...]


def load_sldsc_sumstats(asset: Asset, pre_pipe: DataProcessingPipe) -> pl.DataFrame:
    """SNP, Z, N (Float64, possibly null) for every usable row of the sumstats."""
    data = pre_pipe.process(narwhals.from_native(read_sumstats(asset).data).lazy()).collect().to_polars()
    assert GWASLAB_SAMPLE_SIZE_COLUMN in data.columns, "sumstats need an N column for S-LDSC"
    z = pl.col(GWASLAB_Z_COL) if GWASLAB_Z_COL in data.columns else pl.col(GWASLAB_BETA_COL) / pl.col(GWASLAB_SE_COL)
    alleles = [c for c in (GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL) if c in data.columns]
    return (
        data.select(
            pl.col(GWASLAB_RSID_COL).cast(pl.String).alias(SNP_COL),
            z.cast(pl.Float64).alias(GWASLAB_Z_COL),
            pl.col(GWASLAB_SAMPLE_SIZE_COLUMN).cast(pl.Float64).fill_nan(None),
            *[pl.col(c) for c in alleles],
        )
        .drop_nulls(subset=[SNP_COL, GWASLAB_Z_COL, *alleles])
        .filter(pl.col(GWASLAB_Z_COL).is_not_nan())
        .select(SNP_COL, GWASLAB_Z_COL, GWASLAB_SAMPLE_SIZE_COLUMN)
    )


def build_regression_inputs(
    sumstats: pl.DataFrame, baseline: PreparedBaseline, annotation_set: PreparedAnnotationSet
) -> AlignedSLDSCRegression:
    regression = _filter_chisq(_cast_n(_join_baseline(sumstats, baseline)))
    baseline_rows = regression[_ROW_COL].to_numpy()
    annotation_rows = _annotation_rows(regression, annotation_set.ld)
    z = regression[GWASLAB_Z_COL].to_numpy()
    return AlignedSLDSCRegression(
        inputs=SLDSCRegressionInputs(
            chi2=z * z,
            n=regression[GWASLAB_SAMPLE_SIZE_COLUMN].to_numpy().astype(np.float64),
            w_ld=regression[W_LD_COL].to_numpy().astype(np.float64),
            baseline_ld=baseline.ld.matrix[baseline_rows].astype(np.float64),
            baseline_m=baseline.ld.m_5_50,
            annotation_ld=annotation_set.ld.matrix[annotation_rows],
            annotation_m=annotation_set.ld.m_5_50,
        ),
        models=tuple(
            SLDSCModel(
                label=m.label,
                target_index=annotation_set.ld.column_index(m.target_column),
                control_indices=tuple(annotation_set.ld.column_index(c) for c in m.control_columns),
            )
            for m in annotation_set.models
        ),
    )


def _join_baseline(sumstats: pl.DataFrame, baseline: PreparedBaseline) -> pl.DataFrame:
    joined = (
        baseline.ld.snps.select(SNP_COL)
        .with_row_index(_ROW_COL)
        .join(sumstats, on=SNP_COL, how="inner")
        .join(baseline.w_ld, on=SNP_COL, how="inner")
        .sort(_ROW_COL)
    )
    assert joined.height > 0, "no sumstats SNPs overlap the S-LDSC baseline"
    assert joined[SNP_COL].is_unique().all(), "sumstats contain duplicate rsIDs among S-LDSC regression SNPs"
    return joined


def _cast_n(joined: pl.DataFrame) -> pl.DataFrame:
    n_missing = joined[GWASLAB_SAMPLE_SIZE_COLUMN].null_count()
    assert n_missing == 0, (
        f"{n_missing} S-LDSC regression SNPs have no sample size; this indicates an upstream "
        "data-processing problem"
    )
    # gwaslab casts N with .apply("int64"): truncation toward zero.
    return joined.with_columns(pl.col(GWASLAB_SAMPLE_SIZE_COLUMN).cast(pl.Int64))


def _filter_chisq(joined: pl.DataFrame) -> pl.DataFrame:
    n_max = joined[GWASLAB_SAMPLE_SIZE_COLUMN].max()
    assert isinstance(n_max, int)
    chisq_max = max(_CHISQ_FILTER_N_COEF * n_max, _CHISQ_FILTER_FLOOR)
    return joined.filter(pl.col(GWASLAB_Z_COL) ** 2 < chisq_max)


def _annotation_rows(regression: pl.DataFrame, ld: LdScoreTable) -> np.ndarray:
    rows = (
        regression.select(SNP_COL)
        .with_row_index("_order")
        .join(ld.snps.select(SNP_COL).with_row_index(_ROW_COL), on=SNP_COL, how="left")
        .sort("_order")[_ROW_COL]
    )
    assert rows.null_count() == 0, (
        f"{rows.null_count()} regression SNPs are missing from the S-LDSC annotation set"
    )
    return rows.to_numpy()
```

- [ ] **Step 4: Implement the Task**

`vectorized_cell_analysis_by_sldsc_task.py`:

```python
"""
Cell-type-specific stratified LD score regression ("Heritability enrichment of specifically
expressed genes identifies disease-relevant tissues and cell types", Finucane et al. 2018) over
every model of a prepared annotation set at once.

weighting="exact" reproduces ldsc's cell_type_specific (as run by CellAnalysisByLDSCTask) to
floating-point tolerance; "shared_controls" is a faster approximation (see batched_sldsc_cts).
The annotation set may come from any source with an adapter to the prepared S-LDSC contract.
"""

from pathlib import Path, PurePath

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.sldsc.batched_sldsc_cts import (
    BatchedSLDSCResult,
    batched_sldsc_cts,
)
from mecfs_bio.build_system.task.sldsc.prepared_sldsc_reference import (
    read_prepared_annotation_set,
    read_prepared_baseline,
)
from mecfs_bio.build_system.task.sldsc.sldsc_regression_inputs import (
    build_regression_inputs,
    load_sldsc_sumstats,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.sldsc_constants import (
    SLDSC_CTS_COEF_COL,
    SLDSC_CTS_COEF_SE_COL,
    SLDSC_CTS_NAME_COL,
    SLDSC_CTS_P_COL,
    SLDSCWeighting,
)

_RESULT_FILENAME = "ldsc_h2_cts.csv"


@frozen(slots=True)
class VectorizedCellAnalysisBySLDSCTask(Task):
    """
    source_sumstats_task: gwaslab sumstats (pickle).
    baseline_task: prepared S-LDSC baseline (PrepareSLDSCBaselineTask).
    annotation_set_task: prepared S-LDSC annotation set (e.g. PrepareLdctsAnnotationSetTask).
    weighting: exact (default) or shared_controls.
    pre_pipe: applied to the sumstats before analysis.
    """

    meta: Meta
    source_sumstats_task: Task
    baseline_task: Task
    annotation_set_task: Task
    weighting: SLDSCWeighting = "exact"
    pre_pipe: DataProcessingPipe = IdentityPipe()

    @property
    def deps(self) -> list[Task]:
        return [self.source_sumstats_task, self.baseline_task, self.annotation_set_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        sumstats = load_sldsc_sumstats(fetch(self.source_sumstats_task.asset_id), self.pre_pipe)
        baseline = read_prepared_baseline(_directory(fetch, self.baseline_task))
        annotation_set = read_prepared_annotation_set(_directory(fetch, self.annotation_set_task))
        aligned = build_regression_inputs(sumstats, baseline, annotation_set)
        result = batched_sldsc_cts(aligned.inputs, aligned.models, weighting=self.weighting)
        out_path = scratch_dir / _RESULT_FILENAME
        _result_frame(result).write_csv(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        source_sumstats_task: Task,
        baseline_task: Task,
        annotation_set_task: Task,
        weighting: SLDSCWeighting = "exact",
        pre_pipe: DataProcessingPipe = IdentityPipe(),
    ) -> "VectorizedCellAnalysisBySLDSCTask":
        sumstats_meta = source_sumstats_task.meta
        assert isinstance(sumstats_meta, GWASLabSumStatsMeta)
        return cls(
            meta=ResultTableMeta(
                id=asset_id,
                trait=sumstats_meta.trait,
                project=sumstats_meta.project,
                sub_dir=PurePath("analysis"),
                extension=".csv",
                read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
            ),
            source_sumstats_task=source_sumstats_task,
            baseline_task=baseline_task,
            annotation_set_task=annotation_set_task,
            weighting=weighting,
            pre_pipe=pre_pipe,
        )


def _directory(fetch: Fetch, task: Task) -> Path:
    asset = fetch(task.asset_id)
    assert isinstance(asset, DirectoryAsset)
    return asset.path


def _result_frame(result: BatchedSLDSCResult) -> pl.DataFrame:
    return pl.DataFrame(
        {
            SLDSC_CTS_NAME_COL: list(result.labels),
            SLDSC_CTS_COEF_COL: result.coef,
            SLDSC_CTS_COEF_SE_COL: result.coef_se,
            SLDSC_CTS_P_COL: result.p,
        }
    ).sort(SLDSC_CTS_P_COL, maintain_order=True)
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/sldsc/test_vectorized_cell_analysis_by_sldsc_task.py -v`
Expected: all pass. If the exact oracle fails while Task 3's kernel oracle passes, the bug is in alignment: compare the regression SNP order and count against gwaslab's (run the old task with a debugger at `cell_type_specific` and dump `keep_snps`), and check the zero-variance drop and N cast.

- [ ] **Step 6: Green and commit**

```bash
pixi r invoke green 2>&1 | tee /tmp/green_vec_sldsc.log; echo EXIT=$?
git add mecfs_bio/build_system/task/sldsc/sldsc_regression_inputs.py mecfs_bio/build_system/task/sldsc/vectorized_cell_analysis_by_sldsc_task.py test_mecfs_bio/unit/build_system/task/sldsc/test_vectorized_cell_analysis_by_sldsc_task.py
git commit -m "Add VectorizedCellAnalysisBySLDSCTask reproducing gwaslab cell_type_specific

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Kernel, shared_controls mode (oracle: dense WLS with shared weights)

**Files:**
- Modify: `mecfs_bio/build_system/task/sldsc/batched_sldsc_cts.py` (replace `_shared_controls_group`)
- Test: `test_mecfs_bio/unit/build_system/task/sldsc/test_batched_sldsc_cts.py` (add tests)

**Interfaces:**
- Consumes: Task 3 internals (`_SharedDesign`, `_TargetEstimates`, `ldsc_hsq_weights`, `ldsc_block_separators`).
- Produces: `batched_sldsc_cts(..., weighting="shared_controls")` working; no signature change.

Math: per control group, w = ldsc_hsq_weights(ld_tot = shared.ld_sum[:, None], m_tot = [shared.m_sum])[:, 0]. Per block b: G_b = A_bᵀ diag(w_b) A_b, a_y,b = A_bᵀ(w_b y_b); per target chunk: AT_b = A_bᵀ(w_b ∘ T_b) (q, k), tt_b = Σ w T², ty_b = Σ w T y. For any summed (G, a_y, AT, tt, ty): solve G X = [a_y | AT]; β_t = (ty − Σ_q AT·X[:, 0]) / (tt − Σ_q AT ∘ X[:, 1:]). Delete values use totals minus block.

- [ ] **Step 1: Write the failing tests** (append to `test_batched_sldsc_cts.py`; add `ldsc_block_separators, ldsc_hsq_weights` to the kernel import)

```python
def _dense_shared_weight_oracle(
    inputs: SLDSCRegressionInputs, model: SLDSCModel, n_blocks: int
) -> tuple[float, float]:
    """Direct weighted least squares + leave-one-block-out refits, with the shared weight."""
    s = inputs.chi2.shape[0]
    controls = inputs.annotation_ld[:, list(model.control_indices)].astype(np.float64)
    shared_unscaled = np.hstack([controls, inputs.baseline_ld])
    m_shared = inputs.annotation_m[list(model.control_indices)].sum() + inputs.baseline_m.sum()
    w = ldsc_hsq_weights(
        ld_tot=shared_unscaled.sum(axis=1)[:, None], w_ld=inputs.w_ld, n=inputs.n,
        m_tot=np.array([m_shared]), chi2=inputs.chi2,
    )[:, 0]
    n_bar = inputs.n.mean()
    target = inputs.annotation_ld[:, model.target_index].astype(np.float64)
    x = np.column_stack([np.column_stack([target, shared_unscaled]) * (inputs.n / n_bar)[:, None], np.ones(s)])
    sw = np.sqrt(w)

    def fit(rows: np.ndarray) -> float:
        beta, *_ = np.linalg.lstsq(x[rows] * sw[rows, None], inputs.chi2[rows] * sw[rows], rcond=None)
        return float(beta[0])

    nb = min(s, n_blocks)
    seps = ldsc_block_separators(s, nb)
    full = fit(np.arange(s))
    delete = np.array([fit(np.r_[0 : seps[j], seps[j + 1] : s]) for j in range(nb)])
    pseudo = nb * full - (nb - 1) * delete
    return full / n_bar, float(np.sqrt(np.var(pseudo, ddof=1) / nb) / n_bar)


def test_shared_controls_matches_dense_oracle() -> None:
    inputs = _random_inputs(3_000)
    result = batched_sldsc_cts(inputs, _MODELS, weighting="shared_controls", chunk_size=2)
    for i, model in enumerate(_MODELS):
        coef, se = _dense_shared_weight_oracle(inputs, model, n_blocks=200)
        np.testing.assert_allclose(result.coef[i], coef, rtol=1e-8)
        np.testing.assert_allclose(result.coef_se[i], se, rtol=1e-8)


def test_shared_controls_close_to_exact_on_random_data() -> None:
    # Sanity check, not a requirement: shared weights should not move z-scores far.
    inputs = _random_inputs(3_000)
    exact = batched_sldsc_cts(inputs, _MODELS, weighting="exact")
    shared = batched_sldsc_cts(inputs, _MODELS, weighting="shared_controls")
    np.testing.assert_allclose(shared.coef / shared.coef_se, exact.coef / exact.coef_se, atol=0.2)
```

Replace the tuple-returning oracle helpers in this file with a small `@frozen _OracleFit(coef: float, se: float)` when implementing (applies to `_gwaslab_oracle` too).

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/sldsc/test_batched_sldsc_cts.py -v`
Expected: the two new tests FAIL with NotImplementedError; earlier tests pass.

- [ ] **Step 3: Implement shared_controls**

Replace `_shared_controls_group` in `batched_sldsc_cts.py`:

```python
def _shared_controls_group(
    inputs: SLDSCRegressionInputs,
    shared: _SharedDesign,
    targets: list[int],
    scale: np.ndarray,
    separators: np.ndarray,
    chunk_size: int,
) -> _TargetEstimates:
    w = ldsc_hsq_weights(
        ld_tot=shared.ld_sum[:, None],
        w_ld=inputs.w_ld,
        n=inputs.n,
        m_tot=np.array([shared.m_sum]),
        chi2=inputs.chi2,
    )[:, 0]
    a, y = shared.design, inputs.chi2
    blocks = [slice(separators[b], separators[b + 1]) for b in range(len(separators) - 1)]
    gram = np.stack([a[r].T @ (w[r, None] * a[r]) for r in blocks])  # (nb, q, q)
    a_y = np.stack([a[r].T @ (w[r] * y[r]) for r in blocks])  # (nb, q)
    ests, deletes = [], []
    for start in range(0, len(targets), chunk_size):
        chunk = targets[start : start + chunk_size]
        t = inputs.annotation_ld[:, chunk].astype(np.float64) * scale[:, None]
        wt = w[:, None] * t
        a_t = np.stack([a[r].T @ wt[r] for r in blocks])  # (nb, q, k)
        t_t = np.stack([(wt[r] * t[r]).sum(axis=0) for r in blocks])  # (nb, k)
        t_y = np.stack([(wt[r] * y[r, None]).sum(axis=0) for r in blocks])  # (nb, k)
        ests.append(_fwl_target_coef(gram.sum(0), a_y.sum(0), a_t.sum(0), t_t.sum(0), t_y.sum(0)))
        deletes.append(
            _fwl_target_coef(
                gram.sum(0)[None] - gram,
                a_y.sum(0)[None] - a_y,
                a_t.sum(0)[None] - a_t,
                t_t.sum(0)[None] - t_t,
                t_y.sum(0)[None] - t_y,
            )
        )
    return _TargetEstimates(est=np.concatenate(ests), delete=np.concatenate(deletes, axis=1))


def _fwl_target_coef(
    gram: np.ndarray, a_y: np.ndarray, a_t: np.ndarray, t_t: np.ndarray, t_y: np.ndarray
) -> np.ndarray:
    """
    Target coefficient of [t, A] regressions by Frisch-Waugh-Lovell, broadcasting over leading
    axes: gram (..., q, q) = A'WA, a_y (..., q) = A'Wy, a_t (..., q, k) = A'Wt, t_t and t_y (..., k).
    """
    sol = np.linalg.solve(gram, np.concatenate([a_y[..., None], a_t], axis=-1))
    numerator = t_y - (a_t * sol[..., :1]).sum(axis=-2)
    denominator = t_t - (a_t * sol[..., 1:]).sum(axis=-2)
    return numerator / denominator
```

Remove the "implemented in a later step" placeholder.

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/sldsc/test_batched_sldsc_cts.py -v`
Expected: all pass. If the sanity test's atol=0.2 fails on random data, print the max |Δz| and report it rather than loosening silently; the real-data measurement (Task 7) is what matters.

- [ ] **Step 5: Green and commit**

```bash
pixi r invoke green 2>&1 | tee /tmp/green_vec_sldsc.log; echo EXIT=$?
git add mecfs_bio/build_system/task/sldsc/batched_sldsc_cts.py test_mecfs_bio/unit/build_system/task/sldsc/test_batched_sldsc_cts.py
git commit -m "Add shared_controls fast weighting to batched S-LDSC kernel

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Generators and prepared reference asset instances

**Files:**
- Modify: `mecfs_bio/build_system/task_generator/sldsc_task_generator.py`
- Modify: `mecfs_bio/asset_generator/concrete_sldsc_generator.py`
- Create: `mecfs_bio/build_system/task_generator/vectorized_sldsc_task_generator.py`
- Create: `mecfs_bio/assets/reference_data/linkage_disequilibrium_score_reference_data/prepared/__init__.py` (empty)
- Create: `mecfs_bio/assets/reference_data/linkage_disequilibrium_score_reference_data/prepared/prepared_sldsc_reference_sets.py`
- Create: `mecfs_bio/asset_generator/concrete_vectorized_sldsc_generator.py`

**Interfaces:**
- Consumes: Tasks 2 and 4.
- Produces: `build_cell_analysis_task_group(base_name: str, entry_name: str, cell_task: Task, labels: CellOrTissueLabelRecord | None, multiple_testing_alpha: float, multiple_testing_method: Method, id_suffix: str) -> CellAnalysisTaskGroup` (free function in sldsc_task_generator.py); `PreparedAnnotationSetRecord(entry_name: str, annotation_set_task: Task, cell_or_tissue_labels_task: CellOrTissueLabelRecord | None)`; `VectorizedSLDSCTaskGenerator.create(base_name, source_sumstats_task, baseline_task, entries, multiple_testing_alpha, multiple_testing_method, weighting="exact", pre_pipe=IdentityPipe())` with `partitioned_tasks: Mapping[str, CellAnalysisTaskGroup]` and `get_terminal_tasks()`; `standard_vectorized_sldsc_task_generator(sumstats_task: Task, base_name: str, pre_pipe: DataProcessingPipe = IdentityPipe(), weighting: SLDSCWeighting = "exact") -> VectorizedSLDSCTaskGenerator`; label constants `IMMGEN_LABELS`, `MULTI_TISSUE_CHROMATIN_LABELS`, `MULTI_TISSUE_GENE_EXPR_LABELS` in concrete_sldsc_generator.py.

No unit tests for these wiring-only builders (repo convention); import-time construction plus `invoke green` (typecheck, import-linter) is the check.

- [ ] **Step 1: Refactor sldsc_task_generator.py**

Change `CellAnalysisTaskGroup.cell_analysis_task: CellAnalysisByLDSCTask` to `cell_analysis_task: Task`. Move the body of the per-entry loop after `cell_task = ...` into:

```python
VEC_ID_SUFFIX = "_vec"


def build_cell_analysis_task_group(
    base_name: str,
    entry_name: str,
    cell_task: Task,
    labels: CellOrTissueLabelRecord | None,
    multiple_testing_alpha: float,
    multiple_testing_method: Method,
    id_suffix: str,
) -> CellAnalysisTaskGroup:
    """
    Downstream tasks for one S-LDSC result table: multiple-testing correction, a markdown table,
    and (when labels are given) a label join and scatter plot. id_suffix is appended to every
    generated asset id.
    """
    prefix = base_name + "_" + entry_name
    multiple_testing_task = MultipleTestingTableTask.create_from_result_table_task(
        asset_id=prefix + "_s_ldsc_multiple_testing_correction" + id_suffix,
        alpha=multiple_testing_alpha,
        p_value_column=SLDSC_CTS_P_COL,
        source_task=cell_task,
        apply_filter=False,
        method=multiple_testing_method,
    )
    markdown_task = ConvertDataFrameToMarkdownTask.create_from_result_table_task(
        source_task=multiple_testing_task,
        asset_id=prefix + "_s_ldsc_cell_analysis_md_table" + id_suffix,
        pipe=DropColPipe(["_Corrected P Value_", SLDSC_CTS_COEF_SE_COL]),
    )
    if labels is None:
        add_labels_task, plot_task = None, None
    else:
        add_labels_task = JoinDataFramesTask.create_from_result_df(
            asset_id=prefix + "_add_labels" + id_suffix,
            result_df_task=multiple_testing_task,
            reference_df_task=labels.cell_or_tissue_label_task,
            how="left",
            left_on=[labels.left_join_on],
            right_on=[labels.right_join_on],
            df_1_pipe=labels.pipe_left,
            df_2_pipe=labels.pipe_right,
        )
        plot_task = SLDSCScatterPlotTask.create(
            asset_id=prefix + "_cell_analysis_s_ldsc_plot" + id_suffix,
            source_task=add_labels_task,
        )
    return CellAnalysisTaskGroup(
        cell_analysis_task=cell_task,
        multiple_testing_task=multiple_testing_task,
        add_categories_task=add_labels_task,
        plot_task=plot_task,
        multiple_testing_task_markdown=markdown_task,
    )
```

Copy the exact argument values (DropColPipe columns, join settings) from the current loop body; the snippet above must produce identical tasks with `id_suffix=""`. `SLDSCTaskGenerator.create` then builds `cell_task` as today and calls `build_cell_analysis_task_group(..., id_suffix="")`. Import `SLDSC_CTS_P_COL`, `SLDSC_CTS_COEF_SE_COL` from `mecfs_bio.constants.sldsc_constants`.

- [ ] **Step 2: Verify old asset ids are unchanged**

Run a throwaway check (not committed):

```bash
git stash
pixi r python -c "from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_sldsc import DECODE_ME_S_LDSC as G; print(sorted(str(t.asset_id) for t in G.get_terminal_tasks()))" | tail -1 > /tmp/ids_before.txt
git stash pop
pixi r python -c "from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_sldsc import DECODE_ME_S_LDSC as G; print(sorted(str(t.asset_id) for t in G.get_terminal_tasks()))" | tail -1 > /tmp/ids_after.txt
diff /tmp/ids_before.txt /tmp/ids_after.txt && echo IDS_UNCHANGED
```

Expected: `IDS_UNCHANGED`.

- [ ] **Step 3: Hoist label records in concrete_sldsc_generator.py**

Move the three inline `CellOrTissueLabelRecord(...)` values (ImmGen, Multi_tissue_chromatin, Multi_tissue_gene_expr) to module-level constants `IMMGEN_LABELS`, `MULTI_TISSUE_CHROMATIN_LABELS`, `MULTI_TISSUE_GENE_EXPR_LABELS` with unchanged contents, and reference them in `standard_sldsc_task_generator`. Re-run Step 2's check; expected `IDS_UNCHANGED`.

- [ ] **Step 4: Create the prepared reference asset instances**

`prepared_sldsc_reference_sets.py`:

```python
"""
Prepared (parquet) forms of the standard S-LDSC reference data, for the vectorized
cell-type-specific S-LDSC Task: the six LDSC-SEG annotation sets of Finucane et al. 2018 and the
baseline v1.2 LD scores with HapMap3 (no MHC) regression weights.
"""

from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_cahoy_ld_scores_extracted import (
    PARTITIONED_MODEL_CAHOY_LD_SCORES_EXTRACTED_RD,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_corces_atac_ld_scores_extracted import (
    PARTITIONED_MODEL_CORCES_ATAC_LD_SCORES_EXTRACTED_RD,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_gtex_brain_ld_scores_extracted import (
    PARTITIONED_MODEL_GTEX_BRAIN_LD_SCORES_EXTRACTED_RD,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_immgen_ld_scores_extracted import (
    PARTITIONED_MODEL_IMMGEN_LD_SCORES_EXTRACTED_RD,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_multi_tissue_chromatin_ld_scores_extracted import (
    PARTITIONED_MODEL_MULTI_TISSUE_CHROMATIN_LD_SCORES_EXTRACTED_RD,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_multi_tissue_gene_expr_ld_score_extracted import (
    PARTITIONED_MODEL_MULTI_TISSUE_GENE_EXPR_LD_SCORES_EXTRACTED_RD,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_regression_weights_extracted import (
    PARTITIONED_MODEL_REGRESSION_WEIGHTS_EXTRACTED,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.thousand_genome_baseline_ld_extracted import (
    BASE_MODEL_PARTITIONED_LD_SCORES_EXTRACTED,
)
from mecfs_bio.build_system.task.sldsc.prepare_ldcts_annotation_set_task import (
    PrepareLdctsAnnotationSetTask,
)
from mecfs_bio.build_system.task.sldsc.prepare_sldsc_baseline_task import (
    PrepareSLDSCBaselineTask,
)

PREPARED_BASELINE_V1_2_HM3_WEIGHTS = PrepareSLDSCBaselineTask.create(
    asset_id="prepared_sldsc_baseline_v1_2_hm3_no_mhc_weights",
    baseline_directory_task=BASE_MODEL_PARTITIONED_LD_SCORES_EXTRACTED,
    baseline_pattern="baseline_v1.2/baseline.@",
    weights_directory_task=PARTITIONED_MODEL_REGRESSION_WEIGHTS_EXTRACTED,
    weights_pattern="1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.@",
)

PREPARED_CAHOY_ANNOTATION_SET = PrepareLdctsAnnotationSetTask.create(
    asset_id="prepared_sldsc_cahoy_annotation_set",
    source_directory_task=PARTITIONED_MODEL_CAHOY_LD_SCORES_EXTRACTED_RD,
    ldcts_filename="Cahoy.ldcts",
)
PREPARED_GTEX_BRAIN_ANNOTATION_SET = PrepareLdctsAnnotationSetTask.create(
    asset_id="prepared_sldsc_gtex_brain_annotation_set",
    source_directory_task=PARTITIONED_MODEL_GTEX_BRAIN_LD_SCORES_EXTRACTED_RD,
    ldcts_filename="GTEx_brain.ldcts",
)
PREPARED_IMMGEN_ANNOTATION_SET = PrepareLdctsAnnotationSetTask.create(
    asset_id="prepared_sldsc_immgen_annotation_set",
    source_directory_task=PARTITIONED_MODEL_IMMGEN_LD_SCORES_EXTRACTED_RD,
    ldcts_filename="ImmGen.ldcts",
)
PREPARED_CORCES_ATAC_ANNOTATION_SET = PrepareLdctsAnnotationSetTask.create(
    asset_id="prepared_sldsc_corces_atac_annotation_set",
    source_directory_task=PARTITIONED_MODEL_CORCES_ATAC_LD_SCORES_EXTRACTED_RD,
    ldcts_filename="Corces_ATAC.ldcts",
)
PREPARED_MULTI_TISSUE_CHROMATIN_ANNOTATION_SET = PrepareLdctsAnnotationSetTask.create(
    asset_id="prepared_sldsc_multi_tissue_chromatin_annotation_set",
    source_directory_task=PARTITIONED_MODEL_MULTI_TISSUE_CHROMATIN_LD_SCORES_EXTRACTED_RD,
    ldcts_filename="Multi_tissue_chromatin.ldcts",
)
PREPARED_MULTI_TISSUE_GENE_EXPR_ANNOTATION_SET = PrepareLdctsAnnotationSetTask.create(
    asset_id="prepared_sldsc_multi_tissue_gene_expr_annotation_set",
    source_directory_task=PARTITIONED_MODEL_MULTI_TISSUE_GENE_EXPR_LD_SCORES_EXTRACTED_RD,
    ldcts_filename="Multi_tissue_gene_expr.ldcts",
)
```

- [ ] **Step 5: Create the vectorized generator**

`vectorized_sldsc_task_generator.py`:

```python
"""
Task generator applying vectorized cell-type-specific S-LDSC to one GWAS across several prepared
annotation sets, with the same downstream tables and plots as SLDSCTaskGenerator. Agnostic to the
source of the annotation sets; every generated asset id carries the _vec suffix.
"""

from typing import Mapping, Sequence

from attrs import frozen

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.multiple_testing_table_task import Method
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.sldsc.vectorized_cell_analysis_by_sldsc_task import (
    VectorizedCellAnalysisBySLDSCTask,
)
from mecfs_bio.build_system.task_generator.sldsc_task_generator import (
    VEC_ID_SUFFIX,
    CellAnalysisTaskGroup,
    CellOrTissueLabelRecord,
    build_cell_analysis_task_group,
)
from mecfs_bio.constants.sldsc_constants import SLDSCWeighting


@frozen(slots=True)
class PreparedAnnotationSetRecord:
    entry_name: str
    annotation_set_task: Task
    cell_or_tissue_labels_task: CellOrTissueLabelRecord | None


@frozen(slots=True)
class VectorizedSLDSCTaskGenerator:
    partitioned_tasks: Mapping[str, CellAnalysisTaskGroup]

    def get_terminal_tasks(self) -> Sequence[Task]:
        result: list[Task] = []
        for group in self.partitioned_tasks.values():
            result.extend(group.terminal_tasks())
        return result

    @classmethod
    def create(
        cls,
        base_name: str,
        source_sumstats_task: Task,
        baseline_task: Task,
        entries: Sequence[PreparedAnnotationSetRecord],
        multiple_testing_alpha: float,
        multiple_testing_method: Method,
        weighting: SLDSCWeighting = "exact",
        pre_pipe: DataProcessingPipe = IdentityPipe(),
    ) -> "VectorizedSLDSCTaskGenerator":
        return cls(
            {
                entry.entry_name: build_cell_analysis_task_group(
                    base_name=base_name,
                    entry_name=entry.entry_name,
                    cell_task=VectorizedCellAnalysisBySLDSCTask.create(
                        asset_id=base_name + "_" + entry.entry_name + "_cell_analysis_by_s_ldsc" + VEC_ID_SUFFIX,
                        source_sumstats_task=source_sumstats_task,
                        baseline_task=baseline_task,
                        annotation_set_task=entry.annotation_set_task,
                        weighting=weighting,
                        pre_pipe=pre_pipe,
                    ),
                    labels=entry.cell_or_tissue_labels_task,
                    multiple_testing_alpha=multiple_testing_alpha,
                    multiple_testing_method=multiple_testing_method,
                    id_suffix=VEC_ID_SUFFIX,
                )
                for entry in entries
            }
        )
```

`concrete_vectorized_sldsc_generator.py`:

```python
"""
Asset generator applying vectorized S-LDSC with the six standard LDSC-SEG annotation sets
(Finucane et al. 2018), mirroring standard_sldsc_task_generator.
"""

from mecfs_bio.asset_generator.concrete_sldsc_generator import (
    IMMGEN_LABELS,
    MULTI_TISSUE_CHROMATIN_LABELS,
    MULTI_TISSUE_GENE_EXPR_LABELS,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.prepared.prepared_sldsc_reference_sets import (
    PREPARED_BASELINE_V1_2_HM3_WEIGHTS,
    PREPARED_CAHOY_ANNOTATION_SET,
    PREPARED_CORCES_ATAC_ANNOTATION_SET,
    PREPARED_GTEX_BRAIN_ANNOTATION_SET,
    PREPARED_IMMGEN_ANNOTATION_SET,
    PREPARED_MULTI_TISSUE_CHROMATIN_ANNOTATION_SET,
    PREPARED_MULTI_TISSUE_GENE_EXPR_ANNOTATION_SET,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task_generator.vectorized_sldsc_task_generator import (
    PreparedAnnotationSetRecord,
    VectorizedSLDSCTaskGenerator,
)
from mecfs_bio.constants.sldsc_constants import SLDSCWeighting


def standard_vectorized_sldsc_task_generator(
    sumstats_task: Task,
    base_name: str,
    pre_pipe: DataProcessingPipe = IdentityPipe(),
    weighting: SLDSCWeighting = "exact",
) -> VectorizedSLDSCTaskGenerator:
    return VectorizedSLDSCTaskGenerator.create(
        base_name=base_name,
        source_sumstats_task=sumstats_task,
        baseline_task=PREPARED_BASELINE_V1_2_HM3_WEIGHTS,
        entries=[
            PreparedAnnotationSetRecord("cahoy_cns", PREPARED_CAHOY_ANNOTATION_SET, None),
            PreparedAnnotationSetRecord("gtex_brain", PREPARED_GTEX_BRAIN_ANNOTATION_SET, None),
            PreparedAnnotationSetRecord("immgen", PREPARED_IMMGEN_ANNOTATION_SET, IMMGEN_LABELS),
            PreparedAnnotationSetRecord("corces_atac", PREPARED_CORCES_ATAC_ANNOTATION_SET, None),
            PreparedAnnotationSetRecord(
                "multi_tissue_chromatin", PREPARED_MULTI_TISSUE_CHROMATIN_ANNOTATION_SET, MULTI_TISSUE_CHROMATIN_LABELS
            ),
            PreparedAnnotationSetRecord(
                "multi_tissue_gene_expression", PREPARED_MULTI_TISSUE_GENE_EXPR_ANNOTATION_SET, MULTI_TISSUE_GENE_EXPR_LABELS
            ),
        ],
        multiple_testing_alpha=0.01,
        multiple_testing_method="fdr_bh",
        weighting=weighting,
        pre_pipe=pre_pipe,
    )
```

Entry names, label records, alpha and method must match `standard_sldsc_task_generator` exactly; check against it while implementing.

- [ ] **Step 6: Import check, green, commit**

```bash
pixi r python -c "from mecfs_bio.asset_generator.concrete_vectorized_sldsc_generator import standard_vectorized_sldsc_task_generator as g; from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_sldsc import DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS as s; print(len(g(s, 'probe').get_terminal_tasks()))"
pixi r invoke green 2>&1 | tee /tmp/green_vec_sldsc.log; echo EXIT=$?
git add mecfs_bio/build_system/task_generator/sldsc_task_generator.py mecfs_bio/asset_generator/concrete_sldsc_generator.py mecfs_bio/build_system/task_generator/vectorized_sldsc_task_generator.py mecfs_bio/assets/reference_data/linkage_disequilibrium_score_reference_data/prepared/__init__.py mecfs_bio/assets/reference_data/linkage_disequilibrium_score_reference_data/prepared/prepared_sldsc_reference_sets.py mecfs_bio/asset_generator/concrete_vectorized_sldsc_generator.py
git commit -m "Add vectorized S-LDSC generators and prepared reference assets

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

If `DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS` is not importable from that module, import it from wherever decode_me_sldsc.py imports it. Expected print: the same count as the old generator's terminal tasks (labelled sets contribute a plot as well as a markdown table).

---

### Task 7: Real-data experiment (profiling + DecodeME validation)

**Files:**
- Create: `experiments/claude/vectorized_sldsc/vectorized_sldsc_experiment.py`

**Interfaces:**
- Consumes: everything above; `DEFAULT_RUNNER` (`mecfs_bio/analysis/runner/default_runner.py`), `DECODE_ME_S_LDSC` (old results), `DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS`, `ComputePIfNeededPipe`.
- Produces: logs (uncommitted) and a summary reported back to the user.

- [ ] **Step 1: Write the script**

```python
"""
Vectorized S-LDSC validation on DecodeME.

Modes (one per process so peak RSS is per configuration):
    profile-old   cProfile the existing CellAnalysisByLDSCTask on GTEx_brain (13 models).
    run-exact     build prepared references (first run only) and vectorized exact results.
    run-shared    vectorized shared_controls results.
    compare       compare old vs exact vs shared for all six sets.

Usage:
    pixi r python experiments/claude/vectorized_sldsc/vectorized_sldsc_experiment.py <mode> 2>&1 | tee experiments/claude/vectorized_sldsc/<mode>.log
"""

import cProfile
import pstats
import resource
import sys
import time
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.asset_generator.concrete_vectorized_sldsc_generator import (
    standard_vectorized_sldsc_task_generator,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_sldsc import (
    DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
    DECODE_ME_S_LDSC,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.pipes.compute_p_if_needed_pipe import ComputePIfNeededPipe
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.sldsc_constants import (
    SLDSC_CTS_COEF_COL,
    SLDSC_CTS_COEF_SE_COL,
    SLDSC_CTS_NAME_COL,
    SLDSC_CTS_P_COL,
)

_FDR_ALPHA = 0.01
_EXACT = standard_vectorized_sldsc_task_generator(
    DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS, "decode_me_gwas_1", pre_pipe=ComputePIfNeededPipe()
)
_SHARED = standard_vectorized_sldsc_task_generator(
    DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
    "decode_me_gwas_1_shared_controls",
    pre_pipe=ComputePIfNeededPipe(),
    weighting="shared_controls",
)


def _peak_rss_gb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2


def _cell_tasks(generator) -> dict[str, object]:  # noqa: ANN001 - replace with a Protocol when implementing
    return {name: group.cell_analysis_task for name, group in generator.partitioned_tasks.items()}


def profile_old() -> None:
    task = DECODE_ME_S_LDSC.partitioned_tasks["gtex_brain"].cell_analysis_task
    assets = DEFAULT_RUNNER.run(list(task.deps))
    with TemporaryDirectory() as tmp:
        profiler = cProfile.Profile()
        start = time.perf_counter()
        profiler.runcall(task.execute, scratch_dir=Path(tmp), fetch=assets.__getitem__, wf=make_wf())
        print(f"old GTEx_brain wall time: {time.perf_counter() - start:.1f}s")
    pstats.Stats(profiler).sort_stats("cumulative").print_stats(30)
    print(f"peak RSS {_peak_rss_gb():.2f} GB")


def run(generator, label: str) -> None:  # noqa: ANN001
    for name, task in _cell_tasks(generator).items():
        start = time.perf_counter()
        DEFAULT_RUNNER.run([task])
        print(f"{label} {name}: {time.perf_counter() - start:.1f}s (includes prepared-reference build on first run)")
    print(f"{label} peak RSS {_peak_rss_gb():.2f} GB")


def _bh_reject(p: np.ndarray, alpha: float) -> np.ndarray:
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    passed = np.zeros(len(p), dtype=bool)
    below = np.nonzero(ranked <= alpha)[0]
    if below.size:
        passed[order[: below.max() + 1]] = True
    return passed


def compare() -> None:
    old_tasks = {n: g.cell_analysis_task for n, g in DECODE_ME_S_LDSC.partitioned_tasks.items()}
    assets = DEFAULT_RUNNER.run([*old_tasks.values(), *_cell_tasks(_EXACT).values(), *_cell_tasks(_SHARED).values()])

    def read(task) -> pl.DataFrame:  # noqa: ANN001
        asset = assets[task.asset_id]
        assert isinstance(asset, FileAsset)
        return pl.read_csv(asset.path).sort(SLDSC_CTS_NAME_COL)

    for name, old_task in old_tasks.items():
        old = read(old_task)
        z_old = (old[SLDSC_CTS_COEF_COL] / old[SLDSC_CTS_COEF_SE_COL]).to_numpy()
        for label, generator in (("exact", _EXACT), ("shared", _SHARED)):
            new = read(_cell_tasks(generator)[name])
            assert new[SLDSC_CTS_NAME_COL].to_list() == old[SLDSC_CTS_NAME_COL].to_list()
            z_new = (new[SLDSC_CTS_COEF_COL] / new[SLDSC_CTS_COEF_SE_COL]).to_numpy()
            dz = np.abs(z_new - z_old)
            dcoef = np.abs(new[SLDSC_CTS_COEF_COL].to_numpy() - old[SLDSC_CTS_COEF_COL].to_numpy()) / old[SLDSC_CTS_COEF_SE_COL].to_numpy()
            fdr_changes = int(np.sum(_bh_reject(new[SLDSC_CTS_P_COL].to_numpy(), _FDR_ALPHA) != _bh_reject(old[SLDSC_CTS_P_COL].to_numpy(), _FDR_ALPHA)))
            print(f"{name:32s} {label:6s} n={old.height:4d} max|dz|={dz.max():.2e} median|dz|={np.median(dz):.2e} max|dcoef|/se={dcoef.max():.2e} FDR changes={fdr_changes}")


if __name__ == "__main__":
    mode = sys.argv[1]
    {"profile-old": profile_old, "run-exact": lambda: run(_EXACT, "exact"), "run-shared": lambda: run(_SHARED, "shared"), "compare": compare}[mode]()
```

Implementation notes: replace the `# noqa: ANN001` parameters with real types (`VectorizedSLDSCTaskGenerator`, `Task`); `run` and the generator objects are typed differently for old vs new, so give `_cell_tasks` a small Protocol with `partitioned_tasks: Mapping[str, CellAnalysisTaskGroup]`. If `DEFAULT_RUNNER.run` returns assets keyed differently than `task.asset_id`, adapt `read`. Locate `ComputePIfNeededPipe`'s module from decode_me_sldsc.py's imports. The old results for DecodeME should already be cached; if the runner wants to rebuild them, stop and ask the user before spending hours.

- [ ] **Step 2: Run the modes**

```bash
mkdir -p experiments/claude/vectorized_sldsc
for mode in profile-old run-exact run-shared compare; do
  pixi r python experiments/claude/vectorized_sldsc/vectorized_sldsc_experiment.py $mode 2>&1 | tee experiments/claude/vectorized_sldsc/$mode.log
done
```

Run `profile-old` in the background (it may take tens of minutes for 13 models); the others should take minutes once prepared references exist.

- [ ] **Step 3: Evaluate against the success criteria**

- exact: max |dz| < 1e-4 and FDR changes = 0 for all six sets.
- run-exact per-set times in minutes; peak RSS well below 16 GB.
- profile-old: confirm whether most cumulative time sits in ldsc_parse reading (`_read_chr_split_files`, `read_csv`, merges) versus `Hsq`.
- shared: report max/median |dz| and FDR changes; no pass/fail.

If exact fails on real data while unit tests pass, likely causes are row order (pandas vs polars sort ties at identical CHR/BP) or SNP de-duplication; investigate with the superpowers:systematic-debugging skill before changing code.

- [ ] **Step 4: Commit the script only (never the logs)**

```bash
git add experiments/claude/vectorized_sldsc/vectorized_sldsc_experiment.py
git commit -m "Add vectorized S-LDSC DecodeME validation experiment

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

Report the comparison table, timings, peak RSS and profiling conclusion to the user.
