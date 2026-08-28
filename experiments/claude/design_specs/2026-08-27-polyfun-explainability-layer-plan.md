# Polyfun Explainability Layer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the layer that explains a polyfun-vs-uniform SUSIE fine-mapping result in annotation terms: matched run pairs, prior-lift/contrast tables, a docs-facing display table, an explainability figure (PNG+SVG), and the inner/outer asset generators.

**Architecture:** Everything consumes the existing `SusieRFinemapTask` directory outputs plus the two Spec 1 assets (ridge `weights.parquet` and the 187-annotation parquet). One additive change to `SusieRFinemapTask` (a fail-fast prior-coverage guard + a `prior.parquet` output); then two new tasks (contrast, plot) in a new subpackage; then a new inner+outer asset generator. No existing plot task or generator is touched.

**Tech Stack:** Python, polars, numpy, matplotlib, rpy2/susieR (existing), the repo build system (`Task`, `FileAsset`/`DirectoryAsset`, `Fetch`, `create()` meta-derivation), pixi (`pixi r`), pytest with injected-`fetch` + `FakeTask` unit tests.

**Spec:** `experiments/claude/design_specs/2026-08-27-polyfun-explainability-layer-design.md`

## Global Constraints

- Run everything via pixi: `pixi r <command>`, `pixi r python <script>`, `pixi r invoke green`.
- After each task, run `pixi r invoke green` (ruff lint/format, typos, lychee, import-linter, ty typecheck, actionlint, pytest-testmon). Capture to a logfile (testmon can report "no tests ran" after env-only changes; tail shows only R noise).
- Docstrings: no backticks around inline code; no RST.
- Prefer polars over pandas for new dataframe code. Use `Path`/`PurePath` for filesystem paths. Type enumerable string params as `Literal`. Column names come from constants, never repeated literals.
- Tests are Task-level (the Task is the public API), dependency-injected (no monkeypatch, no skipif); do not assert on log/error-message text; share one constant across creation and assertion.
- Existing constant values relied on (copy verbatim):
  - gwaslab columns: `GWASLAB_CHROM_COL="CHR"`, `GWASLAB_POS_COL="POS"`, `GWASLAB_EFFECT_ALLELE_COL="EA"`, `GWASLAB_NON_EFFECT_ALLELE_COL="NEA"`, `GWASLAB_BETA_COL="BETA"`, `GWASLAB_SE_COL="SE"` (in `mecfs_bio/constants/gwaslab_constants.py`).
  - annotation parquet key columns: `ANNOT_KEY_COLUMNS = ["CHR", "BP", "SNP", "CM"]` (note the annotation position column is `BP`, not `POS`).
  - SUSIE run output filenames/columns (in `mecfs_bio/build_system/task/r_tasks/susie_r_finemap_task.py`): `PIP_FILENAME="pip.parquet"`, `PIP_COLUMN="PIP"`, `FILTERED_GWAS_FILENAME="filtered_gwas.parquet"`, `COMBINED_CS_FILENAME="combined_cs.parquet"`, `CS_COLUMN="cs"` (values like "L1"), `FILTERED_LD_FILENAME="filtered_ld.npy"`.
  - ridge weights (in `ridge_annotation_weights_task.py`): `WEIGHTS_PARQUET_FILENAME="weights.parquet"`, `ANNOTATION_COL="annotation"`, `GAMMA_RAW_COL="gamma_raw"`, `FAMILY_COL="family"`.
  - prior asset (in `polyfun_precomputed_prior.py`): `POLYFUN_PRIOR_COL="prior"`, `create_prior_col_pipe(q)`, `COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS`.
  - assets (in `mecfs_bio/assets/reference_data/polyfun/annotations/`): `BASELINE_LF_ANNOTATION_MATRIX` (annotation parquet task), `BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS` (ridge weights task).

---

## File Structure

- Modify: `mecfs_bio/build_system/task/r_tasks/susie_r_finemap_task.py` — add the prior-coverage guard in `align_data` and a `prior.parquet` output (Task 1).
- Modify: `test_mecfs_bio/unit/build_system/task/test_susie_r_finemap_task.py` — extend for the guard + `prior.parquet` (Task 1).
- Modify: `mecfs_bio/constants/polyfun_annotation_families.py` — add `FAMILY_SHORT_LABELS` (Task 2).
- Create: `mecfs_bio/build_system/task/polyfun_explain/__init__.py` — empty (Task 2).
- Create: `mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_contrast_task.py` — `PolyfunExplainContrastTask` (Task 2).
- Create: `test_mecfs_bio/unit/build_system/task/polyfun_explain/__init__.py` + `test_polyfun_explain_contrast_task.py` (Task 2).
- Create: `mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_plot_task.py` — `PolyfunExplainPlotTask` (Task 3).
- Create: `test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_plot_task.py` (Task 3).
- Create: `mecfs_bio/asset_generator/polyfun_explain_fine_mapping_asset_generator.py` — inner + outer generators (Task 4).
- Create: `test_mecfs_bio/unit/asset_generator/test_polyfun_explain_fine_mapping_asset_generator.py` (Task 4).
- Create: one DecodeME demonstrator analysis module wiring the outer generator (Task 5).

---

## Task 1: SusieRFinemapTask — prior coverage guard + prior.parquet

**Files:**
- Modify: `mecfs_bio/build_system/task/r_tasks/susie_r_finemap_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/test_susie_r_finemap_task.py`

**Interfaces:**
- Consumes: nothing new.
- Produces: module constants `PRIOR_FILENAME = "prior.parquet"` and `PRIOR_WEIGHT_COLUMN = "prior_weight"`. A new `prior.parquet` in every SUSIE run directory: columns `CHR, POS, EA, NEA, prior_weight` (Int64 ok for CHR/POS here; this file is machine-read, not docs-rendered), one row per retained variant in `filtered_gwas.parquet` order. `align_data` now raises `ValueError` when the polyfun prior does not cover every (gwas ∩ ld) variant.

- [ ] **Step 1: Write the failing tests**

Add to `test_mecfs_bio/unit/build_system/task/test_susie_r_finemap_task.py`. Reuse the existing `susie_prerequisite_file_tasks` and `dummy_prior_task` fixtures. Add a new fixture for a prior missing one variant, and import the new constants.

```python
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    PRIOR_FILENAME,
    PRIOR_WEIGHT_COLUMN,
)


@pytest.fixture
def dummy_prior_task_missing_one(tmp_path: Path) -> Iterator[Task]:
    """Prior table covering only 99 of the 100 synthetic variants (POS 0 omitted)."""
    m = 100
    prior_data = pd.DataFrame(
        {
            "CHR": [1] * (m - 1),
            "BP": list(range(1, m)),  # omits BP == 0
            "A1": "A",
            "A2": "C",
            "snpvar": np.linspace(1.0, 2.0, m - 1),
        }
    )
    prior_path = tmp_path / "prior_data_missing"
    prior_data.to_parquet(prior_path)
    yield ExternalFileCopyTask(
        SimpleFileMeta(
            AssetId("prior_data_missing"),
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        ),
        external_path=prior_path,
    )


def _run_susie_task_to_store(tmp_path: Path, susie_tsk: SusieRFinemapTask):
    tasks = find_tasks([susie_tsk])
    wf = make_wf()
    info = VerifyingTraceInfo.empty()
    asset_dir = tmp_path / "asset_dir"
    asset_dir.mkdir(exist_ok=True, parents=True)
    meta_to_path = SimpleMetaToPath(root=asset_dir)
    rebuilder = VerifyingTraceRebuilder(SimpleHasher.md5_hasher())
    store, _ = topological(
        rebuilder=rebuilder,
        tasks=tasks,
        targets=[susie_tsk.asset_id],
        wf=wf,
        info=info,
        meta_to_path=meta_to_path,
    )
    return store


def test_prior_coverage_guard_raises_on_gap(
    tmp_path: Path,
    susie_prerequisite_file_tasks: tuple[Task, Task, Task, list[int]],
    dummy_prior_task_missing_one: Task,
):
    gwas_data_task, ld_labels_task, ld_matrix_task, _ = susie_prerequisite_file_tasks
    susie_tsk = SusieRFinemapTask(
        meta=SimpleDirectoryMeta(AssetId("directory")),
        gwas_data_task=gwas_data_task,
        ld_labels_task=ld_labels_task,
        ld_matrix_source=BroadInstituteFormatLDMatrix(ld_matrix_task),
        effective_sample_size=_susie_n,
        max_credible_sets=10,
        prior_info=PriorInfo(
            prior_task=dummy_prior_task_missing_one, prior_col="snpvar"
        ),
    )
    with pytest.raises(ValueError):
        _run_susie_task_to_store(tmp_path, susie_tsk)


def test_prior_parquet_written_for_polyfun_run(
    tmp_path: Path,
    susie_prerequisite_file_tasks: tuple[Task, Task, Task, list[int]],
    dummy_prior_task: Task,
):
    gwas_data_task, ld_labels_task, ld_matrix_task, _ = susie_prerequisite_file_tasks
    susie_tsk = SusieRFinemapTask(
        meta=SimpleDirectoryMeta(AssetId("directory")),
        gwas_data_task=gwas_data_task,
        ld_labels_task=ld_labels_task,
        ld_matrix_source=BroadInstituteFormatLDMatrix(ld_matrix_task),
        effective_sample_size=_susie_n,
        max_credible_sets=10,
        prior_info=PriorInfo(prior_task=dummy_prior_task, prior_col="snpvar"),
    )
    store = _run_susie_task_to_store(tmp_path, susie_tsk)
    asset = store[susie_tsk.asset_id]
    assert isinstance(asset, DirectoryAsset)
    prior = pl.read_parquet(asset.path / PRIOR_FILENAME)
    assert prior.height == 100
    assert set(["CHR", "POS", "EA", "NEA", PRIOR_WEIGHT_COLUMN]).issubset(
        prior.columns
    )
    weights = prior[PRIOR_WEIGHT_COLUMN].to_numpy()
    assert abs(weights.min() - 1.0) < 1e-6
    assert abs(weights.max() - 2.0) < 1e-6


def test_prior_parquet_is_constant_for_uniform_run(
    tmp_path: Path,
    susie_prerequisite_file_tasks: tuple[Task, Task, Task, list[int]],
):
    gwas_data_task, ld_labels_task, ld_matrix_task, _ = susie_prerequisite_file_tasks
    susie_tsk = SusieRFinemapTask(
        meta=SimpleDirectoryMeta(AssetId("directory")),
        gwas_data_task=gwas_data_task,
        ld_labels_task=ld_labels_task,
        ld_matrix_source=BroadInstituteFormatLDMatrix(ld_matrix_task),
        effective_sample_size=_susie_n,
        max_credible_sets=10,
    )
    store = _run_susie_task_to_store(tmp_path, susie_tsk)
    asset = store[susie_tsk.asset_id]
    assert isinstance(asset, DirectoryAsset)
    prior = pl.read_parquet(asset.path / PRIOR_FILENAME)
    assert (prior[PRIOR_WEIGHT_COLUMN].to_numpy() == 1.0).all()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/test_susie_r_finemap_task.py -k "prior_coverage or prior_parquet" -v`
Expected: FAIL — `PRIOR_FILENAME` import error, and `prior.parquet` absent.

- [ ] **Step 3: Add the constants and the coverage guard**

In `susie_r_finemap_task.py`, near the other filename constants (after `FILTERED_LD_FILENAME`):

```python
PRIOR_FILENAME = "prior.parquet"
PRIOR_WEIGHT_COLUMN = "prior_weight"
```

In `align_data`, replace the prior inner-join block. The current code is:

```python
    if prior is not None:
        prior = prior.with_row_index(name="prior_index") if prior is not None else None
        joined = joined.with_columns(
            unordered_allele_key(
                GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
            ).alias("allele_key")
        ).join(
            prior.with_columns(
                unordered_allele_key(
                    GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
                ).alias("allele_key")
            ),
            on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, "allele_key"],
        )
        prior_out = joined[_PRIOR_COL].to_numpy()
    else:
        prior_out = np.ones(len(joined))
```

Replace it with a left join plus a coverage assertion (fail fast on any gap):

```python
    if prior is not None:
        joined = joined.with_columns(
            unordered_allele_key(
                GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
            ).alias("allele_key")
        )
        n_before = len(joined)
        joined = joined.join(
            prior.with_columns(
                unordered_allele_key(
                    GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
                ).alias("allele_key")
            ),
            on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, "allele_key"],
            how="left",
        )
        missing = joined.filter(pl.col(_PRIOR_COL).is_null())
        if missing.height > 0:
            examples = missing.select(
                GWASLAB_CHROM_COL,
                GWASLAB_POS_COL,
                GWASLAB_EFFECT_ALLELE_COL,
                GWASLAB_NON_EFFECT_ALLELE_COL,
            ).head(5)
            raise ValueError(
                f"polyfun prior does not cover {missing.height} of {n_before} "
                f"(gwas intersect ld) variants; first missing:\n{examples}"
            )
        prior_out = joined[_PRIOR_COL].to_numpy()
    else:
        prior_out = np.ones(len(joined))
```

(This keeps the effective inner-join result when coverage is complete, but raises instead of silently dropping.)

- [ ] **Step 4: Write prior.parquet in execute**

Add a helper near `_save_adjustment`:

```python
def _save_prior(scratch_dir: Path, gwas_table: pl.DataFrame, prior: np.ndarray) -> None:
    """Write the per-variant prior actually passed to susie_rss, aligned to the
    filtered gwas variants. For a uniform run the values are all ones."""
    gwas_table.select(
        GWASLAB_CHROM_COL,
        GWASLAB_POS_COL,
        GWASLAB_EFFECT_ALLELE_COL,
        GWASLAB_NON_EFFECT_ALLELE_COL,
    ).with_columns(
        pl.Series(name=PRIOR_WEIGHT_COLUMN, values=prior)
    ).write_parquet(scratch_dir / PRIOR_FILENAME)
```

In `execute`, after `filter_variants_based_on_diagnostics(...)` returns the filtered `gwas_table, ld_matrix, prior` (the call already present, around the `ld_matrix = (1 - adjustment) * ...` line), and before `write_result(...)`, add:

```python
        _save_prior(scratch_dir, gwas_table=gwas_table, prior=prior)
```

`prior` here is the numpy array aligned to the filtered `gwas_table` (same length), so the `pl.Series` assignment lines up row-for-row.

- [ ] **Step 5: Run tests to verify they pass**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/test_susie_r_finemap_task.py -v`
Expected: PASS (existing tests still green; new ones pass).

- [ ] **Step 6: invoke green + commit**

```bash
pixi r invoke green > /tmp/green.log 2>&1; tail -5 /tmp/green.log
git add mecfs_bio/build_system/task/r_tasks/susie_r_finemap_task.py test_mecfs_bio/unit/build_system/task/test_susie_r_finemap_task.py
git commit -m "feat: fail-fast prior coverage guard + prior.parquet on SusieRFinemapTask

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: PolyfunExplainContrastTask

**Files:**
- Modify: `mecfs_bio/constants/polyfun_annotation_families.py`
- Create: `mecfs_bio/build_system/task/polyfun_explain/__init__.py`
- Create: `mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_contrast_task.py`
- Create: `test_mecfs_bio/unit/build_system/task/polyfun_explain/__init__.py`
- Test: `test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_contrast_task.py`

**Interfaces:**
- Consumes: SUSIE run dirs (`pip.parquet`, `filtered_gwas.parquet`, `combined_cs.parquet`, `prior.parquet` from Task 1); ridge `weights.parquet`; annotation parquet (`CHR, BP, SNP, CM, <187>`); `family_for_annotation`, `FAMILY_SHORT_LABELS`.
- Produces:
  - `PolyfunExplainContrastTask` (frozen `Task`) with `create(asset_id, susie_uniform_task, susie_polyfun_task, ridge_weights_task, annotation_parquet_task, n_important_families=3)`.
  - Output `DirectoryAsset` files (module constants): `DISPLAY_TABLE_FILENAME="display_table.parquet"`, `PER_ANNOTATION_CONTRAST_FILENAME="per_annotation_contrast.parquet"`, `PER_FAMILY_CONTRAST_FILENAME="per_family_contrast.parquet"`, `PRIOR_LIFT_FILENAME="prior_lift.parquet"`, `SELECTION_JSON_FILENAME="selection.json"`.
  - Display-table column constants used by the plot task: `DISP_CHR="chr"`, `DISP_POS="pos"`, `DISP_EA="ea"`, `DISP_NEA="nea"`, `DISP_CS_PF="cs_pf"`, `DISP_CS_U="cs_u"`, `DISP_PIP_PF="pip_pf"`, `DISP_PIP_U="pip_u"`, `DISP_LIFT="lift"`.
  - `selection.json` shape: `{"focal_variant": {"chr":.., "pos":.., "ea":.., "nea":..}, "important_families": ["conserved", ...]}` (full family names, not short labels).

- [ ] **Step 1: Add short family labels (with test)**

In `mecfs_bio/constants/polyfun_annotation_families.py`, after the `AnnotationFamily` alias, add:

```python
FAMILY_SHORT_LABELS: dict[AnnotationFamily, str] = {
    "non_synonymous": "nonsyn",
    "coding": "cod",
    "conserved": "cons",
    "promoter_or_enhancer": "prom_enh",
    "histone_marks": "hist",
    "repressed": "repr",
    "open_chromatin": "openchr",
    "maf_bins": "maf",
    "ld_related_continuous": "ld",
    "molecular_qtl": "qtl",
    "other": "other",
}
```

Add to `test_mecfs_bio/unit/constants/test_polyfun_annotation_families.py` (create if absent; if a test module already exists for this constants file, add there instead):

```python
from typing import get_args

from mecfs_bio.constants.polyfun_annotation_families import (
    AnnotationFamily,
    FAMILY_SHORT_LABELS,
)


def test_short_labels_cover_every_family():
    assert set(FAMILY_SHORT_LABELS) == set(get_args(AnnotationFamily))
```

- [ ] **Step 2: Write the failing contrast test**

Create `test_mecfs_bio/unit/build_system/task/polyfun_explain/__init__.py` (empty) and `test_polyfun_explain_contrast_task.py`:

```python
import json
from pathlib import Path, PurePath

import numpy as np
import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    DISP_CS_PF,
    DISP_CS_U,
    DISP_LIFT,
    DISP_PIP_PF,
    DISPLAY_TABLE_FILENAME,
    PER_FAMILY_CONTRAST_FILENAME,
    SELECTION_JSON_FILENAME,
    PolyfunExplainContrastTask,
)
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    COMBINED_CS_FILENAME,
    CS_COLUMN,
    FILTERED_GWAS_FILENAME,
    PIP_COLUMN,
    PIP_FILENAME,
    PRIOR_FILENAME,
    PRIOR_WEIGHT_COLUMN,
)
from mecfs_bio.build_system.wf.base_wf import make_wf

# Two real baseline-LF annotations from different families so family aggregation
# is exercised: Coding_UCSC_common -> coding, GERP.NS -> conserved.
_ANNOT_A = "Coding_UCSC_common"
_ANNOT_B = "GERP.NS"


def _write_run_dir(
    directory: Path,
    variants: pl.DataFrame,  # CHR, POS, EA, NEA
    pip: np.ndarray,
    cs_members: dict[str, list[int]],  # cs label -> row indices
    prior_weights: np.ndarray | None,
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    # pip.parquet: one PIP column, row order == variants
    pl.DataFrame({PIP_COLUMN: pip}).write_parquet(directory / PIP_FILENAME)
    variants.write_parquet(directory / FILTERED_GWAS_FILENAME)
    if prior_weights is not None:
        variants.with_columns(
            pl.Series(name=PRIOR_WEIGHT_COLUMN, values=prior_weights)
        ).write_parquet(directory / PRIOR_FILENAME)
    else:
        variants.with_columns(
            pl.Series(name=PRIOR_WEIGHT_COLUMN, values=np.ones(variants.height))
        ).write_parquet(directory / PRIOR_FILENAME)
    cs_rows = []
    for label, idxs in cs_members.items():
        for i in idxs:
            row = variants.row(i, named=True)
            cs_rows.append({**row, CS_COLUMN: label, PIP_COLUMN: float(pip[i])})
    pl.DataFrame(cs_rows).write_parquet(directory / COMBINED_CS_FILENAME)


def test_contrast_closed_form(tmp_path: Path):
    n = 6
    variants = pl.DataFrame(
        {
            "CHR": [1] * n,
            "POS": [10, 20, 30, 40, 50, 60],
            "EA": ["A"] * n,
            "NEA": ["C"] * n,
        }
    )
    # Annotation values (by CHR/BP). BP == POS.
    a = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # coding: focal only
    b = np.array([2.0, 1.0, 1.0, 1.0, 1.0, 1.0])  # conserved: focal higher
    annot = pl.DataFrame(
        {
            "CHR": [1] * n,
            "BP": [10, 20, 30, 40, 50, 60],
            "SNP": [f"rs{i}" for i in range(n)],
            "CM": [0.0] * n,
            _ANNOT_A: a,
            _ANNOT_B: b,
        }
    )
    annot_path = tmp_path / "annot.parquet"
    annot.write_parquet(annot_path)

    weights = pl.DataFrame(
        {
            "annotation": [_ANNOT_A, _ANNOT_B],
            "gamma_raw": [3.0, 0.5],
            "gamma_standardized": [0.0, 0.0],
            "family": ["coding", "conserved"],
        }
    )
    weights_path = tmp_path / "weights.parquet"
    weights.write_parquet(weights_path)

    # Uniform run: diffuse PIP; polyfun run: concentrated on variant 0 (focal).
    pip_u = np.array([0.2, 0.2, 0.2, 0.2, 0.1, 0.1])
    pip_pf = np.array([0.8, 0.05, 0.05, 0.05, 0.03, 0.02])
    uni_dir = tmp_path / "uniform"
    pf_dir = tmp_path / "polyfun"
    _write_run_dir(uni_dir, variants, pip_u, {"L1": [0, 1, 2, 3]}, prior_weights=None)
    _write_run_dir(
        pf_dir,
        variants,
        pip_pf,
        {"L1": [0]},
        prior_weights=np.array([8.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
    )

    uni_task = FakeTask(
        ResultDirectoryMeta(id=AssetId("uni"), trait="t", project="p")
    )
    pf_task = FakeTask(ResultDirectoryMeta(id=AssetId("pf"), trait="t", project="p"))
    weights_task = FakeTask(
        ReferenceFileMeta(
            group="polyfun",
            sub_group="annotations",
            sub_folder=PurePath("ridge"),
            id=AssetId("weights"),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
    )
    annot_task = FakeTask(
        SimpleFileMeta("annot", read_spec=DataFrameReadSpec(DataFrameParquetFormat()))
    )

    task = PolyfunExplainContrastTask.create(
        asset_id="contrast",
        susie_uniform_task=uni_task,
        susie_polyfun_task=pf_task,
        ridge_weights_task=weights_task,
        annotation_parquet_task=annot_task,
        n_important_families=2,
    )

    def fetch(asset_id: AssetId) -> Asset:
        mapping = {
            "uni": DirectoryAsset(uni_dir),
            "pf": DirectoryAsset(pf_dir),
            "weights": FileAsset(weights_path),
            "annot": FileAsset(annot_path),
        }
        return mapping[str(asset_id)]

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)

    # Focal variant is the max-PIP-polyfun variant (POS 10).
    selection = json.loads((result.path / SELECTION_JSON_FILENAME).read_text())
    assert selection["focal_variant"]["pos"] == 10

    # abar_c (uniform PIP-weighted over ALL variants):
    #   sum(pip_u)=1.0; abar_A = 0.2*1/1.0 = 0.2; abar_B = (0.2*2 + 0.8*1)/1.0 = 1.2
    # Focal contrast: A: 3.0*(1-0.2)=2.4 ; B: 0.5*(2-1.2)=0.4. Coding wins.
    assert selection["important_families"][0] == "coding"

    per_family = pl.read_parquet(result.path / PER_FAMILY_CONTRAST_FILENAME)
    focal_coding = per_family.filter(
        (pl.col("POS") == 10) & (pl.col("family") == "coding")
    )["family_contrast"][0]
    assert abs(focal_coding - 2.4) < 1e-9

    # Display table: union of both CS (rows 0,1,2,3), sorted desc by pip_pf.
    display = pl.read_parquet(result.path / DISPLAY_TABLE_FILENAME)
    assert display.height == 4
    assert display[DISP_PIP_PF].to_list() == sorted(
        display[DISP_PIP_PF].to_list(), reverse=True
    )
    assert display.schema["chr"] == pl.Int32
    assert display.schema["pos"] == pl.Int32
    # Focal lift = m * pi = 6 * (8/13) ~= 3.692
    focal_lift = display.filter(pl.col("pos") == 10)[DISP_LIFT][0]
    assert abs(focal_lift - 6.0 * (8.0 / 13.0)) < 1e-6
    # cs columns present, focal in both runs' CS
    focal_row = display.filter(pl.col("pos") == 10)
    assert focal_row[DISP_CS_PF][0] == 1
    assert focal_row[DISP_CS_U][0] == 1
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_contrast_task.py -v`
Expected: FAIL — module `polyfun_explain_contrast_task` does not exist.

- [ ] **Step 4: Implement the task**

Create `mecfs_bio/build_system/task/polyfun_explain/__init__.py` (empty).

Create `mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_contrast_task.py`:

```python
"""Explain a polyfun-vs-uniform SUSIE result in annotation terms.

Computes the prior lift m*pi_i and the local annotation contrast
C_c(i) = gamma_raw_c * (a_ic - abar_c), where abar_c is the uniform-run
PIP-weighted mean of annotation c over all locus variants. Aggregates the
contrast to families, selects the top families at the focal (max-PIP-polyfun)
variant, and writes the docs-facing display table plus detail tables.
"""

import json
from pathlib import Path, PurePath

import numpy as np
import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.annotation_weights.build_baseline_lf_annotation_parquet_task import (
    ANNOT_KEY_COLUMNS,
)
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    ANNOTATION_COL,
    FAMILY_COL,
    GAMMA_RAW_COL,
    WEIGHTS_PARQUET_FILENAME,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    COMBINED_CS_FILENAME,
    CS_COLUMN,
    FILTERED_GWAS_FILENAME,
    PIP_COLUMN,
    PIP_FILENAME,
    PRIOR_FILENAME,
    PRIOR_WEIGHT_COLUMN,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)
from mecfs_bio.constants.polyfun_annotation_families import FAMILY_SHORT_LABELS

DISPLAY_TABLE_FILENAME = "display_table.parquet"
PER_ANNOTATION_CONTRAST_FILENAME = "per_annotation_contrast.parquet"
PER_FAMILY_CONTRAST_FILENAME = "per_family_contrast.parquet"
PRIOR_LIFT_FILENAME = "prior_lift.parquet"
SELECTION_JSON_FILENAME = "selection.json"

DISP_CHR = "chr"
DISP_POS = "pos"
DISP_EA = "ea"
DISP_NEA = "nea"
DISP_CS_PF = "cs_pf"
DISP_CS_U = "cs_u"
DISP_PIP_PF = "pip_pf"
DISP_PIP_U = "pip_u"
DISP_LIFT = "lift"

FAMILY_CONTRAST_COL = "family_contrast"
FAMILY_SCALED_COL = "family_scaled"  # sum_c gamma_raw_c * a_ic (NOT the contrast)
CONTRAST_COL = "contrast"
_ANNOT_BP_COL = "BP"

_KEY = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
]


@frozen
class PolyfunExplainContrastTask(Task):
    meta: Meta
    susie_uniform_task: Task
    susie_polyfun_task: Task
    ridge_weights_task: Task
    annotation_parquet_task: Task
    n_important_families: int = 3

    @property
    def deps(self) -> list["Task"]:
        return [
            self.susie_uniform_task,
            self.susie_polyfun_task,
            self.ridge_weights_task,
            self.annotation_parquet_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        uni_dir = _dir(fetch, self.susie_uniform_task)
        pf_dir = _dir(fetch, self.susie_polyfun_task)

        uni_variants = _load_run_variants(uni_dir)
        pf_variants = _load_run_variants(pf_dir)
        pf_prior = pl.read_parquet(pf_dir / PRIOR_FILENAME)

        weights = _load_weights(fetch, self.ridge_weights_task)
        annot_cols = weights[ANNOTATION_COL].to_list()
        gamma = dict(zip(weights[ANNOTATION_COL], weights[GAMMA_RAW_COL]))
        family = dict(zip(weights[ANNOTATION_COL], weights[FAMILY_COL]))

        chrom = int(pf_variants[GWASLAB_CHROM_COL][0])
        bp_min = int(pf_variants[GWASLAB_POS_COL].min())
        bp_max = int(pf_variants[GWASLAB_POS_COL].max())
        annot = _load_annotations(
            fetch, self.annotation_parquet_task, chrom, bp_min, bp_max, annot_cols
        )

        # abar_c: uniform PIP-weighted mean of each annotation over all uniform vars.
        uni_annot = uni_variants.join(annot, on=_KEY, how="inner")
        w = uni_annot[PIP_COLUMN].to_numpy()
        abar = {
            c: float(np.average(uni_annot[c].to_numpy(), weights=w)) for c in annot_cols
        }

        # prior lift on the polyfun run.
        prior_w = pf_prior[PRIOR_WEIGHT_COLUMN].to_numpy()
        m = pf_variants.height
        lift = m * prior_w / prior_w.sum()
        pf_variants = pf_variants.with_columns(pl.Series(name=DISP_LIFT, values=lift))

        # attribution row set: union of the two runs' credible-set variants.
        cs_pf = _load_cs_numbers(pf_dir)
        cs_u = _load_cs_numbers(uni_dir)
        union_keys = pl.concat(
            [cs_pf.select(_KEY), cs_u.select(_KEY)], how="vertical"
        ).unique()

        pf_annot = pf_variants.join(annot, on=_KEY, how="inner")

        per_annot, per_family = _contrasts(
            pf_annot, union_keys, annot_cols, gamma, family, abar
        )
        family_scaled = _family_scaled(pf_annot, annot_cols, gamma, family)

        focal = pf_variants.sort(PIP_COLUMN, descending=True).head(1)
        focal_key = {k: focal[k][0] for k in _KEY}
        focal_families = _select_families(
            per_family, focal_key, self.n_important_families
        )

        display = _display_table(
            union_keys=union_keys,
            pf_variants=pf_variants,
            uni_variants=uni_variants,
            cs_pf=cs_pf,
            cs_u=cs_u,
            family_scaled=family_scaled,
            focal_families=focal_families,
        )

        per_annot.write_parquet(scratch_dir / PER_ANNOTATION_CONTRAST_FILENAME)
        per_family.write_parquet(scratch_dir / PER_FAMILY_CONTRAST_FILENAME)
        pf_variants.select(*_KEY, DISP_LIFT).write_parquet(
            scratch_dir / PRIOR_LIFT_FILENAME
        )
        display.write_parquet(scratch_dir / DISPLAY_TABLE_FILENAME)
        (scratch_dir / SELECTION_JSON_FILENAME).write_text(
            json.dumps(
                {
                    "focal_variant": {
                        "chr": int(focal_key[GWASLAB_CHROM_COL]),
                        "pos": int(focal_key[GWASLAB_POS_COL]),
                        "ea": focal_key[GWASLAB_EFFECT_ALLELE_COL],
                        "nea": focal_key[GWASLAB_NON_EFFECT_ALLELE_COL],
                    },
                    "important_families": focal_families,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        susie_uniform_task: Task,
        susie_polyfun_task: Task,
        ridge_weights_task: Task,
        annotation_parquet_task: Task,
        n_important_families: int = 3,
    ) -> "PolyfunExplainContrastTask":
        source_meta = susie_polyfun_task.meta
        if not isinstance(source_meta, ResultDirectoryMeta):
            raise ValueError(f"Unknown meta for polyfun susie task: {source_meta}")
        meta = ResultDirectoryMeta(
            id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=PurePath("analysis"),
        )
        return cls(
            meta=meta,
            susie_uniform_task=susie_uniform_task,
            susie_polyfun_task=susie_polyfun_task,
            ridge_weights_task=ridge_weights_task,
            annotation_parquet_task=annotation_parquet_task,
            n_important_families=n_important_families,
        )


def _dir(fetch: Fetch, task: Task) -> Path:
    asset = fetch(task.asset_id)
    assert isinstance(asset, DirectoryAsset)
    return asset.path


def _load_run_variants(run_dir: Path) -> pl.DataFrame:
    """filtered_gwas keyed rows + the run's PIP, in the same order."""
    gwas = pl.read_parquet(run_dir / FILTERED_GWAS_FILENAME).select(_KEY)
    pip = pl.read_parquet(run_dir / PIP_FILENAME).select(PIP_COLUMN)
    return gwas.hstack(pip)


def _load_cs_numbers(run_dir: Path) -> pl.DataFrame:
    """One row per credible-set variant with its 1-based L-index (lowest if many)."""
    cs = pl.read_parquet(run_dir / COMBINED_CS_FILENAME)
    if cs.height == 0:
        return pl.DataFrame(schema={**{k: cs.schema.get(k, pl.Int64) for k in _KEY}})
    return (
        cs.with_columns(
            pl.col(CS_COLUMN).str.replace("L", "").cast(pl.Int32).alias("cs_number")
        )
        .group_by(_KEY)
        .agg(pl.col("cs_number").min())
    )


def _load_weights(fetch: Fetch, task: Task) -> pl.DataFrame:
    asset = fetch(task.asset_id)
    assert isinstance(asset, (FileAsset, DirectoryAsset))
    path = asset.path if isinstance(asset, FileAsset) else asset.path / WEIGHTS_PARQUET_FILENAME
    return pl.read_parquet(path)


def _load_annotations(
    fetch: Fetch,
    task: Task,
    chrom: int,
    bp_min: int,
    bp_max: int,
    annot_cols: list[str],
) -> pl.DataFrame:
    """Locus-windowed annotation slice, keyed like the run variants (BP -> POS)."""
    asset = fetch(task.asset_id)
    frame = (
        scan_dataframe_asset(asset, task.meta)
        .filter(
            (pl.col("CHR") == chrom)
            & (pl.col(_ANNOT_BP_COL) >= bp_min)
            & (pl.col(_ANNOT_BP_COL) <= bp_max)
        )
        .select("CHR", _ANNOT_BP_COL, *annot_cols)
        .collect()
        .to_polars()
    )
    # annotation has no alleles; broadcast to the run's alleles via the join key by
    # renaming BP->POS and letting the (CHR,POS) match carry EA/NEA from the run.
    return frame.rename({_ANNOT_BP_COL: GWASLAB_POS_COL})


def _contrasts(
    pf_annot: pl.DataFrame,
    union_keys: pl.DataFrame,
    annot_cols: list[str],
    gamma: dict[str, float],
    family: dict[str, str],
    abar: dict[str, float],
) -> tuple[pl.DataFrame, pl.DataFrame]:
    rows = pf_annot.join(union_keys, on=_KEY, how="inner")
    long = rows.unpivot(
        on=annot_cols, index=_KEY, variable_name=ANNOTATION_COL, value_name="a_ic"
    ).with_columns(
        (
            pl.col(ANNOTATION_COL).replace_strict(gamma)
            * (pl.col("a_ic") - pl.col(ANNOTATION_COL).replace_strict(abar))
        ).alias(CONTRAST_COL),
        pl.col(ANNOTATION_COL).replace_strict(family).alias(FAMILY_COL),
    )
    per_annotation = long.select(*_KEY, ANNOTATION_COL, FAMILY_COL, CONTRAST_COL)
    per_family = (
        long.group_by([*_KEY, FAMILY_COL])
        .agg(pl.col(CONTRAST_COL).sum().alias(FAMILY_CONTRAST_COL))
        .sort([*_KEY, FAMILY_COL])
    )
    return per_annotation, per_family


def _family_scaled(
    pf_annot: pl.DataFrame,
    annot_cols: list[str],
    gamma: dict[str, float],
    family: dict[str, str],
) -> pl.DataFrame:
    """Per variant per family: sum_c gamma_raw_c * a_ic (raw scaled value)."""
    long = pf_annot.unpivot(
        on=annot_cols, index=_KEY, variable_name=ANNOTATION_COL, value_name="a_ic"
    ).with_columns(
        (pl.col(ANNOTATION_COL).replace_strict(gamma) * pl.col("a_ic")).alias(
            "scaled"
        ),
        pl.col(ANNOTATION_COL).replace_strict(family).alias(FAMILY_COL),
    )
    return (
        long.group_by([*_KEY, FAMILY_COL])
        .agg(pl.col("scaled").sum().alias(FAMILY_SCALED_COL))
        .sort([*_KEY, FAMILY_COL])
    )


def _select_families(
    per_family: pl.DataFrame, focal_key: dict, n: int
) -> list[str]:
    focal = per_family
    for k, v in focal_key.items():
        focal = focal.filter(pl.col(k) == v)
    return (
        focal.sort(FAMILY_CONTRAST_COL, descending=True)
        .head(n)[FAMILY_COL]
        .to_list()
    )


def _display_table(
    union_keys: pl.DataFrame,
    pf_variants: pl.DataFrame,
    uni_variants: pl.DataFrame,
    cs_pf: pl.DataFrame,
    cs_u: pl.DataFrame,
    family_scaled: pl.DataFrame,
    focal_families: list[str],
) -> pl.DataFrame:
    out = (
        union_keys.join(
            pf_variants.select(*_KEY, pl.col(PIP_COLUMN).alias(DISP_PIP_PF), DISP_LIFT),
            on=_KEY,
            how="left",
        )
        .join(
            uni_variants.select(*_KEY, pl.col(PIP_COLUMN).alias(DISP_PIP_U)),
            on=_KEY,
            how="left",
        )
        .join(
            cs_pf.select(*_KEY, pl.col("cs_number").alias(DISP_CS_PF)),
            on=_KEY,
            how="left",
        )
        .join(
            cs_u.select(*_KEY, pl.col("cs_number").alias(DISP_CS_U)),
            on=_KEY,
            how="left",
        )
    )
    for fam in focal_families:
        col = FAMILY_SHORT_LABELS[fam]
        fam_col = (
            family_scaled.filter(pl.col(FAMILY_COL) == fam)
            .select(*_KEY, pl.col(FAMILY_SCALED_COL).alias(col))
        )
        out = out.join(fam_col, on=_KEY, how="left")
    out = out.rename(
        {
            GWASLAB_CHROM_COL: DISP_CHR,
            GWASLAB_POS_COL: DISP_POS,
            GWASLAB_EFFECT_ALLELE_COL: DISP_EA,
            GWASLAB_NON_EFFECT_ALLELE_COL: DISP_NEA,
        }
    ).with_columns(
        pl.col(DISP_CHR).cast(pl.Int32), pl.col(DISP_POS).cast(pl.Int32)
    )
    ordered = [
        DISP_CHR, DISP_POS, DISP_EA, DISP_NEA, DISP_CS_PF, DISP_CS_U,
        DISP_PIP_PF, DISP_PIP_U, DISP_LIFT,
    ] + [FAMILY_SHORT_LABELS[f] for f in focal_families]
    return out.select(ordered).sort(DISP_PIP_PF, descending=True, nulls_last=True)
```

Note on `replace_strict`: it maps every value via the dict and raises on an unmapped key — the fail-fast behavior we want if a weights/annotation mismatch ever occurs. `unpivot` is the current polars name for melt.

- [ ] **Step 5: Run test to verify it passes**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_contrast_task.py test_mecfs_bio/unit/constants/test_polyfun_annotation_families.py -v`
Expected: PASS. If a closed-form value is off, recompute abar_c/contrast by hand against the fixture before touching the code.

- [ ] **Step 6: invoke green + commit**

```bash
pixi r invoke green > /tmp/green.log 2>&1; tail -5 /tmp/green.log
git add mecfs_bio/constants/polyfun_annotation_families.py mecfs_bio/build_system/task/polyfun_explain/ test_mecfs_bio/unit/build_system/task/polyfun_explain/ test_mecfs_bio/unit/constants/test_polyfun_annotation_families.py
git commit -m "feat: PolyfunExplainContrastTask (prior lift + annotation/family contrast + display table)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: PolyfunExplainPlotTask

**Files:**
- Create: `mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_plot_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_plot_task.py`

**Interfaces:**
- Consumes: the two SUSIE run dirs, the contrast task dir (`selection.json`, and per-variant family-scaled values it recomputes from the annotation parquet + weights for the plot), the annotation parquet (for CM and family-scaled tracks), and a gene-info task.
- Produces: `PolyfunExplainPlotTask` with `create(asset_id, susie_uniform_task, susie_polyfun_task, contrast_task, annotation_parquet_task, gene_info_task, gene_info_pipe=IdentityPipe(), n_family_panels=3)`; output `DirectoryAsset` with module constants `PLOT_PNG_FILENAME="explain_plot.png"` and `PLOT_SVG_FILENAME="explain_plot.svg"`.

The plot reuses the family-scaled computation from Task 2 by importing `_family_scaled`, `_load_annotations`, `_load_run_variants`, `FAMILY_SCALED_COL`, `FAMILY_COL`, and reads `selection.json` for the focal variant + families so the figure matches the tables.

- [ ] **Step 1: Write the failing smoke test**

Create `test_polyfun_explain_plot_task.py`:

```python
from pathlib import Path

import numpy as np
import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_plot_task import (
    PLOT_PNG_FILENAME,
    PLOT_SVG_FILENAME,
    PolyfunExplainPlotTask,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
    GENE_INFO_STRAND_COL,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from test_mecfs_bio.unit.build_system.task.polyfun_explain.test_polyfun_explain_contrast_task import (
    build_synthetic_explain_inputs,
)


def test_plot_writes_png_and_svg(tmp_path: Path):
    inputs = build_synthetic_explain_inputs(tmp_path)  # helper extracted in Task 3 step 3

    gene_info = pl.DataFrame(
        {
            GENE_INFO_CHROM_COL: [1],
            GENE_INFO_START_COL: [5],
            GENE_INFO_END_COL: [65],
            GENE_INFO_STRAND_COL: ["+"],
            GENE_INFO_NAME_COL: ["GENE1"],
        }
    )
    gene_path = tmp_path / "genes.parquet"
    gene_info.write_parquet(gene_path)
    gene_task = FakeTask(
        SimpleFileMeta("genes", read_spec=DataFrameReadSpec(DataFrameParquetFormat()))
    )

    plot_task = PolyfunExplainPlotTask.create(
        asset_id="plot",
        susie_uniform_task=inputs.uni_task,
        susie_polyfun_task=inputs.pf_task,
        contrast_task=inputs.contrast_task,
        annotation_parquet_task=inputs.annot_task,
        gene_info_task=gene_task,
        gene_info_pipe=IdentityPipe(),
        n_family_panels=2,
    )

    def fetch(asset_id: AssetId) -> Asset:
        mapping = dict(inputs.fetch_map)
        mapping["genes"] = FileAsset(gene_path)
        return mapping[str(asset_id)]

    scratch = tmp_path / "plot_scratch"
    scratch.mkdir()
    result = plot_task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)
    assert (result.path / PLOT_PNG_FILENAME).is_file()
    assert (result.path / PLOT_SVG_FILENAME).is_file()
```

- [ ] **Step 2: Extract a shared fixture builder in the contrast test**

Refactor the Task 2 test so its synthetic setup is reusable. Add to `test_polyfun_explain_contrast_task.py` a helper returning the tasks + a `fetch` map + the already-run contrast task directory (executed once), as a small frozen holder:

```python
from attrs import frozen


@frozen
class _ExplainInputs:
    uni_task: object
    pf_task: object
    weights_task: object
    annot_task: object
    contrast_task: object
    fetch_map: tuple  # tuple of (str asset_id, Asset)


def build_synthetic_explain_inputs(tmp_path: Path) -> _ExplainInputs:
    # (Move the body of test_contrast_closed_form's setup here, execute the
    # contrast task, and return the holder. test_contrast_closed_form then calls
    # this and keeps only its assertions.)
    ...
```

Have `test_contrast_closed_form` call `build_synthetic_explain_inputs` and assert against the returned/executed contrast directory. This keeps one synthetic dataset shared across Tasks 2 and 3 (repo convention: share creation across tests).

- [ ] **Step 3: Run the plot test to verify it fails**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_plot_task.py -v`
Expected: FAIL — plot module does not exist.

- [ ] **Step 4: Implement the plot task**

Create `mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_plot_task.py`:

```python
"""Stacked explainability figure for a polyfun-vs-uniform SUSIE result.

Panels, top to bottom, sharing the genomic-position x-axis:
  1. Manhattan (-log10 p), points colored by LD r^2 with the min-p lead variant,
     with local recombination rate (dCM/dBP from the annotation CM column) on a
     secondary axis.
  2..(1+x). One panel per important family: sum_c gamma_raw_c * a_ic across the
     locus (raw scaled value; the contrast is the profile).
  3. Prior fold m*pi_i (log scale).
  4. PIP, uniform run.
  5. PIP, polyfun run.
  6. Genes.

Writes both explain_plot.png and explain_plot.svg. Inspired by SusieStackPlotTask
but independent of it.
"""

import json
from pathlib import Path, PurePath

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
from attrs import frozen  # noqa: E402
from matplotlib import gridspec  # noqa: E402

from mecfs_bio.build_system.asset.base_asset import Asset  # noqa: E402
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset  # noqa: E402
from mecfs_bio.build_system.asset.file_asset import FileAsset  # noqa: E402
from mecfs_bio.build_system.meta.asset_id import AssetId  # noqa: E402
from mecfs_bio.build_system.meta.meta import Meta  # noqa: E402
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (  # noqa: E402
    scan_dataframe_asset,
)
from mecfs_bio.build_system.meta.result_directory_meta import (  # noqa: E402
    ResultDirectoryMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch  # noqa: E402
from mecfs_bio.build_system.task.base_task import Task  # noqa: E402
from mecfs_bio.build_system.task.pipes.data_processing_pipe import (  # noqa: E402
    DataProcessingPipe,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe  # noqa: E402
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (  # noqa: E402
    DISP_LIFT,
    FAMILY_COL,
    FAMILY_SCALED_COL,
    SELECTION_JSON_FILENAME,
    _family_scaled,
    _load_annotations,
    _load_run_variants,
)
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (  # noqa: E402
    FILTERED_GWAS_FILENAME,
    FILTERED_LD_FILENAME,
    PIP_COLUMN,
    PRIOR_FILENAME,
    PRIOR_WEIGHT_COLUMN,
)
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (  # noqa: E402
    ANNOTATION_COL,
    FAMILY_COL as WEIGHTS_FAMILY_COL,
    GAMMA_RAW_COL,
    WEIGHTS_PARQUET_FILENAME,
)
from mecfs_bio.build_system.wf.base_wf import WF  # noqa: E402
from mecfs_bio.constants.gwaslab_constants import (  # noqa: E402
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (  # noqa: E402
    GENE_INFO_END_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
)

PLOT_PNG_FILENAME = "explain_plot.png"
PLOT_SVG_FILENAME = "explain_plot.svg"
_ANNOT_BP_COL = "BP"


@frozen
class PolyfunExplainPlotTask(Task):
    meta: Meta
    susie_uniform_task: Task
    susie_polyfun_task: Task
    contrast_task: Task
    annotation_parquet_task: Task
    gene_info_task: Task
    gene_info_pipe: DataProcessingPipe = IdentityPipe()
    n_family_panels: int = 3

    @property
    def deps(self) -> list["Task"]:
        return [
            self.susie_uniform_task,
            self.susie_polyfun_task,
            self.contrast_task,
            self.annotation_parquet_task,
            self.gene_info_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        uni_dir = _dir(fetch, self.susie_uniform_task)
        pf_dir = _dir(fetch, self.susie_polyfun_task)
        contrast_dir = _dir(fetch, self.contrast_task)
        selection = json.loads(
            (contrast_dir / SELECTION_JSON_FILENAME).read_text()
        )
        families = selection["important_families"][: self.n_family_panels]

        pf_variants = _load_run_variants(pf_dir).sort(GWASLAB_POS_COL)
        uni_variants = _load_run_variants(uni_dir)
        pf_full = pl.read_parquet(pf_dir / FILTERED_GWAS_FILENAME)
        ld = np.load(pf_dir / FILTERED_LD_FILENAME)
        prior_w = pl.read_parquet(pf_dir / PRIOR_FILENAME)[PRIOR_WEIGHT_COLUMN].to_numpy()

        weights = _load_weights(fetch, self.contrast_task, pf_dir, fetch2=fetch)
        # weights actually come from the ridge asset via the contrast task's dep
        # graph; see _load_weights_for_plot below.

        chrom = int(pf_variants[GWASLAB_CHROM_COL][0])
        bp_min = int(pf_variants[GWASLAB_POS_COL].min())
        bp_max = int(pf_variants[GWASLAB_POS_COL].max())

        # (implementation continues: build panels; see steps below)
        _render(
            scratch_dir=scratch_dir,
            pf_variants=pf_variants,
            uni_variants=uni_variants,
            pf_full=pf_full,
            ld=ld,
            prior_w=prior_w,
            families=families,
            fetch=fetch,
            annotation_parquet_task=self.annotation_parquet_task,
            gene_info_task=self.gene_info_task,
            gene_info_pipe=self.gene_info_pipe,
            chrom=chrom,
            bp_min=bp_min,
            bp_max=bp_max,
        )
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        susie_uniform_task: Task,
        susie_polyfun_task: Task,
        contrast_task: Task,
        annotation_parquet_task: Task,
        gene_info_task: Task,
        gene_info_pipe: DataProcessingPipe = IdentityPipe(),
        n_family_panels: int = 3,
    ) -> "PolyfunExplainPlotTask":
        source_meta = susie_polyfun_task.meta
        if not isinstance(source_meta, ResultDirectoryMeta):
            raise ValueError(f"Unknown meta for polyfun susie task: {source_meta}")
        meta = ResultDirectoryMeta(
            id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=PurePath("analysis"),
        )
        return cls(
            meta=meta,
            susie_uniform_task=susie_uniform_task,
            susie_polyfun_task=susie_polyfun_task,
            contrast_task=contrast_task,
            annotation_parquet_task=annotation_parquet_task,
            gene_info_task=gene_info_task,
            gene_info_pipe=gene_info_pipe,
            n_family_panels=n_family_panels,
        )
```

The renderer and helpers (same module):

```python
def _dir(fetch: Fetch, task: Task) -> Path:
    asset = fetch(task.asset_id)
    assert isinstance(asset, DirectoryAsset)
    return asset.path


def _weights_from_contrast_dep(fetch: Fetch, ridge_weights_task: Task) -> pl.DataFrame:
    asset = fetch(ridge_weights_task.asset_id)
    path = (
        asset.path
        if isinstance(asset, FileAsset)
        else asset.path / WEIGHTS_PARQUET_FILENAME
    )
    return pl.read_parquet(path)


def _render(
    scratch_dir: Path,
    pf_variants: pl.DataFrame,
    uni_variants: pl.DataFrame,
    pf_full: pl.DataFrame,
    ld: np.ndarray,
    prior_w: np.ndarray,
    families: list[str],
    fetch: Fetch,
    annotation_parquet_task: Task,
    gene_info_task: Task,
    gene_info_pipe: DataProcessingPipe,
    chrom: int,
    bp_min: int,
    bp_max: int,
) -> None:
    n_panels = 1 + len(families) + 3 + 1  # manhattan + families + lift + 2 pip + genes
    fig = plt.figure(figsize=(10, 2.0 * n_panels))
    gs = gridspec.GridSpec(nrows=n_panels, ncols=1, hspace=0.4)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(n_panels)]
    x = pf_full[GWASLAB_POS_COL].to_numpy()

    # Panel 1: Manhattan colored by LD with lead (min-p) variant + recomb rate.
    z = (pf_full[GWASLAB_BETA_COL] / pf_full[GWASLAB_SE_COL]).to_numpy()
    neglogp = -np.log10(2.0 * _norm_sf(np.abs(z)))
    lead = int(np.argmax(np.abs(z)))
    r2 = ld[lead, :] ** 2
    ax0 = axes[0]
    sc = ax0.scatter(x, neglogp, c=r2, cmap="viridis", vmin=0, vmax=1, s=12)
    ax0.set_ylabel("-log10 p")
    fig.colorbar(sc, ax=ax0, label="r2 w/ lead")

    annot = _load_annotations(
        fetch, annotation_parquet_task, chrom, bp_min, bp_max, annot_cols=[]
    )  # annot_cols=[] -> CM only path; see note
    cm = _load_cm(fetch, annotation_parquet_task, chrom, bp_min, bp_max)
    ax0b = ax0.twinx()
    _plot_recomb(ax0b, cm)

    # Family panels
    weights = _weights_from_contrast_dep_via_annotation(fetch, annotation_parquet_task)
    # NOTE: weights come from the ridge asset; wired through the generator so the
    # plot task depends on it. See Task 4 wiring; in tests inject it directly.
    scaled = _family_scaled_for_plot(fetch, annotation_parquet_task, pf_variants, weights)
    for k, fam in enumerate(families):
        axk = axes[1 + k]
        fam_series = scaled.filter(pl.col(FAMILY_COL) == fam).sort(GWASLAB_POS_COL)
        axk.plot(
            fam_series[GWASLAB_POS_COL].to_numpy(),
            fam_series[FAMILY_SCALED_COL].to_numpy(),
        )
        axk.set_ylabel(fam)

    # Prior fold (log)
    lift = pf_variants.height * prior_w / prior_w.sum()
    ax_lift = axes[1 + len(families)]
    ax_lift.scatter(pf_full[GWASLAB_POS_COL].to_numpy(), lift, s=12)
    ax_lift.set_yscale("log")
    ax_lift.set_ylabel("prior fold")

    # PIP uniform / polyfun
    ax_u = axes[2 + len(families)]
    ax_u.scatter(
        uni_variants[GWASLAB_POS_COL].to_numpy(),
        uni_variants[PIP_COLUMN].to_numpy(),
        s=12,
    )
    ax_u.set_ylabel("PIP uniform")
    ax_pf = axes[3 + len(families)]
    ax_pf.scatter(
        pf_variants[GWASLAB_POS_COL].to_numpy(),
        pf_variants[PIP_COLUMN].to_numpy(),
        s=12,
    )
    ax_pf.set_ylabel("PIP polyfun")

    # Genes
    genes = (
        gene_info_pipe.process(
            scan_dataframe_asset(fetch(gene_info_task.asset_id), gene_info_task.meta)
        )
        .collect()
        .to_polars()
    )
    _plot_genes(axes[-1], genes, bp_min, bp_max)
    axes[-1].set_xlabel(f"chr{chrom} position (bp)")

    fig.savefig(scratch_dir / PLOT_PNG_FILENAME, dpi=150, bbox_inches="tight")
    fig.savefig(scratch_dir / PLOT_SVG_FILENAME, bbox_inches="tight")
    plt.close(fig)
```

Small helpers in the same module:

```python
def _norm_sf(z: np.ndarray) -> np.ndarray:
    from scipy.stats import norm

    return norm.sf(z)


def _load_cm(
    fetch: Fetch, task: Task, chrom: int, bp_min: int, bp_max: int
) -> pl.DataFrame:
    frame = (
        scan_dataframe_asset(fetch(task.asset_id), task.meta)
        .filter(
            (pl.col("CHR") == chrom)
            & (pl.col(_ANNOT_BP_COL) >= bp_min)
            & (pl.col(_ANNOT_BP_COL) <= bp_max)
        )
        .select("CHR", _ANNOT_BP_COL, "CM")
        .collect()
        .to_polars()
        .sort(_ANNOT_BP_COL)
    )
    return frame


def _plot_recomb(ax, cm: pl.DataFrame) -> None:
    bp = cm[_ANNOT_BP_COL].to_numpy().astype(float)
    cmv = cm["CM"].to_numpy().astype(float)
    if len(bp) >= 2:
        rate = np.gradient(cmv, bp) * 1e6  # cM/Mb
        ax.plot(bp, rate, color="tab:red", alpha=0.4, linewidth=1)
    ax.set_ylabel("cM/Mb", color="tab:red")


def _plot_genes(ax, genes: pl.DataFrame, bp_min: int, bp_max: int) -> None:
    row = 0
    for g in genes.iter_rows(named=True):
        start = max(int(g[GENE_INFO_START_COL]), bp_min)
        end = min(int(g[GENE_INFO_END_COL]), bp_max)
        ax.plot([start, end], [row, row], linewidth=4)
        ax.text((start + end) / 2, row + 0.1, g[GENE_INFO_NAME_COL], ha="center")
        row += 1
    ax.set_yticks([])
    ax.set_ylabel("genes")
```

Wiring the weights into the plot task: give the plot task an explicit `ridge_weights_task` dependency instead of routing through the contrast task. Update the `create()` signature and `deps` to include `ridge_weights_task`, and replace the `_weights_from_contrast_dep_via_annotation` / `_family_scaled_for_plot` sketch with a direct call:

```python
weights = _weights_from_contrast_dep(fetch, self.ridge_weights_task)
annot = _load_annotations(fetch, self.annotation_parquet_task, chrom, bp_min, bp_max, weights[ANNOTATION_COL].to_list())
pf_annot = pf_variants.join(annot, on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, "EA", "NEA"], how="inner")
gamma = dict(zip(weights[ANNOTATION_COL], weights[GAMMA_RAW_COL]))
family = dict(zip(weights[ANNOTATION_COL], weights[WEIGHTS_FAMILY_COL]))
scaled = _family_scaled(pf_annot, weights[ANNOTATION_COL].to_list(), gamma, family)
```

Update the test's `PolyfunExplainPlotTask.create(...)` call and `fetch` map to include `ridge_weights_task=inputs.weights_task` and map `"weights"` to `FileAsset(weights_path)`. (This keeps the plot task's data path identical to the contrast task's, so the family panels match the tables.)

Simplify: remove the `_weights_from_contrast_dep_via_annotation`, `_family_scaled_for_plot`, and the stray `_load_weights`/`weights` lines in the first `execute` sketch; the plot task's real dependencies are `[uniform, polyfun, contrast, ridge_weights, annotation, gene_info]`, and only `selection.json` is read from the contrast dir.

- [ ] **Step 5: Run test to verify it passes**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_plot_task.py -v`
Expected: PASS — both files exist.

- [ ] **Step 6: invoke green + commit**

```bash
pixi r invoke green > /tmp/green.log 2>&1; tail -5 /tmp/green.log
git add mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_plot_task.py test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_plot_task.py test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_contrast_task.py
git commit -m "feat: PolyfunExplainPlotTask (stacked explainability figure, PNG+SVG)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: Inner + Outer asset generators

**Files:**
- Create: `mecfs_bio/asset_generator/polyfun_explain_fine_mapping_asset_generator.py`
- Test: `test_mecfs_bio/unit/asset_generator/test_polyfun_explain_fine_mapping_asset_generator.py`

**Interfaces:**
- Consumes: `SusieRFinemapTask`, `PriorInfo`, `BroadInstituteFormatLDMatrix`, `PolyfunExplainContrastTask`, `PolyfunExplainPlotTask`, `COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS`, `create_prior_col_pipe`, `POLYFUN_PRIOR_COL`, `BASELINE_LF_ANNOTATION_MATRIX`, `BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS`, and the existing shared-setup helpers from `fine_mapping_asset_generator` / `ukbb_broad_ld_matrix_generator`.
- Produces:
  - `RunConfig` frozen (`label: str`, `max_credible_sets: int`, `z_score_filtering_threshold: float = 2.0`).
  - `RUN_CONFIGS: tuple[RunConfig, ...]` = L=1, L=2, L=10, L=10-strict (z=1.0).
  - `SharedFineMapInputs` frozen holding the per-locus shared tasks (harmonized sumstats, renamed ld labels, ld matrix, gene info, effective N, base name, q_factor).
  - `PolyfunExplainGroup` frozen (`susie_uniform`, `susie_polyfun`, `contrast`, `plot`).
  - `generate_polyfun_explain_group(shared, config) -> PolyfunExplainGroup`.
  - `build_explainability_groups(shared, configs=RUN_CONFIGS) -> list[PolyfunExplainGroup]`.
  - `PolyfunExplainOuterGroup` frozen with `groups: list[PolyfunExplainGroup]` and `terminal_tasks() -> list[Task]`.
  - `generate_assets_polyfun_explain_fine_map(...) -> PolyfunExplainOuterGroup` (does the LD-interval lookup + harmonization, then `build_explainability_groups`).

- [ ] **Step 1: Write the failing test (pure builder, FakeTask inputs)**

Create `test_mecfs_bio/unit/asset_generator/test_polyfun_explain_fine_mapping_asset_generator.py`:

```python
from pathlib import PurePath

from mecfs_bio.asset_generator.polyfun_explain_fine_mapping_asset_generator import (
    RUN_CONFIGS,
    SharedFineMapInputs,
    build_explainability_groups,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask


def _shared() -> SharedFineMapInputs:
    harmonized = FakeTask(
        FilteredGWASDataMeta(
            id=AssetId("harmonized"),
            trait="mecfs",
            project="decodeme",
            sub_dir=PurePath("processed_gwas_data"),
        )
    )
    ld_labels = FakeTask(SimpleFileMeta("ld_labels"))
    ld_matrix = FakeTask(SimpleFileMeta("ld_matrix"))
    gene_info = FakeTask(SimpleFileMeta("genes"))
    return SharedFineMapInputs(
        base_name="mecfs_chr1_174",
        harmonized_sumstats_task=harmonized,
        ld_labels_task=ld_labels,
        ld_matrix_task=ld_matrix,
        gene_info_task=gene_info,
        effective_sample_size=10000,
        q_factor=100,
    )


def test_builds_four_groups_eight_susie_runs():
    groups = build_explainability_groups(_shared())
    assert len(groups) == len(RUN_CONFIGS) == 4
    susie_ids = []
    for g in groups:
        susie_ids += [g.susie_uniform.asset_id, g.susie_polyfun.asset_id]
    assert len(susie_ids) == 8
    assert len(set(susie_ids)) == 8  # all distinct


def test_polyfun_run_has_prior_uniform_does_not():
    group = build_explainability_groups(_shared())[0]
    assert group.susie_polyfun.prior_info is not None
    assert group.susie_uniform.prior_info is None
```

Note: `FilteredGWASDataMeta` requires `id, trait, project, sub_dir` (read_spec/extension default); `SusieRFinemapTask.create` reads `.trait` and `.project`.

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi r python -m pytest test_mecfs_bio/unit/asset_generator/test_polyfun_explain_fine_mapping_asset_generator.py -v`
Expected: FAIL — generator module does not exist.

- [ ] **Step 3: Implement the generators**

Create `mecfs_bio/asset_generator/polyfun_explain_fine_mapping_asset_generator.py`:

```python
"""Asset generator for polyfun explainability fine-mapping.

Inner generator: given a locus's shared inputs and one run config, produce a
matched pair of SUSIE runs (polyfun prior + uniform) plus the contrast and plot
tasks explaining the pair. Outer generator: call the inner generator for each of
the four run configs (L=1, L=2, L=10, L=10-strict) -> 8 SUSIE runs per locus.

Inspired by fine_mapping_asset_generator; the existing generator is unchanged.
"""

from pathlib import PurePath

from attrs import frozen

from mecfs_bio.asset_generator.fine_mapping_asset_generator import (
    _build_shared_locus_inputs,  # see Step 3b
)
from mecfs_bio.assets.reference_data.polyfun.annotations.annotation_ridge_weights import (
    BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS,
)
from mecfs_bio.assets.reference_data.polyfun.annotations.baseline_lf_annotations import (
    BASELINE_LF_ANNOTATION_MATRIX,
)
from mecfs_bio.assets.reference_data.polyfun.precomputed_prior.polyfun_precomputed_prior import (
    COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS,
    POLYFUN_PRIOR_COL,
    create_prior_col_pipe,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    PolyfunExplainContrastTask,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_plot_task import (
    PolyfunExplainPlotTask,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    BroadInstituteFormatLDMatrix,
    PriorInfo,
    SusieRFinemapTask,
)


@frozen
class RunConfig:
    label: str
    max_credible_sets: int
    z_score_filtering_threshold: float = 2.0


RUN_CONFIGS: tuple[RunConfig, ...] = (
    RunConfig(label="l1", max_credible_sets=1),
    RunConfig(label="l2", max_credible_sets=2),
    RunConfig(label="l10", max_credible_sets=10),
    RunConfig(label="l10_strict", max_credible_sets=10, z_score_filtering_threshold=1.0),
)


@frozen
class SharedFineMapInputs:
    base_name: str
    harmonized_sumstats_task: Task
    ld_labels_task: Task
    ld_matrix_task: Task
    gene_info_task: Task
    effective_sample_size: int
    q_factor: int = 100


@frozen
class PolyfunExplainGroup:
    susie_uniform: Task
    susie_polyfun: Task
    contrast: Task
    plot: Task


@frozen
class PolyfunExplainOuterGroup:
    groups: list[PolyfunExplainGroup]

    def terminal_tasks(self) -> list[Task]:
        out: list[Task] = []
        for g in self.groups:
            out += [g.susie_uniform, g.susie_polyfun, g.contrast, g.plot]
        return out


def generate_polyfun_explain_group(
    shared: SharedFineMapInputs, config: RunConfig
) -> PolyfunExplainGroup:
    stem = f"{shared.base_name}_{config.label}"
    prior_info = PriorInfo(
        prior_task=COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS,
        prior_pipe=create_prior_col_pipe(shared.q_factor),
        prior_col=POLYFUN_PRIOR_COL,
    )
    susie_uniform = SusieRFinemapTask.create(
        asset_id=f"{stem}_susie_uniform",
        gwas_data_task=shared.harmonized_sumstats_task,
        ld_labels_task=shared.ld_labels_task,
        ld_matrix_source=BroadInstituteFormatLDMatrix(shared.ld_matrix_task),
        effective_sample_size=shared.effective_sample_size,
        max_credible_sets=config.max_credible_sets,
        z_score_filtering_threshold=config.z_score_filtering_threshold,
        prior_info=None,
    )
    susie_polyfun = SusieRFinemapTask.create(
        asset_id=f"{stem}_susie_polyfun",
        gwas_data_task=shared.harmonized_sumstats_task,
        ld_labels_task=shared.ld_labels_task,
        ld_matrix_source=BroadInstituteFormatLDMatrix(shared.ld_matrix_task),
        effective_sample_size=shared.effective_sample_size,
        max_credible_sets=config.max_credible_sets,
        z_score_filtering_threshold=config.z_score_filtering_threshold,
        prior_info=prior_info,
    )
    contrast = PolyfunExplainContrastTask.create(
        asset_id=f"{stem}_explain_contrast",
        susie_uniform_task=susie_uniform,
        susie_polyfun_task=susie_polyfun,
        ridge_weights_task=BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS,
        annotation_parquet_task=BASELINE_LF_ANNOTATION_MATRIX,
    )
    plot = PolyfunExplainPlotTask.create(
        asset_id=f"{stem}_explain_plot",
        susie_uniform_task=susie_uniform,
        susie_polyfun_task=susie_polyfun,
        contrast_task=contrast,
        ridge_weights_task=BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS,
        annotation_parquet_task=BASELINE_LF_ANNOTATION_MATRIX,
        gene_info_task=shared.gene_info_task,
        gene_info_pipe=IdentityPipe(),
    )
    return PolyfunExplainGroup(
        susie_uniform=susie_uniform,
        susie_polyfun=susie_polyfun,
        contrast=contrast,
        plot=plot,
    )


def build_explainability_groups(
    shared: SharedFineMapInputs, configs: tuple[RunConfig, ...] = RUN_CONFIGS
) -> list[PolyfunExplainGroup]:
    return [generate_polyfun_explain_group(shared, c) for c in configs]
```

(Adjust `PolyfunExplainPlotTask.create` to accept `ridge_weights_task` per the Task 3 wiring note.)

- [ ] **Step 3b: The outer generator + shared-setup extraction**

The existing `generate_assets_broad_ukbb_fine_map` does the per-locus shared setup inline (LD interval, renamed labels, harmonize). Extract that setup into a reusable helper `_build_shared_locus_inputs(chrom, pos, build_37_sumstats_task, base_name, sumstats_pipe, sample_size, chrom_range, palindrome_strategy) -> (base_name, harmonized_sumstats_task, ld_labels_task_renamed, ld_matrix_task)` in `fine_mapping_asset_generator.py`, and have the existing generator call it (behavior-preserving refactor; its tests must stay green). Then in the new module:

```python
def generate_assets_polyfun_explain_fine_map(
    chrom: int,
    pos: int,
    build_37_sumstats_task: Task,
    base_name: str,
    sumstats_pipe,
    sample_size_or_effect_sample_size: int,
    gene_info_task: Task,
    q_factor: int = 100,
) -> PolyfunExplainOuterGroup:
    base, harmonized, ld_labels_renamed, ld_matrix = _build_shared_locus_inputs(
        chrom=chrom,
        pos=pos,
        build_37_sumstats_task=build_37_sumstats_task,
        base_name=base_name,
        sumstats_pipe=sumstats_pipe,
        sample_size=sample_size_or_effect_sample_size,
    )
    shared = SharedFineMapInputs(
        base_name=base,
        harmonized_sumstats_task=harmonized,
        ld_labels_task=ld_labels_renamed,
        ld_matrix_task=ld_matrix,
        gene_info_task=gene_info_task,
        effective_sample_size=sample_size_or_effect_sample_size,
        q_factor=q_factor,
    )
    return PolyfunExplainOuterGroup(groups=build_explainability_groups(shared))
```

If extracting `_build_shared_locus_inputs` proves invasive, instead inline the same three steps (LD interval via `get_optimal_ukbb_ld_interval` + `get_ld_labels_and_matrix_task_for_genomic_interval_build_37`, the `PipeDataFrameTask` rename, and `HarmonizeGWASWithReferenceViaAlleles.create`) directly in this function, copying from `generate_assets_broad_ukbb_fine_map`; the outer generator is not unit-tested (it needs real LD-interval data), so correctness is verified in Task 5.

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi r python -m pytest test_mecfs_bio/unit/asset_generator/test_polyfun_explain_fine_mapping_asset_generator.py -v`
Expected: PASS. If the existing generator was refactored, also run its tests: `pixi r python -m pytest test_mecfs_bio/unit/asset_generator/ -v`.

- [ ] **Step 5: invoke green + commit**

```bash
pixi r invoke green > /tmp/green.log 2>&1; tail -5 /tmp/green.log
git add mecfs_bio/asset_generator/ test_mecfs_bio/unit/asset_generator/test_polyfun_explain_fine_mapping_asset_generator.py
git commit -m "feat: polyfun explainability inner+outer asset generators (8 runs/locus)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: Demonstrator locus wiring + verification

**Files:**
- Create: one analysis module under `mecfs_bio/assets/gwas/me_cfs/decode_me/analysis/fine_mapping/` wiring `generate_assets_polyfun_explain_fine_map` at one existing DecodeME polyfun locus (mirror an existing module such as the chr1_174 polyfun one for the sumstats task, gene info task, effective N, chrom/pos).
- Test: a construction test asserting the module builds an 8-run outer group.

**Interfaces:**
- Consumes: the outer generator; an existing DecodeME build-37 sumstats task, `MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW` for gene info.
- Produces: a module-level `POLYFUN_EXPLAIN_<LOCUS> = generate_assets_polyfun_explain_fine_map(...)`.

- [ ] **Step 1: Write the construction test**

```python
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.<demonstrator_module> import (
    POLYFUN_EXPLAIN_CHR1_174,
)


def test_demonstrator_wires_eight_susie_runs():
    assert len(POLYFUN_EXPLAIN_CHR1_174.terminal_tasks()) == 16  # 4*(2 susie+contrast+plot)
    susie = [
        t
        for g in POLYFUN_EXPLAIN_CHR1_174.groups
        for t in (g.susie_uniform, g.susie_polyfun)
    ]
    assert len({t.asset_id for t in susie}) == 8
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi r python -m pytest -k demonstrator_wires -v`
Expected: FAIL — module not created.

- [ ] **Step 3: Create the demonstrator module**

Copy the imports/inputs pattern from the existing chr1_174 polyfun fine-mapping module (sumstats task, sumstats pipe, effective sample size, chrom/pos), and add:

```python
from mecfs_bio.asset_generator.polyfun_explain_fine_mapping_asset_generator import (
    generate_assets_polyfun_explain_fine_map,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.raw.magma_ensembl_gene_location_reference_data_build_37 import (
    MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
)

POLYFUN_EXPLAIN_CHR1_174 = generate_assets_polyfun_explain_fine_map(
    chrom=1,
    pos=174_128_548,
    build_37_sumstats_task=<the existing decode_me build-37 sumstats task>,
    base_name="decode_me_polyfun_explain_",
    sumstats_pipe=<the existing sumstats pipe>,
    sample_size_or_effect_sample_size=<existing effective N>,
    gene_info_task=MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
)
```

- [ ] **Step 4: Run test + green**

Run: `pixi r python -m pytest -k demonstrator_wires -v`
Expected: PASS.

Run: `pixi r invoke green > /tmp/green.log 2>&1; tail -5 /tmp/green.log`
Expected: EXIT 0.

- [ ] **Step 5: Commit**

```bash
git add mecfs_bio/assets/gwas/me_cfs/decode_me/analysis/fine_mapping/ test_mecfs_bio/
git commit -m "feat: demonstrator DecodeME polyfun-explainability locus (chr1_174)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

- [ ] **Step 6: Manual real-data verification (documented, not a unit test)**

This needs the Spec 1 assets built (the 11GB annotation tarball -> annotation parquet -> ridge weights) and the UKBB LD panel. Build one plot task target for the demonstrator locus and inspect:

```bash
# Build the demonstrator plot for L=10 (adjust the target asset_id to the l10 group's plot):
pixi r python -m <the repo's build entrypoint> --target decode_me_polyfun_explain_..._l10_explain_plot
```

Confirm: the run does not raise the prior-coverage guard (or, if it does, that is the real signal to discuss with the user before adding a suppression opt-in); the display table sorts by pip_pf with the family columns present; the figure has 1 + 3 + 3 + 1 = 8 panels and both PNG and SVG are written.

---

## Follow-up (NOT in this plan)

After the demonstrator is verified, roll out to all major DecodeME loci: copy
`mecfs_bio/assets/gwas/me_cfs/decode_me/analysis/fine_mapping/with_palindromes` to a
new sibling folder and rewire each locus module to call
`generate_assets_polyfun_explain_fine_map`. This is instantiation, not new
machinery, and is done as a separate change.

## Self-Review notes (addressed)

- Spec coverage: prior guard + prior.parquet (Task 1); prior lift, per-annotation and per-family contrast, abar_c, focal + family selection, display table with all specified columns/sorting/Int32 (Task 2); 8-panel figure with recomb-from-CM, LD-colored Manhattan, family panels, PNG+SVG (Task 3); inner/outer generators = 8 runs (Task 4); demonstrator + verification (Task 5); rollout deferred (Follow-up).
- Type consistency: display column constants (DISP_*), filenames, and `_KEY` are defined once in the contrast module and imported by the plot task and tests. `RunConfig`/`SharedFineMapInputs`/`PolyfunExplainGroup` names are consistent across generator and tests.
- Resolved seams: `BASELINE_LF_ANNOTATION_MATRIX` lives in `mecfs_bio/assets/reference_data/polyfun/annotations/baseline_lf_annotations.py`; `FilteredGWASDataMeta(id, trait, project, sub_dir, ...)`. One remaining choice, with a concrete in-task fallback: extract `_build_shared_locus_inputs` from the existing generator, or inline the same LD-interval/harmonize steps in the outer generator (the outer generator is verified in Task 5, not unit-tested).
