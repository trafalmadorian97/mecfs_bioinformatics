import json
from pathlib import Path, PurePath

import numpy as np
import polars as pl
import pytest
from attrs import frozen

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
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    WEIGHTS_PARQUET_FILENAME,
)
from mecfs_bio.build_system.task.base_task import Task
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
    FILTERED_LD_FILENAME,
    PIP_COLUMN,
    PIP_FILENAME,
    PRIOR_FILENAME,
    PRIOR_WEIGHT_COLUMN,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import GWASLAB_BETA_COL, GWASLAB_SE_COL

# Two real baseline-LF annotations from different families so family aggregation
# is exercised: Coding_UCSC_common -> coding, GERP.NS -> conserved.
_ANNOT_A = "Coding_UCSC_common"
_ANNOT_B = "GERP.NS"

_N_VARIANTS = 6
# Diffuse uniform-run PIP used by the closed-form test.
_DEFAULT_UNIFORM_PIP: tuple[float, ...] = (0.2, 0.2, 0.2, 0.2, 0.1, 0.1)
# Polyfun-run PIP is concentrated on variant 0 (focal), for both tests.
_POLYFUN_PIP: tuple[float, ...] = (0.8, 0.05, 0.05, 0.05, 0.03, 0.02)
_PRIOR_WEIGHTS: tuple[float, ...] = (8.0, 1.0, 1.0, 1.0, 1.0, 1.0)
# BETA/SE so the focal variant (row 0) is the Manhattan-panel lead; used by
# the plot task (Task 3), harmless extra columns for the contrast task.
_BETAS: tuple[float, ...] = (2.0, 0.1, 0.1, 0.1, 0.1, 0.1)
_SES: tuple[float, ...] = (0.1, 0.1, 0.1, 0.1, 0.1, 0.1)


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
    gwas = variants.with_columns(
        pl.Series(name=GWASLAB_BETA_COL, values=_BETAS),
        pl.Series(name=GWASLAB_SE_COL, values=_SES),
    )
    gwas.write_parquet(directory / FILTERED_GWAS_FILENAME)
    # Identity LD matrix (only the plot task, Task 3, reads this).
    np.save(directory / FILTERED_LD_FILENAME, np.eye(variants.height))
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
            # A failed uniform run still lists its credible-set variants even
            # when their PIP is 0 - the combined_cs writer never drops rows on
            # PIP value, so the fixture must not either.
            cs_rows.append({**row, CS_COLUMN: label, PIP_COLUMN: float(pip[i])})
    pl.DataFrame(cs_rows).write_parquet(directory / COMBINED_CS_FILENAME)


def _make_contrast_fixture(
    tmp_path: Path, uniform_pip: tuple[float, ...] = _DEFAULT_UNIFORM_PIP
) -> tuple[Path, Path, Path, Path]:
    """Build a two-run (uniform + polyfun), two-annotation locus fixture on disk.

    Returns (uni_dir, pf_dir, weights_dir, annot_path). Shared by every test in
    this module (and available for other tests in this package to build on) so
    the fixture shape stays in one place. weights_dir mirrors the real shape
    RidgeAnnotationWeightsTask.execute produces (a DirectoryAsset containing
    weights.parquet), not a bare file, so the production DirectoryAsset branch
    of _load_weights is what the tests actually exercise.
    """
    assert len(uniform_pip) == _N_VARIANTS
    variants = pl.DataFrame(
        {
            "CHR": [1] * _N_VARIANTS,
            "POS": [10, 20, 30, 40, 50, 60],
            "EA": ["A"] * _N_VARIANTS,
            "NEA": ["C"] * _N_VARIANTS,
        }
    )
    # Annotation values (by CHR/BP). BP == POS.
    a = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # coding: focal only
    b = np.array([2.0, 1.0, 1.0, 1.0, 1.0, 1.0])  # conserved: focal higher
    annot = pl.DataFrame(
        {
            "CHR": [1] * _N_VARIANTS,
            "BP": [10, 20, 30, 40, 50, 60],
            "SNP": [f"rs{i}" for i in range(_N_VARIANTS)],
            "CM": [0.0] * _N_VARIANTS,
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
    weights_dir = tmp_path / "weights"
    weights_dir.mkdir()
    weights.write_parquet(weights_dir / WEIGHTS_PARQUET_FILENAME)

    pip_u = np.array(uniform_pip)
    pip_pf = np.array(_POLYFUN_PIP)
    uni_dir = tmp_path / "uniform"
    pf_dir = tmp_path / "polyfun"
    _write_run_dir(uni_dir, variants, pip_u, {"L1": [0, 1, 2, 3]}, prior_weights=None)
    _write_run_dir(
        pf_dir,
        variants,
        pip_pf,
        {"L1": [0]},
        prior_weights=np.array(_PRIOR_WEIGHTS),
    )
    return uni_dir, pf_dir, weights_dir, annot_path


def _build_contrast_task_and_fetch_map(
    uni_dir: Path,
    pf_dir: Path,
    weights_dir: Path,
    annot_path: Path,
    n_important_families: int = 2,
) -> tuple[PolyfunExplainContrastTask, dict[str, Asset]]:
    """Build a PolyfunExplainContrastTask plus its {asset_id: Asset} fetch map
    over a fixture produced by _make_contrast_fixture. Shared by every test in
    this module and by build_synthetic_explain_inputs."""
    uni_task = FakeTask(ResultDirectoryMeta(id=AssetId("uni"), trait="t", project="p"))
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
        n_important_families=n_important_families,
    )

    fetch_map: dict[str, Asset] = {
        "uni": DirectoryAsset(uni_dir),
        "pf": DirectoryAsset(pf_dir),
        # Mirrors RidgeAnnotationWeightsTask.execute's real return type.
        "weights": DirectoryAsset(weights_dir),
        "annot": FileAsset(annot_path),
    }
    return task, fetch_map


def _run_contrast_task(
    tmp_path: Path,
    uni_dir: Path,
    pf_dir: Path,
    weights_dir: Path,
    annot_path: Path,
    n_important_families: int = 2,
) -> DirectoryAsset:
    """Build and execute a PolyfunExplainContrastTask over a fixture produced by
    _make_contrast_fixture. Shared by every test in this module."""
    task, fetch_map = _build_contrast_task_and_fetch_map(
        uni_dir, pf_dir, weights_dir, annot_path, n_important_families
    )

    def fetch(asset_id: AssetId) -> Asset:
        return fetch_map[str(asset_id)]

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)
    return result


@frozen
class _ExplainInputs:
    """Shared synthetic fixture for the polyfun-explain contrast and plot
    tasks: the task objects (so a caller can wire them as deps into a
    downstream task) plus a fetch map for the contrast task's own dep set and
    the already-executed contrast task directory."""

    uni_task: Task
    pf_task: Task
    weights_task: Task
    annot_task: Task
    contrast_task: Task
    fetch_map: tuple  # tuple of (str asset_id, Asset)


def build_synthetic_explain_inputs(tmp_path: Path) -> _ExplainInputs:
    """Build the two-run, two-annotation synthetic locus fixture, run the
    contrast task once over it, and return the tasks + fetch map + executed
    contrast directory so Task 3's plot task can be wired against the exact
    same data the contrast task's tables were computed from."""
    uni_dir, pf_dir, weights_dir, annot_path = _make_contrast_fixture(tmp_path)
    contrast_task, fetch_map = _build_contrast_task_and_fetch_map(
        uni_dir, pf_dir, weights_dir, annot_path
    )

    def fetch(asset_id: AssetId) -> Asset:
        return fetch_map[str(asset_id)]

    scratch = tmp_path / "contrast_scratch"
    scratch.mkdir()
    contrast_dir = contrast_task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(contrast_dir, DirectoryAsset)

    full_fetch_map = {**fetch_map, "contrast": contrast_dir}
    return _ExplainInputs(
        uni_task=contrast_task.susie_uniform_task,
        pf_task=contrast_task.susie_polyfun_task,
        weights_task=contrast_task.ridge_weights_task,
        annot_task=contrast_task.annotation_parquet_task,
        contrast_task=contrast_task,
        fetch_map=tuple(full_fetch_map.items()),
    )


def test_contrast_closed_form(tmp_path: Path):
    inputs = build_synthetic_explain_inputs(tmp_path)
    result = dict(inputs.fetch_map)["contrast"]

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


def test_contrast_uniform_all_zero_pip_uses_equal_weights(tmp_path: Path):
    # When the uniform run finds no signal (all PIPs 0), abar_c falls back to an
    # unweighted mean, so contrasts are still well-defined (not NaN).
    uni_dir, pf_dir, weights_dir, annot_path = _make_contrast_fixture(
        tmp_path, uniform_pip=(0.0,) * _N_VARIANTS
    )
    result = _run_contrast_task(tmp_path, uni_dir, pf_dir, weights_dir, annot_path)

    per_family = pl.read_parquet(result.path / PER_FAMILY_CONTRAST_FILENAME)
    assert per_family["family_contrast"].is_nan().sum() == 0
    # abar_A now = unweighted mean of a = (1,0,0,0,0,0) over ALL 6 locus
    # variants = 1/6, so focal A contrast = 3*(1 - 1/6) = 2.5.
    focal_coding = per_family.filter(
        (pl.col("POS") == 10) & (pl.col("family") == "coding")
    )["family_contrast"][0]
    assert abs(focal_coding - 2.5) < 1e-9


def test_contrast_raises_on_duplicate_annotation_position(tmp_path: Path):
    # The annotation parquet is deduped only on SNP (rsID), so a multi-allelic
    # site can carry two rows at the same (CHR, BP). Since the annotation join
    # keys on (CHR, POS) only (annotations carry no alleles), an undetected
    # duplicate position would silently cross-multiply variant rows into
    # doubled/misattributed contrast values. The task must fail fast instead.
    uni_dir, pf_dir, weights_dir, annot_path = _make_contrast_fixture(tmp_path)
    annot = pl.read_parquet(annot_path)
    dup_row = annot.filter(pl.col("BP") == 10).with_columns(
        pl.lit("rs0_dup").alias("SNP")
    )
    pl.concat([annot, dup_row], how="vertical").write_parquet(annot_path)

    with pytest.raises(ValueError):
        _run_contrast_task(tmp_path, uni_dir, pf_dir, weights_dir, annot_path)
