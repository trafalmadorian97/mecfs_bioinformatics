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

_N_VARIANTS = 6
# Diffuse uniform-run PIP used by the closed-form test.
_DEFAULT_UNIFORM_PIP: tuple[float, ...] = (0.2, 0.2, 0.2, 0.2, 0.1, 0.1)
# Polyfun-run PIP is concentrated on variant 0 (focal), for both tests.
_POLYFUN_PIP: tuple[float, ...] = (0.8, 0.05, 0.05, 0.05, 0.03, 0.02)
_PRIOR_WEIGHTS: tuple[float, ...] = (8.0, 1.0, 1.0, 1.0, 1.0, 1.0)


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
            # A failed uniform run still lists its credible-set variants even
            # when their PIP is 0 - the combined_cs writer never drops rows on
            # PIP value, so the fixture must not either.
            cs_rows.append({**row, CS_COLUMN: label, PIP_COLUMN: float(pip[i])})
    pl.DataFrame(cs_rows).write_parquet(directory / COMBINED_CS_FILENAME)


def _make_contrast_fixture(
    tmp_path: Path, uniform_pip: tuple[float, ...] = _DEFAULT_UNIFORM_PIP
) -> tuple[Path, Path, Path, Path]:
    """Build a two-run (uniform + polyfun), two-annotation locus fixture on disk.

    Returns (uni_dir, pf_dir, weights_path, annot_path). Shared by every test in
    this module (and available for other tests in this package to build on) so
    the fixture shape stays in one place.
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
    weights_path = tmp_path / "weights.parquet"
    weights.write_parquet(weights_path)

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
    return uni_dir, pf_dir, weights_path, annot_path


def _run_contrast_task(
    tmp_path: Path,
    uni_dir: Path,
    pf_dir: Path,
    weights_path: Path,
    annot_path: Path,
    n_important_families: int = 2,
) -> DirectoryAsset:
    """Build and execute a PolyfunExplainContrastTask over a fixture produced by
    _make_contrast_fixture. Shared by every test in this module."""
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
    return result


def test_contrast_closed_form(tmp_path: Path):
    uni_dir, pf_dir, weights_path, annot_path = _make_contrast_fixture(tmp_path)
    result = _run_contrast_task(tmp_path, uni_dir, pf_dir, weights_path, annot_path)

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
    uni_dir, pf_dir, weights_path, annot_path = _make_contrast_fixture(
        tmp_path, uniform_pip=(0.0,) * _N_VARIANTS
    )
    result = _run_contrast_task(tmp_path, uni_dir, pf_dir, weights_path, annot_path)

    per_family = pl.read_parquet(result.path / PER_FAMILY_CONTRAST_FILENAME)
    assert per_family["family_contrast"].is_nan().sum() == 0
    # abar_A now = unweighted mean of a = (1,0,0,0,0,0) over ALL 6 locus
    # variants = 1/6, so focal A contrast = 3*(1 - 1/6) = 2.5.
    focal_coding = per_family.filter(
        (pl.col("POS") == 10) & (pl.col("family") == "coding")
    )["family_contrast"][0]
    assert abs(focal_coding - 2.5) < 1e-9
