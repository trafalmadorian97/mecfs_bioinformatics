"""Unit tests for AcatTask and the ACAT p-value combination kernel."""

from pathlib import Path, PurePath

import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.task.acat_task import ACAT_P_COL, AcatTask, _acat_combine
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import make_wf


def _run_acat_task(
    tmp_path: Path,
    frame: pl.DataFrame,
    excluded_models: tuple[str, ...] = ("syn",),
) -> pl.DataFrame:
    source_meta = FilteredGWASDataMeta(
        id=AssetId("synthetic_phewas"),
        trait="me_cfs",
        project="test",
        sub_dir=PurePath("raw"),
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    )
    source_path = tmp_path / "source.parquet"
    frame.write_parquet(source_path)
    source_task = FakeTask(meta=source_meta)
    assets: dict[AssetId, Asset] = {source_meta.id: FileAsset(source_path)}

    task = AcatTask.create(
        source_task=source_task,
        asset_id="synthetic_acat",
        group_by=["Gene"],
        p_value_col="P value",
        model_col="Collapsing model",
        excluded_models=excluded_models,
    )
    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir(exist_ok=True)
    result = task.execute(
        scratch_dir=scratch_dir,
        fetch=lambda asset_id: assets[asset_id],
        wf=make_wf(),
    )
    assert isinstance(result, FileAsset)
    return pl.read_parquet(result.path)


def test_acat_task_combines_groups_and_preserves_signal_with_unit_and_tiny_pvalues(
    tmp_path: Path,
) -> None:
    frame = pl.DataFrame(
        {
            "Gene": [
                "GENE_UNIT_POLE",
                "GENE_UNIT_POLE",
                "GENE_UNIT_POLE",
                "GENE_UNIT_POLE",
                "GENE_UNIT_POLE",
                "GENE_TINY",
                "GENE_TINY",
                "GENE_TINY",
                "GENE_TINY",
                "GENE_ALL_ONE",
                "GENE_ALL_ONE",
            ],
            "Collapsing model": [
                "ptv",
                "rec",
                "flexdmg",
                "raredmg",
                "syn",
                "ptv",
                "rec",
                "flexdmg",
                "raredmg",
                "ptv",
                "rec",
            ],
            "P value": [
                1e-12,
                1.0,
                0.5,
                0.5,
                0.8,
                1e-20,
                0.5,
                0.5,
                0.5,
                1.0,
                1.0,
            ],
        }
    )

    out = _run_acat_task(tmp_path, frame)
    by_gene = {row["Gene"]: row for row in out.to_dicts()}

    # Excluded synonymous model is omitted from the group summary.
    assert by_gene["GENE_UNIT_POLE"]["n_models"] == 4
    assert by_gene["GENE_UNIT_POLE"]["models_used"] == "flexdmg, ptv, raredmg, rec"
    # A single p=1.0 collapsing model must not cancel a 1e-12 signal to ~1.0.
    assert by_gene["GENE_UNIT_POLE"][ACAT_P_COL] == pytest.approx(
        4e-12, rel=1e-3, abs=0.0
    )

    # Sub-1e-15 signals must not be clamped to ~1e-15 or underflow to 0.0.
    assert by_gene["GENE_TINY"][ACAT_P_COL] == pytest.approx(4e-20, rel=1e-3, abs=0.0)

    # All p=1.0 models yield 1.0.
    assert by_gene["GENE_ALL_ONE"][ACAT_P_COL] == pytest.approx(1.0)


def test_acat_combine_numerical_precision_and_boundaries() -> None:
    # Single sub-1e-15 p-values round-trip accurately without clamping or cancellation to 0.0.
    assert _acat_combine(pl.Series([1e-50])) == pytest.approx(1e-50, rel=1e-6, abs=0.0)
    assert _acat_combine(pl.Series([1e-200])) == pytest.approx(
        1e-200, rel=1e-6, abs=0.0
    )
    assert _acat_combine(pl.Series([1e-50])) < _acat_combine(pl.Series([1e-20]))

    # Exact zero returns 0.0; exact single one returns 1.0.
    assert _acat_combine(pl.Series([0.0, 0.5])) == 0.0
    assert _acat_combine(pl.Series([1.0])) == 1.0

    # Ordinary interior p-values match the closed-form Cauchy median.
    assert _acat_combine(pl.Series([0.5, 0.5])) == pytest.approx(0.5, abs=1e-15)
    assert _acat_combine(pl.Series([0.25])) == pytest.approx(0.25, abs=1e-15)


def test_acat_combine_rejects_out_of_range_pvalues() -> None:
    with pytest.raises(AssertionError):
        _acat_combine(pl.Series([-0.01, 0.5]))
    with pytest.raises(AssertionError):
        _acat_combine(pl.Series([0.5, 1.01]))
    with pytest.raises(AssertionError):
        _acat_combine(pl.Series([float("nan"), 0.5]))
