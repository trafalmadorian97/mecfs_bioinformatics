import pickle
from pathlib import Path

import gwaslab
import gwaslab as gl
import pandas as pd
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
    _validate_eaf_in_range,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF
from mecfs_bio.constants.gwaslab_constants import GWASLAB_STATUS_COL


def test_gwaslab_sumstats(
    tmp_path: Path,
):
    scratch_loc = tmp_path / "scratch"
    scratch_loc.mkdir(exist_ok=True, parents=True)
    task = GWASLabCreateSumstatsTask(
        target_asset_id=AssetId("sumstats_task"),
        basic_check=True,
        df_source_task=FakeTask(
            meta=FilteredGWASDataMeta(
                AssetId("input"),
                "dummy_trait",
                "dummy_project",
                "dummy_Dir",
                read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=" ")),
            )
        ),
        genome_build="38",
    )

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(
            Path("test_mecfs_bio/unit/build_system/task/dummy_data.regenie")
        )

    asset_result = task.execute(scratch_dir=scratch_loc, fetch=fetch, wf=SimpleWF())
    with open(asset_result.path, "rb") as f:
        loaded = pickle.load(
            f,
        )
    assert isinstance(loaded, gwaslab.Sumstats)
    assert loaded.data[GWASLAB_STATUS_COL].dtype == pd.Int64Dtype(), (
        f"Expected gwaslab STATUS column to be Int64, got {loaded.data[GWASLAB_STATUS_COL].dtype}. "
        "If gwaslab changed this dtype, update the STATUS column handling in "
        "gwaslab_create_sumstats_task.py (_do_harmonization and _sumstats_raise_on_error)."
    )
    scratch_loc_2 = tmp_path / "scratch_2"
    scratch_loc_2.mkdir(exist_ok=True, parents=True)

    def fetch_2(asset_id: AssetId) -> Asset:
        return asset_result

    task_2 = GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=task, asset_id="table_task", sub_dir="processed"
    )
    asset_2 = task_2.execute(scratch_dir=scratch_loc_2, wf=SimpleWF(), fetch=fetch_2)
    scan_dataframe_asset(asset=asset_2, meta=task_2.meta)


def _regenie_sumstats_with_a1freq(a1freq: list[float]) -> gl.Sumstats:
    n = len(a1freq)
    df = pd.DataFrame(
        {
            "CHROM": [1] * n,
            "GENPOS": list(range(100, 100 + n)),
            "ID": [f"rs{i}" for i in range(n)],
            "ALLELE0": ["A"] * n,
            "ALLELE1": ["T"] * n,
            "A1FREQ": a1freq,
            "BETA": [0.1] * n,
            "SE": [0.01] * n,
            "LOG10P": [2.0] * n,
            "N": [1000] * n,
        }
    )
    return gl.Sumstats(df, fmt="regenie", verbose=False)


def test_validate_eaf_in_range_rejects_percentage():
    """
    EAF supplied as a percentage (0-100) must fail fast, even via a named format where
    gwaslab renames the source column (A1FREQ) to its standard EAF column.
    """
    sumstats = _regenie_sumstats_with_a1freq([0.1, 25.0, 50.0])
    with pytest.raises(AssertionError, match="percentages"):
        _validate_eaf_in_range(sumstats)


def test_validate_eaf_in_range_accepts_fraction():
    sumstats = _regenie_sumstats_with_a1freq([0.001, 0.25, 0.5, 0.999])
    _validate_eaf_in_range(sumstats)


def test_validate_eaf_in_range_no_eaf_column_is_noop():
    """A sumstats with no EAF column must not raise."""
    df = pd.DataFrame(
        {
            "CHROM": [1, 1],
            "GENPOS": [100, 200],
            "ID": ["rs1", "rs2"],
            "ALLELE0": ["A", "C"],
            "ALLELE1": ["T", "G"],
            "BETA": [0.1, 0.1],
            "SE": [0.01, 0.01],
            "LOG10P": [2.0, 3.0],
            "N": [1000, 1000],
        }
    )
    sumstats = gl.Sumstats(df, fmt="regenie", verbose=False)
    _validate_eaf_in_range(sumstats)
