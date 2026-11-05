from pathlib import Path

import polars as pl

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.protein_coloc.coloc_display_table_task import (
    ColocDisplayTableTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.protein_coloc_constants import (
    DISPLAY_LEAD_VARIANT_COL,
    DISPLAY_LOCUS_COL,
    GENE_COL,
    PP_H4_COL,
)
from test_mecfs_bio.unit.build_system.task.protein_coloc.coloc_result_fixtures import (
    COLOC_ID,
    coloc_task,
    write_coloc_dir,
)


def _display(tmp_path: Path, h4_threshold: float) -> pl.DataFrame:
    coloc_dir = write_coloc_dir(tmp_path / "coloc")
    task = ColocDisplayTableTask.create(
        asset_id="display", coloc_task=coloc_task(), h4_threshold=h4_threshold
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    assets = {AssetId(COLOC_ID): DirectoryAsset(coloc_dir)}
    result = task.execute(
        scratch, fetch=lambda asset_id: assets[asset_id], wf=make_wf()
    )
    assert isinstance(result, FileAsset)
    return pl.read_parquet(result.path)


def test_display_keeps_only_rows_above_threshold_sorted(tmp_path: Path) -> None:
    display = _display(tmp_path, h4_threshold=0.8)
    assert display[GENE_COL].to_list() == ["AAA", "BBB"]
    assert display[PP_H4_COL].to_list() == [0.95, 0.85]
    assert display[DISPLAY_LOCUS_COL].to_list() == ["chr1:1000-3000"] * 2
    assert display[DISPLAY_LEAD_VARIANT_COL].to_list() == ["chr1:2000:G:A"] * 2


def test_display_integer_columns_are_int32(tmp_path: Path) -> None:
    display = _display(tmp_path, h4_threshold=0.8)
    assert pl.Int64 not in display.schema.dtypes()


def test_empty_display_when_nothing_passes(tmp_path: Path) -> None:
    display = _display(tmp_path, h4_threshold=0.99)
    assert display.height == 0
    assert GENE_COL in display.columns and PP_H4_COL in display.columns
