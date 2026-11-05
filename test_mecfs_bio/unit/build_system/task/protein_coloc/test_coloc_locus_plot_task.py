from pathlib import Path, PurePath

import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.protein_coloc.coloc_display_table_task import (
    ColocDisplayTableTask,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_locus_plot_task import (
    PLOT_INDEX_FILENAME,
    ColocLocusPlotTask,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
    GENE_INFO_STRAND_COL,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.protein_coloc_constants import (
    PLOT_PNG_FILENAME_COL,
    PLOT_SVG_FILENAME_COL,
    PROTEIN_ID_COL,
)
from test_mecfs_bio.unit.build_system.task.protein_coloc.coloc_result_fixtures import (
    COLOC_ID,
    coloc_task,
    write_coloc_dir,
)


def _gene_info_task() -> FakeTask:
    return FakeTask(
        ResultTableMeta(
            id=AssetId("genes"),
            trait="synthetic",
            project="test",
            sub_dir=PurePath("reference"),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
    )


def _plots(tmp_path: Path, h4_threshold: float) -> Path:
    coloc_dir = write_coloc_dir(tmp_path / "coloc")
    genes_path = tmp_path / "genes.parquet"
    pl.DataFrame(
        {
            GENE_INFO_CHROM_COL: ["1", "2"],
            GENE_INFO_START_COL: [1500, 1500],
            GENE_INFO_END_COL: [2500, 2500],
            GENE_INFO_STRAND_COL: ["+", "-"],
            GENE_INFO_NAME_COL: ["AAA", "CCC"],
        }
    ).write_parquet(genes_path)
    display = ColocDisplayTableTask.create(
        asset_id="display", coloc_task=coloc_task(), h4_threshold=h4_threshold
    )
    task = ColocLocusPlotTask.create(
        asset_id="plots",
        coloc_task=coloc_task(),
        display_task=display,
        gene_info_task=_gene_info_task(),
    )
    assets: dict[AssetId, Asset] = {
        AssetId(COLOC_ID): DirectoryAsset(coloc_dir),
        AssetId("genes"): FileAsset(genes_path),
    }
    display_scratch = tmp_path / "display_scratch"
    display_scratch.mkdir()
    assets[display.asset_id] = display.execute(
        display_scratch, fetch=lambda asset_id: assets[asset_id], wf=make_wf()
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(
        scratch, fetch=lambda asset_id: assets[asset_id], wf=make_wf()
    )
    assert isinstance(result, DirectoryAsset)
    return result.path


def test_plot_files_and_index_for_each_display_row(tmp_path: Path) -> None:
    out_dir = _plots(tmp_path, h4_threshold=0.8)
    index = pl.read_parquet(out_dir / PLOT_INDEX_FILENAME)
    assert sorted(index[PROTEIN_ID_COL].to_list()) == ["OID1", "OID2"]
    for png, svg in index.select(PLOT_PNG_FILENAME_COL, PLOT_SVG_FILENAME_COL).rows():
        assert (out_dir / png).stat().st_size > 0
        assert (out_dir / svg).stat().st_size > 0


def test_no_display_rows_gives_empty_index_and_no_images(tmp_path: Path) -> None:
    out_dir = _plots(tmp_path, h4_threshold=0.99)
    assert pl.read_parquet(out_dir / PLOT_INDEX_FILENAME).height == 0
    assert not any(out_dir.glob("*.png"))
