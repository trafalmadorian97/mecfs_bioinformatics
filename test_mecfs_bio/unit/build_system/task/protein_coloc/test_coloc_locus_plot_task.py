from pathlib import Path, PurePath

import polars as pl
import pytest

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
from mecfs_bio.build_system.task.protein_coloc.coloc_locus_plot_task import (
    PLOT_INDEX_FILENAME,
    ColocLocusPlotTask,
    ColocPair,
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


def _plots(tmp_path: Path, pairs: tuple[ColocPair, ...]) -> Path:
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
    task = ColocLocusPlotTask.create(
        asset_id="plots",
        coloc_task=coloc_task(),
        gene_info_task=_gene_info_task(),
        pairs=pairs,
    )
    assets: dict[AssetId, Asset] = {
        AssetId(COLOC_ID): DirectoryAsset(coloc_dir),
        AssetId("genes"): FileAsset(genes_path),
    }
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(
        scratch, fetch=lambda asset_id: assets[asset_id], wf=make_wf()
    )
    assert isinstance(result, DirectoryAsset)
    return result.path


def test_plot_files_and_index_for_each_requested_pair(tmp_path: Path) -> None:
    # OID3 has a low PP.H4: requested pairs are plotted whatever their posterior.
    out_dir = _plots(
        tmp_path,
        pairs=(
            ColocPair(locus_id="chr1:1000-3000", protein_id="OID1"),
            ColocPair(locus_id="chr2:1000-3000", protein_id="OID3"),
        ),
    )
    index = pl.read_parquet(out_dir / PLOT_INDEX_FILENAME)
    assert sorted(index[PROTEIN_ID_COL].to_list()) == ["OID1", "OID3"]
    for png, svg in index.select(PLOT_PNG_FILENAME_COL, PLOT_SVG_FILENAME_COL).rows():
        assert (out_dir / png).stat().st_size > 0
        assert (out_dir / svg).stat().st_size > 0


def test_pair_absent_from_results_asserts(tmp_path: Path) -> None:
    # OID1 was tested, but at the chr1 locus, not the chr2 one.
    with pytest.raises(AssertionError):
        _plots(
            tmp_path, pairs=(ColocPair(locus_id="chr2:1000-3000", protein_id="OID1"),)
        )


def test_pair_not_tested_asserts(tmp_path: Path) -> None:
    # OID4 was dropped by the per-locus budget, so it has no aligned window.
    with pytest.raises(AssertionError):
        _plots(
            tmp_path, pairs=(ColocPair(locus_id="chr2:1000-3000", protein_id="OID4"),)
        )


def test_no_pairs_asserts() -> None:
    with pytest.raises(AssertionError):
        ColocLocusPlotTask.create(
            asset_id="plots",
            coloc_task=coloc_task(),
            gene_info_task=_gene_info_task(),
            pairs=(),
        )
