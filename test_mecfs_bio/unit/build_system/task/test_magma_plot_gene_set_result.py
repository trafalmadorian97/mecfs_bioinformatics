from pathlib import Path

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.magma.magma_plot_gene_set_result import (
    MAGMAPlotGeneSetResult,
)
from mecfs_bio.build_system.wf.base_wf import make_wf


def test_plot_gene_set_result(tmp_path: Path):
    scratch = tmp_path / "scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    plot_task = MAGMAPlotGeneSetResult(
        meta=SimpleDirectoryMeta(AssetId("my_plot")),
        gene_set_analysis_task=FakeTask(SimpleDirectoryMeta(AssetId("gene_set"))),
    )

    def fetch(asset_id: AssetId):
        return DirectoryAsset(
            Path("test_mecfs_bio/unit/dummy_data/fake_gene_set_result")
        )

    result = plot_task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)
