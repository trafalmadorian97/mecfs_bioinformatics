from pathlib import Path

import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.external_file_copy_task import ExternalFileCopyTask
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import (
    MagmaGeneSetAnalysisTask,
    ModelParams,
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF

_dummy_spec_data_1 = """GENE
1
10
100"""


_dummy_spec_data_2 = """GENE  Cluster132
1   0.002768
10  0.004922
100 0.002771 """


@pytest.mark.parametrize(
    "dummy_spec_data,model_params",
    [
        [_dummy_spec_data_1, None],
        [
            _dummy_spec_data_2,
            ModelParams(direction_covar="greater", condition_hide=[], joint_pairs=True),
        ],
    ],
)
def test_gene_set_analysis_empty_gene_sets(
    tmp_path: Path, dummy_spec_data: str, model_params: ModelParams | None
):
    """
    Verify that we can gracefully handle the case in which we are asked to perform joint gene-set analysis without any gene sets
    This situation arises in MAGMA when no gene sets are significant
    """
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    dummy_spec_matrix_path = tmp_path / "dummy_spec_matrix.txt"
    dummy_spec_matrix_path.write_text(dummy_spec_data)
    gene_set_covar_analysis_task = MagmaGeneSetAnalysisTask(
        meta=SimpleDirectoryMeta("GeneSetAnalysis"),
        magma_binary_task=FakeTask(meta=SimpleFileMeta("MagmaBinary")),
        magma_gene_analysis_task=FakeTask(meta=SimpleFileMeta("MagmaGeneAnalysis")),
        gene_set_or_covar_task=ExternalFileCopyTask(
            SimpleFileMeta("gene_sets_task"), external_path=dummy_spec_matrix_path
        ),
        set_or_covar="covar",
        model_params=model_params,
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "gene_sets_task":
            return FileAsset(dummy_spec_matrix_path)
        else:
            raise ValueError(
                f"Test failure!  Task tried to access asset: {asset_id}.  It should have terminated before this point"
            )

    result = gene_set_covar_analysis_task.execute(
        scratch_dir=scratch, fetch=fetch, wf=SimpleWF()
    )
    assert isinstance(result, DirectoryAsset)
