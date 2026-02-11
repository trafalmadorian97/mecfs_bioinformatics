from pathlib import Path

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
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF

#
# _dummy_gene_analysis_data = """# VERSION = 110
# # COVAR = NSAMP MAC
# 148398 1 824993 889961 211 28 275488 136.915 1.59826
# 339451 1 860967 911099 151 29 275488 95.3311 0.754161 0.590117
# 84069 1 866872 920488 158 29 275488 118.139 0.747723 0.479151 0.879137
# 26155 1 869583 929679 184 32 275488 111.136 0.770145 0.378889 0.861104 0.922526
# 84808 1 900579 952473 171 22 275488 151.38 0.789199 0.164544 0.312259 0.592093 0.750633
# 9636 1 913847 959920 156 18 275488 172.744 0.454433 0.113036 0.210681 0.427217 0.686265 0.955037
# 375790 1 920503 1001499 226 25 275488 160.27 -0.507389 0.100859 0.202208 0.3747 0.592698 0.852268 0.970024
# 57801 1 924342 971608 136 16 275488 182.522 -0.201964 0.0907966 0.18144 0.355704 0.590505 0.896288 1 0.93935
# """
#
_dummy_spec_data = """GENE
1
10
100"""


def test_gene_set_analysis_empty_gene_sets(tmp_path: Path):
    """
    Verify that we can gracefully handle the case in which we are asked to perform joint gene-set analysis without any gene sets
    This situation arises in MAGMA when no gene sets are significant
    """
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    dummy_spec_matrix_path = tmp_path / "dummy_spec_matrix.txt"
    dummy_spec_matrix_path.write_text(_dummy_spec_data)
    # gene_set_analysis_task= FakeTask(
    #     SimpleDirectoryMeta("GeneSetAnalysisTask"),
    # )
    gene_set_covar_analysis_task = MagmaGeneSetAnalysisTask(
        meta=SimpleDirectoryMeta("GeneSetAnalysis"),
        magma_binary_task=FakeTask(meta=SimpleFileMeta("MagmaBinary")),
        magma_gene_analysis_task=FakeTask(meta=SimpleFileMeta("MagmaGeneAnalysis")),
        gene_set_or_covar_task=ExternalFileCopyTask(
            SimpleFileMeta("gene_sets_task"), external_path=dummy_spec_matrix_path
        ),
        set_or_covar="covar",
        model_params=None,
    )

    # stepwise_task=MagmaForwardStepwiseSelectTask(
    #     meta=SimpleFileMeta("MagmaForwardStepwiseSelect"),
    #     magma_marginal_output_task=gene_set_covar_analysis_task,
    #     magma_conditional_output_task=gene_set_analysis_task,
    # )
    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "gene_sets_task":
            return FileAsset(dummy_spec_matrix_path)
        else:
            raise ValueError("Unknown asset id")

    result = gene_set_covar_analysis_task.execute(
        scratch_dir=scratch, fetch=fetch, wf=SimpleWF()
    )
    assert isinstance(result, DirectoryAsset)
