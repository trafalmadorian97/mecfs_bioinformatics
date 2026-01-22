from pathlib import Path, PurePath
from typing import Literal, Sequence

import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.processed_gwas_data_directory_meta import (
    ProcessedGwasDataDirectoryMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.magma.magma_gene_analysis_task import (
    GENE_ANALYSIS_OUTPUT_STEM_NAME,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()

GENE_SET_ANALYSIS_OUTPUT_STEM_NAME = "gene_set_analysis_output"


@frozen
class DirectoryGeneSetSpec:
    gene_set_task: Task
    path_in_dir: PurePath


SetOrCovar = Literal["set", "covar"]


@frozen()
class ModelParams:
    direction_covar: str | None
    condition_hide: Sequence[str]
    joint_pairs: bool = False

    def __attrs_post_init__(self):
        assert not isinstance(self.condition_hide, str)

    def prep_command(self) -> list[str]:
        result = ["--model"]
        if self.joint_pairs:
            result += ["joint-pairs"]
        if self.direction_covar:
            result += [f"direction-covar={self.direction_covar}"]
        if len(self.condition_hide) > 0:
            result += ["condition-hide=" + ",".join(self.condition_hide)]

        return result


@frozen
class MagmaGeneSetAnalysisTask(Task):
    """
    The final step in the canonical MAGMA pipeline.
    See page 18 of the manual here: https://vu.data.surfsara.nl/s/MUiv3y1SFRePnyG?dir=/&editing=false&openfile=true
    """

    _meta: Meta
    magma_binary_task: Task
    magma_gene_analysis_task: Task
    gene_set_or_covar_task: Task | DirectoryGeneSetSpec
    set_or_covar: SetOrCovar
    model_params: ModelParams | None

    @property
    def _magma_binary_id(self) -> AssetId:
        return self.magma_binary_task.asset_id

    def _magma_gene_analysis_id(self) -> AssetId:
        return self.magma_gene_analysis_task.asset_id

    @property
    def _gene_set_or_covar_task_id(self) -> AssetId:
        if isinstance(self.gene_set_or_covar_task, Task):
            return self.gene_set_or_covar_task.asset_id
        return self.gene_set_or_covar_task.gene_set_task.asset_id

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        deps = [self.magma_binary_task, self.magma_gene_analysis_task]
        if isinstance(self.gene_set_or_covar_task, Task):
            deps.append(self.gene_set_or_covar_task)
        else:
            deps.append(self.gene_set_or_covar_task.gene_set_task)
        return deps

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        binary_asset = fetch(self._magma_binary_id)
        gene_analysis_asset = fetch(self._magma_gene_analysis_id())
        gene_set_or_covar_asset = fetch(self._gene_set_or_covar_task_id)
        assert isinstance(binary_asset, FileAsset)
        assert isinstance(gene_analysis_asset, DirectoryAsset)
        if isinstance(self.gene_set_or_covar_task, Task):
            assert isinstance(gene_set_or_covar_asset, FileAsset)
            gene_set_or_covar_path = gene_set_or_covar_asset.path
        else:
            assert isinstance(self.gene_set_or_covar_task, DirectoryGeneSetSpec)
            assert isinstance(gene_set_or_covar_asset, DirectoryAsset)
            gene_set_or_covar_path = (
                gene_set_or_covar_asset.path / self.gene_set_or_covar_task.path_in_dir
            )

        binary_path = binary_asset.path
        gene_analysis_path_root_path = gene_analysis_asset.path
        gene_analysis_full_path = gene_analysis_path_root_path / (
            GENE_ANALYSIS_OUTPUT_STEM_NAME + ".genes.raw"
        )
        out_dir = scratch_dir / "gene_set_analysis_dir"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_base_path = out_dir / GENE_SET_ANALYSIS_OUTPUT_STEM_NAME

        set_or_covar_command = (
            "--set-annot" if self.set_or_covar == "set" else "--gene-covar"
        )

        cmd = [
            str(binary_path),
            "--gene-results",
            str(gene_analysis_full_path),
            set_or_covar_command,  # "--set-annot"  or "--gene-covar"
            str(gene_set_or_covar_path),
        ]
        if self.model_params is not None:
            cmd.extend(self.model_params.prep_command())
        cmd += [
            "--out",
            str(out_base_path),
        ]
        logger.debug(f"Running command: {' '.join(cmd)}")
        execute_command(cmd)
        return DirectoryAsset(out_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        magma_gene_analysis_task: Task,
        magma_binary_task: Task,
        gene_set_task: Task | DirectoryGeneSetSpec,
        set_or_covar: SetOrCovar,
        model_params: ModelParams | None,
        sub_dir: PurePath = PurePath("analysis_results") / "magma",
    ):
        gene_analysis_meta = magma_gene_analysis_task.meta
        assert isinstance(gene_analysis_meta, ProcessedGwasDataDirectoryMeta)
        meta = ProcessedGwasDataDirectoryMeta(
            short_id=AssetId(asset_id),
            trait=gene_analysis_meta.trait,
            project=gene_analysis_meta.project,
            sub_dir=sub_dir,
        )
        return cls(
            magma_binary_task=magma_binary_task,
            gene_set_or_covar_task=gene_set_task,
            magma_gene_analysis_task=magma_gene_analysis_task,
            meta=meta,
            set_or_covar=set_or_covar,
            model_params=model_params,
        )
