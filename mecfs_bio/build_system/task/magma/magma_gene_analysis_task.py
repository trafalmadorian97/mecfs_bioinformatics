from pathlib import Path, PurePath
from typing import Literal

import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.processed_gwas_data_directory_meta import (
    ProcessedGwasDataDirectoryMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()

GENE_ANALYSIS_OUTPUT_STEM_NAME = "gene_analysis_output"
SynonymMode = Literal["skip", "drop", "drop-dup"]


@frozen
class MagmaGeneAnalysisTask(Task):
    _meta: Meta
    magma_binary_task: Task
    magma_annotation_task: Task
    magma_p_value_task: Task
    magma_ld_ref_task: Task
    ld_ref_file_stem: str
    sample_size: int
    synonym_mode: SynonymMode = "drop-dup"

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def p_value_meta(self) -> Meta:
        return self.magma_p_value_task.meta

    @property
    def p_value_id(self) -> AssetId:
        return self.magma_p_value_task.asset_id

    @property
    def magma_binary_id(self) -> AssetId:
        return self.magma_binary_task.asset_id

    @property
    def magma_annotation_id(self) -> AssetId:
        return self.magma_annotation_task.asset_id

    @property
    def magma_ld_ref_id(self) -> AssetId:
        return self.magma_ld_ref_task.asset_id

    @property
    def deps(self) -> list["Task"]:
        return [
            self.magma_binary_task,
            self.magma_annotation_task,
            self.magma_p_value_task,
            self.magma_ld_ref_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        binary_asset = fetch(self.magma_binary_id)
        annotation_asset = fetch(self.magma_annotation_id)
        p_value_asset = fetch(self.p_value_id)
        ld_ref_asset = fetch(self.magma_ld_ref_id)
        assert isinstance(binary_asset, FileAsset)
        assert isinstance(p_value_asset, FileAsset)
        assert isinstance(annotation_asset, FileAsset)
        assert isinstance(ld_ref_asset, DirectoryAsset)
        binary_path = binary_asset.path
        annotation_path = annotation_asset.path
        p_value_path = p_value_asset.path
        ld_ref_dir_path = ld_ref_asset.path
        out_dir = scratch_dir / "gene_analysis_dir"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_base_path = out_dir / GENE_ANALYSIS_OUTPUT_STEM_NAME
        cmd = [
            str(binary_path),
            "--bfile",
            str(ld_ref_dir_path / self.ld_ref_file_stem),
            f"synonym-dup={self.synonym_mode}",
            "--pval",
            str(p_value_path),
            f"N={self.sample_size}",
            "--gene-annot",
            str(annotation_path),
            "--out",
            str(out_base_path),
        ]
        logger.debug(f"Running command: {' '.join(cmd)}")
        # result = subprocess.run(cmd, check=True, text=True, capture_output=True)
        result = execute_command(cmd)
        # logger.debug(
        #     f"Command produced result: \n\n{result.stdout}\n{result.stderr}\n\n"
        # )
        # out_full_path = Path(str(out_base_path) + ".genes.raw")
        return DirectoryAsset(out_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        magma_annotation_task: Task,
        magma_p_value_task: Task,
        magma_binary_task: Task,
        magma_ld_ref_task: Task,
        ld_ref_file_stem: str,
        sample_size: int,
    ):
        annotation_meta = magma_annotation_task.meta  # magma_p_value_task.meta
        assert isinstance(annotation_meta, FilteredGWASDataMeta)
        meta = ProcessedGwasDataDirectoryMeta(
            short_id=AssetId(asset_id),
            trait=annotation_meta.trait,
            project=annotation_meta.project,
            sub_dir=PurePath(annotation_meta.sub_dir),
        )
        return cls(
            magma_annotation_task=magma_annotation_task,
            magma_p_value_task=magma_p_value_task,
            magma_binary_task=magma_binary_task,
            magma_ld_ref_task=magma_ld_ref_task,
            ld_ref_file_stem=ld_ref_file_stem,
            sample_size=sample_size,
            meta=meta,
        )
