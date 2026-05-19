from pathlib import Path, PurePath
from typing import Literal

import pandas as pd
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
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameWhiteSpaceSepTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()

GENE_ANALYSIS_OUTPUT_STEM_NAME = "gene_analysis_output"
SynonymMode = Literal["skip", "drop", "drop-dup"]
DuplicateMode = Literal["first", "last", "error"]


@frozen
class MagmaGeneAnalysisTask(Task):
    meta: Meta
    magma_binary_task: Task
    magma_annotation_task: Task
    magma_p_value_task: Task
    magma_ld_ref_task: Task
    ld_ref_file_stem: str
    sample_size: int
    synonym_mode: SynonymMode = "drop-dup"
    duplicate_mode: DuplicateMode | None = "first"

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
            f"duplicate={self.duplicate_mode}",
            f"N={self.sample_size}",
            "--gene-annot",
            str(annotation_path),
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
            id=AssetId(asset_id),
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

    @classmethod
    def create_with_prebuilt_annotation(
        cls,
        asset_id: str,
        magma_annotation_task: Task,
        magma_p_value_task: Task,
        magma_binary_task: Task,
        magma_ld_ref_task: Task,
        ld_ref_file_stem: str,
        sample_size: int,
        sub_dir_suffix: PurePath | None = None,
    ):
        """Create a MagmaGeneAnalysisTask that consumes a pre-built annotation file
        (e.g. an H-MAGMA tissue-specific .genes.annot reference file) instead of
        an annotation produced by an upstream MagmaAnnotateTask. trait/project are
        derived from ``magma_p_value_task.meta`` since the annotation has no
        GWAS-derived metadata. ``sub_dir_suffix`` is appended to the p-value
        task's ``sub_dir`` to isolate the output (e.g. ``PurePath("h_magma")``).
        """
        p_value_meta = magma_p_value_task.meta
        assert isinstance(p_value_meta, FilteredGWASDataMeta)
        sub_dir = PurePath(p_value_meta.sub_dir)
        if sub_dir_suffix is not None:
            sub_dir = sub_dir / sub_dir_suffix
        meta = ProcessedGwasDataDirectoryMeta(
            id=AssetId(asset_id),
            trait=p_value_meta.trait,
            project=p_value_meta.project,
            sub_dir=sub_dir,
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


def read_magma_gene_analysis_result(result_dir: Path) -> pd.DataFrame:
    return (
        scan_dataframe(
            result_dir / str(GENE_ANALYSIS_OUTPUT_STEM_NAME + ".genes.out"),
            spec=DataFrameReadSpec(DataFrameWhiteSpaceSepTextFormat(comment_code="#")),
        )
        .collect()
        .to_pandas()
    )
