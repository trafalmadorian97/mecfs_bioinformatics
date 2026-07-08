"""
Task to run POPs against MAGMA gene-level results, producing per-gene priority
scores.

This wraps POPs' pops.py. It consumes a munged feature directory (from
PopsMungeFeatureDirectoryTask) and a MAGMA gene-analysis DirectoryAsset (which
contains the .genes.raw and .genes.out files POPs reads via --magma_prefix), and
writes a directory containing <stem>.preds, <stem>.coefs, and <stem>.marginals.
"""

from pathlib import Path, PurePath

import structlog
from attrs import field, frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
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
from mecfs_bio.build_system.task.pops.pops_utils import (
    CONTROL_FEATURES_RELATIVE_PATH,
    GENE_ANNOT_RELATIVE_PATH,
    MUNGED_FEATURES_PREFIX_NAME,
    POPS_OUTPUT_STEM_NAME,
    POPS_SCRIPT_NAME,
    count_feature_chunks,
    invoke_pops_script,
)
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()

DEFAULT_POPS_SUB_DIR = PurePath("analysis") / "pops"


@frozen
class PopsRunTask(Task):
    """Fit the POPs model against MAGMA gene-level results."""

    meta: Meta
    pops_source_task: Task
    munged_features_task: Task
    magma_gene_analysis_task: Task
    magma_stem_name: str = GENE_ANALYSIS_OUTPUT_STEM_NAME
    munged_features_prefix_name: str = MUNGED_FEATURES_PREFIX_NAME
    extra_args: tuple[str, ...] = field(factory=tuple)

    @property
    def deps(self) -> list["Task"]:
        return [
            self.pops_source_task,
            self.munged_features_task,
            self.magma_gene_analysis_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self.pops_source_task.asset_id)
        munged_asset = fetch(self.munged_features_task.asset_id)
        magma_asset = fetch(self.magma_gene_analysis_task.asset_id)
        assert isinstance(source_asset, DirectoryAsset)
        assert isinstance(munged_asset, DirectoryAsset)
        assert isinstance(magma_asset, DirectoryAsset)
        source_dir = source_asset.path
        munged_dir = munged_asset.path

        num_feature_chunks = count_feature_chunks(
            munged_dir, self.munged_features_prefix_name
        )
        gene_annot_path = source_dir / GENE_ANNOT_RELATIVE_PATH
        control_features_path = source_dir / CONTROL_FEATURES_RELATIVE_PATH
        feature_mat_prefix = munged_dir / self.munged_features_prefix_name
        magma_prefix = magma_asset.path / self.magma_stem_name

        out_dir = scratch_dir / "pops_output"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_prefix = out_dir / POPS_OUTPUT_STEM_NAME

        args = [
            "--gene_annot_path",
            str(gene_annot_path),
            "--feature_mat_prefix",
            str(feature_mat_prefix),
            "--num_feature_chunks",
            str(num_feature_chunks),
            "--magma_prefix",
            str(magma_prefix),
            "--control_features_path",
            str(control_features_path),
            "--out_prefix",
            str(out_prefix),
            *self.extra_args,
        ]
        invoke_pops_script(source_dir, POPS_SCRIPT_NAME, args)
        return DirectoryAsset(out_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        pops_source_task: Task,
        munged_features_task: Task,
        magma_gene_analysis_task: Task,
        sub_dir: PurePath = DEFAULT_POPS_SUB_DIR,
        extra_args: tuple[str, ...] = (),
    ) -> "PopsRunTask":
        magma_meta = magma_gene_analysis_task.meta
        assert isinstance(magma_meta, ProcessedGwasDataDirectoryMeta)
        meta = ProcessedGwasDataDirectoryMeta(
            id=AssetId(asset_id),
            trait=magma_meta.trait,
            project=magma_meta.project,
            sub_dir=sub_dir,
        )
        return cls(
            meta=meta,
            pops_source_task=pops_source_task,
            munged_features_task=munged_features_task,
            magma_gene_analysis_task=magma_gene_analysis_task,
            extra_args=tuple(extra_args),
        )
