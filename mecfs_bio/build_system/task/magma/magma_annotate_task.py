from pathlib import Path, PurePath

import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()


@frozen
class MagmaAnnotateTask(Task):
    """
    Perform the annotate step of the MAGMA pipeline.
    This step associates SNPs with genes.
    See page 5 of the MAGMA manual here: https://vu.data.surfsara.nl/s/MUiv3y1SFRePnyG?dir=/&editing=false&openfile=true

    The window option is used to optionally associate SNPs outside the transcription region of a gene with the gene.
    """

    _meta: Meta
    magma_binary_task: Task
    snp_loc_file_task: Task
    gene_loc_file_task: Task
    window: tuple[int, int] | None = None

    @property
    def magma_binary_id(self) -> AssetId:
        return self.magma_binary_task.asset_id

    @property
    def snp_loc_id(self) -> AssetId:
        return self.snp_loc_file_task.asset_id

    @property
    def gene_loc_id(self) -> AssetId:
        return self.gene_loc_file_task.asset_id

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.magma_binary_task, self.snp_loc_file_task, self.gene_loc_file_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        binary_asset = fetch(self.magma_binary_id)
        snp_loc_asset = fetch(self.snp_loc_id)
        gene_loc_asset = fetch(self.gene_loc_id)
        assert isinstance(binary_asset, FileAsset)
        assert isinstance(snp_loc_asset, FileAsset)
        assert isinstance(gene_loc_asset, FileAsset)
        binary_path = binary_asset.path
        snp_loc_path = snp_loc_asset.path
        gene_loc_path = gene_loc_asset.path
        out_base_path = scratch_dir / "out_base"
        cmd = [
            str(binary_path),
            "--annotate",
        ]
        if self.window is not None:
            cmd.extend([f"window={self.window[0]},{self.window[1]}"])
        cmd += [
            "--snp-loc",
            str(snp_loc_path),
            "--gene-loc",
            str(gene_loc_path),
            "--out",
            str(out_base_path),
        ]
        logger.debug(f"Running command: {' '.join(cmd)}")
        result = execute_command(cmd)
        logger.debug(f"Command produced result: \n\n{result}\n\n\n")
        out_full_path = Path(str(out_base_path) + ".genes.annot")
        return FileAsset(out_full_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        magma_binary_task: Task,
        snp_loc_file_task: Task,
        gene_loc_file_task: Task,
        window: tuple[int, int] | None = None,
    ) -> "MagmaAnnotateTask":
        snp_loc_meta = snp_loc_file_task.meta
        assert isinstance(snp_loc_meta, FilteredGWASDataMeta)
        meta = FilteredGWASDataMeta(
            id=AssetId(asset_id),
            trait=snp_loc_meta.trait,
            project=snp_loc_meta.project,
            sub_dir=PurePath(snp_loc_meta.sub_dir) / "magma",
            extension=".genes.annot",
        )
        return cls(
            meta=meta,
            magma_binary_task=magma_binary_task,
            snp_loc_file_task=snp_loc_file_task,
            gene_loc_file_task=gene_loc_file_task,
            window=window,
        )
