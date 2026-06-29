"""
Task that trains polypwas protein weights from pQTL summary statistics.

polypwas train builds per-protein SBayesRC weights; it shells out to SBayesRC via
Rscript (see polypwas_common.polypwas_setup).  Reference:
Functionally informed cis and trans proteome-wide association studies prioritize
disease-critical genes (https://www.medrxiv.org/content/10.64898/2026.04.24.26351667v1).
"""

from pathlib import Path, PurePath

from attrs import frozen
from structlog import get_logger

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.polypwas.polypwas_common import (
    derive_trait_project,
    polypwas_setup,
    run_polypwas,
)
from mecfs_bio.build_system.task.sbayesrc.sbayesrc_data_source import (
    SBayesRCSource,
    prepare_cojo_ma_input_file,
)
from mecfs_bio.build_system.task.sbayesrc.sbayesrc_utils import (
    write_docker_rscript_shim,
)
from mecfs_bio.build_system.wf.base_wf import WF

logger = get_logger()

POLYPWAS_WEIGHTS_EXTENSION = ".wgts.gz"


@frozen
class PolypwasTrainTask(Task):
    """
    Train polypwas protein weights for a single protein from its pQTL summary
    statistics, an LD reference and (optionally) a functional annotation.
    """

    meta: Meta
    pqtl_source: SBayesRCSource
    ld_reference_directory_task: Task
    annotation_file_task: Task | None = None
    threads: int = 4

    @property
    def deps(self) -> list["Task"]:
        result: list[Task] = [
            self.pqtl_source.task,
            self.ld_reference_directory_task,
        ]
        if self.annotation_file_task is not None:
            result.append(self.annotation_file_task)
        return result

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        ld_asset = fetch(self.ld_reference_directory_task.asset_id)
        assert isinstance(ld_asset, DirectoryAsset)

        annot_args: list[str] = []
        if self.annotation_file_task is not None:
            annot_asset = fetch(self.annotation_file_task.asset_id)
            assert isinstance(annot_asset, FileAsset)
            annot_args = ["--annot", str(annot_asset.path)]

        pqtl_ma = prepare_cojo_ma_input_file(self.pqtl_source, fetch, scratch_dir)
        out_path = scratch_dir / f"weights{POLYPWAS_WEIGHTS_EXTENSION}"

        # polypwas train shells out to SBayesRC via Rscript; forward those calls
        # into the SBayesRC Docker image, mounting the LD reference so the
        # container can read it.
        shim_path = write_docker_rscript_shim(
            extra_mount_dirs=[ld_asset.path], shim_path=scratch_dir / "rscript_shim.sh"
        )
        polypwas_setup(shim_path, threads=self.threads)
        run_polypwas(
            [
                "train",
                "--pqtl",
                str(pqtl_ma),
                "--ldm-dir",
                str(ld_asset.path),
                "--threads",
                str(self.threads),
                *annot_args,
                "--out",
                str(out_path),
            ],
            threads=self.threads,
        )
        assert out_path.is_file(), f"polypwas weights not found at {out_path}"
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        pqtl_source: SBayesRCSource,
        ld_reference_directory_task: Task,
        annotation_file_task: Task | None = None,
        threads: int = 4,
    ) -> "PolypwasTrainTask":
        trait, project = derive_trait_project(pqtl_source.task.meta)
        meta = FilteredGWASDataMeta(
            id=AssetId(asset_id),
            trait=trait,
            project=project,
            sub_dir=PurePath("analysis") / "polypwas",
            extension=POLYPWAS_WEIGHTS_EXTENSION,
        )
        return cls(
            meta=meta,
            pqtl_source=pqtl_source,
            ld_reference_directory_task=ld_reference_directory_task,
            annotation_file_task=annotation_file_task,
            threads=threads,
        )
