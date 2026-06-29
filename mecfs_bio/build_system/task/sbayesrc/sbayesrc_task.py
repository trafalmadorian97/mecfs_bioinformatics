"""
Task that runs SBayesRC to produce polygenic SNP weights from GWAS summary
statistics.

SBayesRC (Zheng et al., Nature Genetics 2024, doi:10.1038/s41588-024-01704-y)
incorporates functional genomic annotations with high-density SNPs for polygenic
prediction, requiring only GWAS summary statistics in COJO format plus an
eigen-decomposed LD reference.  This task runs the standard impute then sbayes RC
sequence through the SBayesRC Docker image; see sbayesrc_utils for why we use the
image rather than a local R install.
"""

import os
import shutil
import tempfile
from pathlib import Path, PurePath

from attrs import frozen
from structlog import get_logger

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.sbayesrc.sbayesrc_data_source import (
    SBayesRCSource,
    prepare_cojo_ma_input_file,
)
from mecfs_bio.build_system.task.sbayesrc.sbayesrc_utils import invoke_sbayesrc
from mecfs_bio.build_system.wf.base_wf import WF

logger = get_logger()

SBAYESRC_OUTPUT_PREFIX = "sbrc"
SBAYESRC_WEIGHTS_FILENAME = f"{SBAYESRC_OUTPUT_PREFIX}.txt"

CONTAINER_LD_DIR = Path("/ld_ref")
CONTAINER_ANNOT_DIR = Path("/annot")


def _derive_result_directory_meta(
    asset_id: str, source_meta: Meta, sub_dir: PurePath
) -> Meta:
    """Derive the output directory meta, inheriting trait/project where available."""
    if isinstance(
        source_meta,
        (FilteredGWASDataMeta, GWASLabSumStatsMeta, GWASSummaryDataFileMeta),
    ):
        return ResultDirectoryMeta(
            id=asset_id,
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=sub_dir,
        )
    if isinstance(source_meta, SimpleDirectoryMeta):
        return SimpleDirectoryMeta(id=AssetId(asset_id))
    raise ValueError(f"Cannot derive result meta from source meta {source_meta}")


@frozen
class SBayesRCTask(Task):
    """
    Run SBayesRC (impute -> sbayes RC) on a GWAS summary statistics source.

    The output directory holds the SNP weights (sbrc.txt) plus the SBayesRC log
    and parameter sidecar files.
    """

    meta: Meta
    gwas_source: SBayesRCSource
    ld_reference_directory_task: Task
    annotation_file_task: Task | None = None
    threads: int = 4
    impute: bool = True

    @property
    def deps(self) -> list["Task"]:
        result: list[Task] = [
            self.gwas_source.task,
            self.ld_reference_directory_task,
        ]
        if self.annotation_file_task is not None:
            result.append(self.annotation_file_task)
        return result

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        ld_asset = fetch(self.ld_reference_directory_task.asset_id)
        assert isinstance(ld_asset, DirectoryAsset)
        extra_mounts: dict[Path, Path] = {ld_asset.path.resolve(): CONTAINER_LD_DIR}

        annot_args: list[str] = []
        if self.annotation_file_task is not None:
            annot_asset = fetch(self.annotation_file_task.asset_id)
            assert isinstance(annot_asset, FileAsset)
            extra_mounts[annot_asset.path.parent.resolve()] = CONTAINER_ANNOT_DIR
            annot_args = ["--annot", str(CONTAINER_ANNOT_DIR / annot_asset.path.name)]

        common_args = [
            "--ldm-eigen",
            str(CONTAINER_LD_DIR),
            "--threads",
            str(self.threads),
        ]

        # The Docker working directory is the host cwd (mounted at /home), so we
        # stage inputs and outputs in a temp dir under cwd and refer to them by
        # relative path, mirroring the MiXeR task.
        with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            tmp_path = Path(tmpdir).relative_to(os.getcwd())
            ma_path = prepare_cojo_ma_input_file(self.gwas_source, fetch, tmp_path)
            out_prefix = tmp_path / SBAYESRC_OUTPUT_PREFIX

            if self.impute:
                invoke_sbayesrc(
                    common_args
                    + [
                        "--gwas-summary",
                        str(ma_path),
                        "--impute-summary",
                        "--out",
                        str(out_prefix),
                    ],
                    extra_mounts=extra_mounts,
                )
                sbayes_input = f"{out_prefix}.imputed.ma"
            else:
                sbayes_input = str(ma_path)

            invoke_sbayesrc(
                common_args
                + [
                    "--gwas-summary",
                    sbayes_input,
                    "--sbayes",
                    "RC",
                    "--out",
                    str(out_prefix),
                ]
                + annot_args,
                extra_mounts=extra_mounts,
            )

            weights_path = tmp_path / SBAYESRC_WEIGHTS_FILENAME
            assert weights_path.is_file(), (
                f"SBayesRC weights not found at {weights_path}"
            )
            for produced in tmp_path.iterdir():
                if produced.name.startswith(SBAYESRC_OUTPUT_PREFIX):
                    shutil.move(str(produced), str(scratch_dir / produced.name))

        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        gwas_source: SBayesRCSource,
        ld_reference_directory_task: Task,
        annotation_file_task: Task | None = None,
        threads: int = 4,
        impute: bool = True,
    ) -> "SBayesRCTask":
        meta = _derive_result_directory_meta(
            asset_id=asset_id,
            source_meta=gwas_source.task.meta,
            sub_dir=PurePath("analysis") / "sbayesrc",
        )
        return cls(
            meta=meta,
            gwas_source=gwas_source,
            ld_reference_directory_task=ld_reference_directory_task,
            annotation_file_task=annotation_file_task,
            threads=threads,
            impute=impute,
        )
