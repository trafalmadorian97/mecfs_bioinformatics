"""
Task that computes polypwas PWAS Z-scores from trained protein weights and GWAS
summary statistics.

polypwas assoc turns weights plus a GWAS into cis/trans PWAS Z-scores; it needs no
R.  Output is a TSV with columns ID, CIS_Z, TRANS_Z.
"""

from pathlib import Path, PurePath

from attrs import frozen
from structlog import get_logger

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.base_meta import DirMeta, FileMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.polypwas.polypwas_common import (
    derive_trait_project,
    run_polypwas,
)
from mecfs_bio.build_system.task.sbayesrc.sbayesrc_data_source import (
    SBayesRCSource,
    prepare_cojo_ma_input_file,
)
from mecfs_bio.build_system.wf.base_wf import WF

logger = get_logger()

POLYPWAS_ASSOC_OUTPUT_FILENAME = "pwas_z.tsv"


def _fetch_file_path(task: Task, fetch: Fetch, filename: str | None = None) -> Path:
    """Resolve a task's asset to a file path, optionally inside a directory asset."""
    asset = fetch(task.asset_id)
    if isinstance(asset, FileAsset):
        return asset.path
    if isinstance(asset, DirectoryAsset):
        assert filename is not None, "filename required for a directory asset"
        path = asset.path / filename
        assert path.is_file(), f"{path} not found"
        return path
    raise ValueError(f"Unexpected asset type: {type(asset)}")


@frozen
class PolypwasAssocTask(Task):
    """
    Compute PWAS cis/trans Z-scores from polypwas weights and GWAS summary stats.
    """

    meta: Meta
    weights_task: Task
    gwas_source: SBayesRCSource
    ld_reference_directory_task: Task
    gene_info_task: Task
    threads: int = 4
    weights_filename: str | None = None

    def __attrs_post_init__(self) -> None:
        # weights_filename selects the weights file only when the weights task
        # produces a directory; it must be absent when the task is the file itself.
        weights_meta = self.weights_task.meta
        if isinstance(weights_meta, DirMeta):
            assert self.weights_filename is not None, (
                "weights_filename is required when the weights task produces a "
                "directory asset"
            )
        elif isinstance(weights_meta, FileMeta):
            assert self.weights_filename is None, (
                "weights_filename must be None when the weights task produces a "
                "file asset"
            )

    @property
    def deps(self) -> list["Task"]:
        return [
            self.weights_task,
            self.gwas_source.task,
            self.ld_reference_directory_task,
            self.gene_info_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        ld_asset = fetch(self.ld_reference_directory_task.asset_id)
        assert isinstance(ld_asset, DirectoryAsset)

        weights_path = _fetch_file_path(self.weights_task, fetch, self.weights_filename)
        gene_info_path = _fetch_file_path(self.gene_info_task, fetch)
        gwas_ma = prepare_cojo_ma_input_file(self.gwas_source, fetch, scratch_dir)
        out_path = scratch_dir / POLYPWAS_ASSOC_OUTPUT_FILENAME

        run_polypwas(
            [
                "assoc",
                "--weights",
                str(weights_path),
                "--gwas",
                str(gwas_ma),
                "--ldm-dir",
                str(ld_asset.path),
                "--gene-info",
                str(gene_info_path),
                "--out",
                str(out_path),
            ],
            threads=self.threads,
        )
        assert out_path.is_file(), f"polypwas assoc output not found at {out_path}"
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        weights_task: Task,
        gwas_source: SBayesRCSource,
        ld_reference_directory_task: Task,
        gene_info_task: Task,
        threads: int = 4,
        weights_filename: str | None = None,
    ) -> "PolypwasAssocTask":
        trait, project = derive_trait_project(gwas_source.task.meta)
        meta = ResultTableMeta(
            id=asset_id,
            trait=trait,
            project=project,
            extension=".tsv",
            sub_dir=PurePath("analysis") / "polypwas",
            read_spec=DataFrameReadSpec(DataFrameTextFormat(separator="\t")),
        )
        return cls(
            meta=meta,
            weights_task=weights_task,
            gwas_source=gwas_source,
            ld_reference_directory_task=ld_reference_directory_task,
            gene_info_task=gene_info_task,
            threads=threads,
            weights_filename=weights_filename,
        )
