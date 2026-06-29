"""
Data sources describing how a GWAS summary statistics asset is turned into the
COJO format .ma file consumed by SBayesRC and polypwas.

This mirrors the MiXeR data-source abstraction (see mixer_task.py): a
SBayesRCDataSource wraps a task that produces a gwaslab-format table and converts
it on the fly, while a PreformattedSBayesRCDataSource wraps a task whose asset is
already a COJO .ma file.
"""

import shutil
from pathlib import Path

from attrs import frozen

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import DirMeta, FileMeta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    PhenotypeInfo,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.sbayesrc.sbayesrc_cojo import write_cojo_ma


@frozen
class SBayesRCDataSource:
    """
    A source of GWAS summary statistics for SBayesRC / polypwas.

    The wrapped task must provide a dataframe in gwaslab format, which will be
    converted to COJO format (column renaming + N from the phenotype info).
    """

    task: Task
    alias: str
    sample_info: PhenotypeInfo
    pipe: DataProcessingPipe = IdentityPipe()

    @property
    def asset_id(self) -> AssetId:
        return self.task.asset_id


@frozen
class PreformattedSBayesRCDataSource:
    """
    A source of GWAS summary statistics already written in COJO .ma format.

    No gwaslab-to-COJO conversion is performed.  The wrapped task should provide a
    DirectoryAsset (in which case filename names the .ma file within it) or a
    FileAsset (in which case filename must be None, since the file is the asset).
    """

    task: Task
    filename: str | None
    alias: str

    def __attrs_post_init__(self) -> None:
        # Make invalid states unrepresentable: filename is meaningful only for a
        # directory-producing task, and must be absent for a file-producing one.
        meta = self.task.meta
        if isinstance(meta, DirMeta):
            assert self.filename is not None, (
                "filename is required when the wrapped task produces a directory asset"
            )
        elif isinstance(meta, FileMeta):
            assert self.filename is None, (
                "filename must be None when the wrapped task produces a file asset"
            )

    @property
    def asset_id(self) -> AssetId:
        return self.task.asset_id


SBayesRCSource = SBayesRCDataSource | PreformattedSBayesRCDataSource


def prepare_cojo_ma_input_file(
    source: SBayesRCSource,
    fetch: Fetch,
    temp_dir: Path,
) -> Path:
    """Materialize a COJO .ma file for source inside temp_dir and return its path."""
    if isinstance(source, PreformattedSBayesRCDataSource):
        source_asset = fetch(source.task.asset_id)
        if isinstance(source_asset, DirectoryAsset):
            assert source.filename is not None
            source_file = source_asset.path / source.filename
        elif isinstance(source_asset, FileAsset):
            source_file = source_asset.path
        else:
            raise ValueError(f"Unexpected asset type: {type(source_asset)}")
        assert source_file.is_file(), f"Source file not found: {source_file}"
        # Preserve the source name (and its suffix) so downstream tools detect the
        # format / compression correctly.
        dest = temp_dir / source_file.name
        shutil.copy(str(source_file), str(dest))
        return dest
    if isinstance(source, SBayesRCDataSource):
        asset = fetch(source.task.asset_id)
        frame = source.pipe.process(scan_dataframe_asset(asset, meta=source.task.meta))
        out_path = temp_dir / f"{source.alias}.ma"
        return write_cojo_ma(frame, source.sample_info, out_path)
    raise ValueError(f"Unexpected source type: {type(source)}")
