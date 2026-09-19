"""
Decompress a gzipped genome FASTA and index it with faidx.

Memory-mapped reference lookups need the uncompressed file plus its .fai, so the
output is a directory holding both. Wrap instances in DiscardDepsWrapper so that the
gzipped download is not also kept in the asset store.
"""

import gzip
import shutil
from pathlib import Path, PurePath

import pysam
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.fasta_meta import FASTAMeta
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    FASTA_FILENAME,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild

_COPY_CHUNK_BYTES = 16 * 1024 * 1024


@frozen
class IndexedFastaTask(Task):
    meta: FASTAMeta
    fasta_gz_task: Task

    @property
    def deps(self) -> list[Task]:
        return [self.fasta_gz_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source = fetch(self.fasta_gz_task.asset_id)
        assert isinstance(source, FileAsset), (
            f"expected {self.fasta_gz_task.asset_id} to be a FileAsset, got {type(source).__name__}"
        )
        out_dir = scratch_dir / "indexed_fasta"
        out_dir.mkdir()
        fasta_path = out_dir / FASTA_FILENAME
        with (
            gzip.open(source.path, "rb") as compressed,
            open(fasta_path, "wb") as plain,
        ):
            shutil.copyfileobj(compressed, plain, length=_COPY_CHUNK_BYTES)
        pysam.faidx(str(fasta_path))
        return DirectoryAsset(out_dir)

    @classmethod
    def create(
        cls, fasta_gz_task: Task, asset_id: str, build: GenomeBuild
    ) -> "IndexedFastaTask":
        source_meta = fasta_gz_task.meta
        assert isinstance(source_meta, ReferenceFileMeta), (
            f"expected a ReferenceFileMeta source for {asset_id}, got {type(source_meta).__name__}"
        )
        return cls(
            meta=FASTAMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                build=build,
            ),
            fasta_gz_task=fasta_gz_task,
        )
