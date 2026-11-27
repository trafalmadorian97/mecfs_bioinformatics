"""
A local copy of one deCODE SomaScan per-protein sumstats file, as published (gzipped text,
about 970 MB), checked against the manifest's size and md5. The database tasks discard their
downloads; this keeps one file for exploratory checks and brute-force comparisons.
"""

from pathlib import Path, PurePath

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.task.decode_database.construct_decode_variant_index_task import (
    download_protein_file,
)
from mecfs_bio.build_system.task.decode_database.decode_somascan_protein import (
    DecodeSomascanProtein,
)
from mecfs_bio.build_system.wf.base_wf import WF

_GZ_SUFFIX = ".txt.gz"


@frozen(slots=True)
class DownloadDecodeSomascanProteinTask(GeneratingTask):
    """
    protein: the manifest row of the file.
    n_connections: parallel range requests for the download.
    """

    meta: ReferenceFileMeta
    protein: DecodeSomascanProtein
    n_connections: int = 128

    @property
    def deps(self) -> list[Task]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        out_path = scratch_dir / PurePath(self.protein.key).name
        download_protein_file(self.protein, out_path, wf, self.n_connections)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls, asset_id: str, protein: DecodeSomascanProtein
    ) -> "DownloadDecodeSomascanProteinTask":
        name = PurePath(protein.key).name
        assert name.endswith(_GZ_SUFFIX), f"unexpected deCODE file name {name}"
        return cls(
            meta=ReferenceFileMeta(
                group="decode_somascan",
                sub_group="raw_protein_files",
                sub_folder=PurePath(protein.key).parent,
                id=AssetId(asset_id),
                filename=name.removesuffix(_GZ_SUFFIX),
                extension=_GZ_SUFFIX,
                read_spec=DataFrameReadSpec(
                    DataFrameTextFormat(separator="\t", null_values=["."])
                ),
            ),
            protein=protein,
        )
