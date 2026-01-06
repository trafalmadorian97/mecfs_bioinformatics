"""
Metadata describing an executable file
"""

from pathlib import PurePath

from attrs import field, frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta


@frozen
class ExecutableMeta(FileMeta):
    """
    Metadata for executable files
    """

    group: str
    sub_folder: PurePath
    extension: str | None
    filename: str | None
    _asset_id: AssetId = field(converter=AssetId)

    @property
    def asset_id(self) -> AssetId:
        return self._asset_id

    @classmethod
    def create(
        cls,
        group: str,
        sub_folder: PurePath,
        asset_id: str,
        filename: str | None = None,
        extension: str | None = None,
    ):
        return cls(
            group=group,
            sub_folder=sub_folder,
            filename=filename,
            extension=extension,
            asset_id=AssetId(asset_id),
        )
