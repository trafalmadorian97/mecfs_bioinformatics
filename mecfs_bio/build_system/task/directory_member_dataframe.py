"""
A dataframe stored as one named file inside a dependency's directory asset.

A DataFrameReadSpec describes a file format only, and a directory asset can hold several
dataframes, so neither the read spec nor the directory's meta says which file to read. A
consumer that reads one of them names it in its own settings with a DirectoryMemberDataFrame.
"""

from pathlib import PurePath

import narwhals as nw
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (
    ValidBackend,
    scan_dataframe,
)


@frozen(slots=True)
class DirectoryMemberDataFrame:
    """
    file_name: the file inside the dependency's directory asset (a single path component).
    read_spec: the file's format.
    """

    file_name: PurePath
    read_spec: DataFrameReadSpec

    def __attrs_post_init__(self) -> None:
        assert not self.file_name.is_absolute() and len(self.file_name.parts) == 1, (
            f"file_name must be one relative path component: {self.file_name}"
        )


def scan_directory_member(
    asset: Asset,
    member: DirectoryMemberDataFrame,
    parquet_backend: ValidBackend = "polars",
) -> nw.LazyFrame:
    """Scan the member file of a directory asset."""
    assert isinstance(asset, DirectoryAsset), (
        f"expected a DirectoryAsset holding {member.file_name}, got {type(asset).__name__}"
    )
    path = asset.path / member.file_name
    assert path.is_file(), f"{member.file_name} not found in {asset.path}"
    return scan_dataframe(path, member.read_spec, parquet_backend)
