"""
An Asset representing a file.
"""

from pathlib import Path

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset


@frozen
class FileAsset(Asset):
    """
    A materialized file output of a Task.
    """

    path: Path

    def __attrs_post_init__(self):
        assert self.path.is_file()
