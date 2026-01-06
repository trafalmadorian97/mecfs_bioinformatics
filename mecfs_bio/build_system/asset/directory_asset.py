"""
An Asset representing a directory,
"""

from pathlib import Path

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset


@frozen
class DirectoryAsset(Asset):
    """
    A materialized directory output of a Task.
    """

    path: Path

    def __attrs_post_init__(self):
        assert self.path.is_dir()
