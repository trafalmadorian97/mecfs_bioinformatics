from pathlib import Path

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset


def read_file_asset_path(asset: Asset) -> Path:
    assert isinstance(asset, FileAsset)
    return asset.path
