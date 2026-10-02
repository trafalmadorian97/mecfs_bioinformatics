import typeguard
from typeguard import TypeCheckError

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.base_meta import DirMeta, FileMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.metadata_to_path.base_meta_to_path import (
    MetaToPath,
)


def get_asset_if_exists(
    meta: Meta,
    meta_to_path: MetaToPath,
) -> Asset | None:
    # asset_type = meta.asset_type
    # if asset_type == FileAsset:
    if _isinstance(meta, FileMeta):
        path = meta_to_path(meta)
        if path.is_file():
            return FileAsset(path)
        return None
    if _isinstance(meta, DirMeta):
        path = meta_to_path(meta)
        if path.is_dir():
            return DirectoryAsset(path)
        return None
    raise ValueError("Unknown asset type")


def _isinstance(val, expected_type):
    try:
        typeguard.check_type(val, expected_type)
    except TypeCheckError:
        return False
    return True
