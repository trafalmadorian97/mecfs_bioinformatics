import gwaslab

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
import structlog
logger = structlog.get_logger()

def read_sumstats(asset: Asset) -> gwaslab.Sumstats:
    assert isinstance(asset, FileAsset)
    logger.debug(f"reading sumstats from {asset.path}")
    return gwaslab.load_pickle(asset.path)
