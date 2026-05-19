import gwaslab
import structlog

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset

logger = structlog.get_logger()


def read_sumstats(asset: Asset) -> gwaslab.Sumstats:
    assert isinstance(asset, FileAsset)
    logger.debug(f"reading sumstats from {asset.path}")
    result = gwaslab.load_pickle(str(asset.path))
    assert isinstance(result, gwaslab.Sumstats)
    return result
