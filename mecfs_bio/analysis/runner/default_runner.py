"""
Code configuring the default runner, which is used to run the main standard analysis

Some of the default runner options can bet overridden by adding a yaml file at
_DEFAULT_RUNNER_CONFIG_PATH.yaml
"""

import functools
from pathlib import Path

import structlog
import yaml

from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner

# MD5_INFO_STORE_PATH = Path("build_system")  /"verifying_trace_md5_info.yaml"
# IMO_32_INFO_STORE_PATH = Path("build_system") / "verifying_trace_imo_xxh_info.yaml"
IMO_128_INFO_STORE_PATH = Path("build_system") / "verifying_trace_imo_xxh_128_info.yaml"
ASSET_ROOT = Path("assets") / "base_asset_store"
# _imo_hasher_32 = ImoHasher.with_xxhash_32()
_imo_hasher_128 = ImoHasher.with_xxhash_128()
# _md5_hash = SimpleHasher.md5_hasher()

logger = structlog.get_logger(__name__)

_DEFAULT_RUNNER_CONFIG_PATH = Path("default_runner_config.yaml")

_ASSET_ROOT_KEY = "asset_root"
_INFO_STORE_KEY = "info_store"


@functools.cache
def load_runner_config() -> dict | None:
    if not _DEFAULT_RUNNER_CONFIG_PATH.exists():
        logger.debug(
            f"No default runner config file found at {_DEFAULT_RUNNER_CONFIG_PATH}"
        )
        return None
    with open(_DEFAULT_RUNNER_CONFIG_PATH) as infile:
        config = yaml.load(infile, Loader=yaml.FullLoader)
        logger.debug(
            f"Loading default runner config from {_DEFAULT_RUNNER_CONFIG_PATH}"
        )
        logger.debug(f"config: \n {config}")
        return config


def _get_asset_root_path() -> Path:
    config = load_runner_config()
    if config is None:
        return ASSET_ROOT
    return Path(config[_ASSET_ROOT_KEY])


def _get_info_store_path() -> Path:
    config = load_runner_config()
    if config is None:
        return IMO_128_INFO_STORE_PATH
    return Path(config[_INFO_STORE_KEY])


DEFAULT_RUNNER = SimpleRunner(
    tracer=_imo_hasher_128,  # _imo_hasher_32,#SimpleHasher.md5_hasher(),
    info_store=_get_info_store_path(),  # IMO_32_INFO_STORE_PATH,#MD5_INFO_STORE_PATH,
    asset_root=_get_asset_root_path(),
)
