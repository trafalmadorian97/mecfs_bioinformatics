"""
Code configuring the default runner, which is used to run the main standard analysis

Some of the default runner options can bet overridden by adding a yaml file at
_DEFAULT_RUNNER_CONFIG_PATH.yaml

The config is machine-local and is not checked in.  Every key is independently optional,
so a config may set only the ones it wants to change:

    asset_root:  root of the asset store
    info_store:  path of the persistent verifying-trace cache
    path_remap:  rules routing selected store subtrees to other filesystems, so that one
                 logical store can span disks.  See the remapping_meta_to_path module for
                 the schema and for guidance on which subtrees are worth remapping.  Every
                 configured root must be present when a run starts, so a detached drive
                 aborts the run instead of silently rebuilding its assets elsewhere.
    remote_region:     AWS region used for remote compute and the S3 reference bucket.
    remote_s3_bucket:  bucket holding the staged GWFM reference / scratch data.
    gctb_image:        public GCTB container image reference.
"""

import functools
from pathlib import Path

import structlog
import yaml

from mecfs_bio.build_system.rebuilder.metadata_to_path.remapping_meta_to_path import (
    PathRemapRule,
)
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

CONFIG_FILE_NAME = "default_runner_config.yaml"
_DEFAULT_RUNNER_CONFIG_PATH = Path(CONFIG_FILE_NAME)

_ASSET_ROOT_KEY = "asset_root"
_INFO_STORE_KEY = "info_store"
_PATH_REMAP_KEY = "path_remap"
_REMOTE_REGION_KEY = "remote_region"
_REMOTE_BUCKET_KEY = "remote_s3_bucket"
_GCTB_IMAGE_KEY = "gctb_image"


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
    if config is None or _ASSET_ROOT_KEY not in config:
        return ASSET_ROOT
    return Path(config[_ASSET_ROOT_KEY])


def _get_info_store_path() -> Path:
    config = load_runner_config()
    if config is None or _INFO_STORE_KEY not in config:
        return IMO_128_INFO_STORE_PATH
    return Path(config[_INFO_STORE_KEY])


def get_path_remap_rules() -> tuple[PathRemapRule, ...]:
    config = load_runner_config()
    if config is None or _PATH_REMAP_KEY not in config:
        return ()
    return PathRemapRule.tuple_from_config(config[_PATH_REMAP_KEY])


def get_remote_region() -> str | None:
    config = load_runner_config()
    if config is None or _REMOTE_REGION_KEY not in config:
        return None
    return config[_REMOTE_REGION_KEY]


def get_remote_s3_bucket() -> str | None:
    config = load_runner_config()
    if config is None or _REMOTE_BUCKET_KEY not in config:
        return None
    return config[_REMOTE_BUCKET_KEY]


def get_gctb_image() -> str | None:
    config = load_runner_config()
    if config is None or _GCTB_IMAGE_KEY not in config:
        return None
    return config[_GCTB_IMAGE_KEY]


DEFAULT_RUNNER = SimpleRunner(
    tracer=_imo_hasher_128,  # _imo_hasher_32,#SimpleHasher.md5_hasher(),
    info_store=_get_info_store_path(),  # IMO_32_INFO_STORE_PATH,#MD5_INFO_STORE_PATH,
    asset_root=_get_asset_root_path(),
    path_remap=get_path_remap_rules(),
)
