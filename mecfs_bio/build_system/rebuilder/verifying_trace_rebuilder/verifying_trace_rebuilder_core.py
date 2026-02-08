from pathlib import Path
from typing import Callable

import structlog
from attrs import frozen
from loguru import logger

logger = structlog.get_logger()

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.rebuilder.base_rebuilder import Rebuilder
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.rebuilder.fetch.restricted_fetch import RestrictedFetch
from mecfs_bio.build_system.rebuilder.metadata_to_path.base_meta_to_path import (
    MetaToPath,
)
from mecfs_bio.build_system.rebuilder.tracking_sandboxed_execute import (
    tracking_sandboxed_execute,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.base_tracer import (
    Tracer,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.verifying_trace_info import (
    VerifyingTraceInfo,
    update_verifying_trace_info_in_place,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class VerifyingTraceRebuilder(Rebuilder[VerifyingTraceInfo]):
    """
    A rebuilder that calculates traces for the assets it manages, and uses these traces
    to decide when to rebuild.
    Based on:
    Mokhov, Andrey, Neil Mitchell, and Simon Peyton Jones.
    "Build systems à la carte: Theory and practice."
    Journal of Functional Programming 30 (2020): e11.
    """

    tracer: Tracer

    def rebuild(
        self,
        task: Task,
        asset: Asset | None,
        fetch: Fetch,
        wf: WF,
        info: VerifyingTraceInfo,
        meta_to_path: MetaToPath,
    ) -> tuple[Asset, VerifyingTraceInfo]:
        must_rebuild = task.asset_id in info.must_rebuild

        def fetch_trace(asset_id: AssetId) -> str:
            return self.tracer(fetch(asset_id))

        if not must_rebuild and asset is not None:
            logger.debug(f"Attempting to verify the trace of asset {task.asset_id}...")
            old_value_trace = self.tracer(asset)
            if verify_trace(
                asset_id=task.asset_id,
                value_trace=old_value_trace,
                fetch_trace=fetch_trace,
                info=info,
            ):
                logger.debug(
                    f"Successfully verified the trace of asset {task.asset_id}."
                )
                return asset, info
            logger.debug(f"Failed to verify the trace of asset {task.asset_id}.")
        logger.debug(f"Materializing asset {task.asset_id}....")
        new_value, deps = tracking_sandboxed_execute(
            task=task,
            meta_to_path=meta_to_path,
            wf=wf,
            fetch=RestrictedFetch.from_task(fetch=fetch, task=task),
        )
        deps_traced = [(k, self.tracer(v)) for k, v in deps]
        new_trace = self.tracer(new_value)
        logger.debug(f"Trace of asset {task.asset_id} is:\n {new_trace}")
        update_verifying_trace_info_in_place(
            verifying_trace_info=info,
            meta=task.meta,
            new_value_trace=new_trace,
            deps_traced=deps_traced,
        )
        return new_value, info

    @classmethod
    def save_info(cls, info: VerifyingTraceInfo, path: Path):
        info.serialize(path)


def verify_trace(
    asset_id: AssetId,
    value_trace: str,
    fetch_trace: Callable[[AssetId], str],
    info: VerifyingTraceInfo,
) -> bool:
    if asset_id not in info.trace_store:
        return False
    recorded_trace, recorded_deps = info.trace_store[asset_id]
    if recorded_trace != value_trace:
        return False
    for dep, dep_hash in recorded_deps:
        if fetch_trace(dep) != dep_hash:
            return False
    return True
