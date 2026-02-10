from pathlib import Path
from typing import Mapping, Sequence

import structlog
from rich.pretty import pprint

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset

logger = structlog.get_logger()
import attrs
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.rebuilder.metadata_to_path.simple_meta_to_path import (
    SimpleMetaToPath,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.base_tracer import (
    Tracer,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.simple_hasher import (
    SimpleHasher,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.verifying_trace_info import (
    VerifyingTraceInfo,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.verifying_trace_rebuilder_core import (
    VerifyingTraceRebuilder,
)
from mecfs_bio.build_system.scheduler.topological_scheduler import (
    dependees_of_targets_from_tasks,
    topological,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.tasks.simple_tasks import find_tasks
from mecfs_bio.build_system.wf.base_wf import RobustDownloadWF


@frozen
class SimpleRunner:
    """
    Simple wrapper class that orchestrates
    an execution of the workflow build system with a topological scheduler and verifying trace rebuilder

    info_store: path at which to store persistent cash for build-system internal information
    asset_root: root under which to create the asset store
    tracer: algorithm uses to calculate verifying traces of assets.  An example would be a hashing algorithm.  changing this forces all assets to be rebuilt
    """

    info_store: Path
    asset_root: Path
    tracer: Tracer = SimpleHasher.md5_hasher()

    @property
    def meta_to_path(self) -> SimpleMetaToPath:
        return SimpleMetaToPath(root=self.asset_root)

    def run(
        self,
        targets: list[Task],
        must_rebuild_transitive: Sequence[Task] = tuple(),
        incremental_save: bool = False,
    ) -> Mapping[AssetId, Asset]:
        """
        Targets: the ultimate targets we aim to produce.  All transitive dependencies of these targets will either be rebuilt, or fetched (determined according to that status of their trace)
        must_rebuild_transitive: list of tasks that the rebuilder will be forced to rebuild, regardless of the status of their trace.
           - This is particularly useful when you have changed the code that generates and asset, and so want it and its depndees to be regenerated.
        returns:
        mapping from asset id to file system information for all assets that were built or retrieved as part of the execution of the scheduler
        """
        if self.info_store.is_file():
            info = VerifyingTraceInfo.deserialize(self.info_store)
        else:
            info = VerifyingTraceInfo.empty()
        if incremental_save:
            logger.debug("incremental save is enabled")
        rebuilder = VerifyingTraceRebuilder(self.tracer)
        wf = RobustDownloadWF()
        # wf= SimpleWF()
        meta_to_path = self.meta_to_path
        tasks = find_tasks(targets)
        must_rebuild_graph = dependees_of_targets_from_tasks(
            tasks=tasks,
            targets=[task.asset_id for task in must_rebuild_transitive],
        )
        info = attrs.evolve(info, must_rebuild=set(must_rebuild_graph.nodes))
        incremental_save_path = self.info_store if incremental_save else None
        store, info = topological(
            rebuilder=rebuilder,
            tasks=tasks,
            info=info,
            wf=wf,
            meta_to_path=meta_to_path,
            targets=[target.asset_id for target in targets],
            incremental_save_path=incremental_save_path,
        )
        info.serialize(self.info_store)
        _print_target_locs(store, targets=targets)
        return store


def _print_target_locs(store: Mapping[AssetId, Asset], targets: list[Task]):
    logger.info("Locations of materialized targets:\n")
    target_locs: dict = {}
    for item in targets:
        asset = store[item.asset_id]
        if isinstance(asset, (FileAsset, DirectoryAsset)):
            target_locs[item.asset_id] = str(asset.path)
        else:
            target_locs[item.asset_id] = asset
    pprint(target_locs)
