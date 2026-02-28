from copy import deepcopy

import attrs
import structlog

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.verifying_trace_info import (
    TraceRecord,
    VerifyingTraceInfo,
)
from mecfs_bio.build_system.tasks.simple_tasks import SimpleTasks

logger = structlog.get_logger()


def info_purge(info: VerifyingTraceInfo, tasks: SimpleTasks) -> VerifyingTraceInfo:
    """
    For each Task stored in Tasks, check whether its dependencies as expressed in info  are a subset of its declared task ependencies
    If this is not the case, purge it from info
    """
    info_dict = deepcopy(info.trace_store)

    for task in tasks.values():
        if task.asset_id in info_dict:
            recorded_deps = _get_deps_from_trace_record(info_dict[task.asset_id])
            actual_deps = set([item.asset_id for item in task.deps])
            if not (recorded_deps <= actual_deps):
                logger.debug(
                    f"Dependencies of asset {task.asset_id} appear to be outdated.\n Recorded deps:{recorded_deps}.  Actual deps: {actual_deps}.  Puring from info store."
                )
                info_dict.pop(task.asset_id)
    return attrs.evolve(info, trace_store=info_dict)


def _get_deps_from_trace_record(trace_record: TraceRecord) -> set[AssetId]:
    return set([item[0] for item in trace_record[1]])
