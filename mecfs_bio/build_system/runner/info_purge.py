from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.verifying_trace_info import VerifyingTraceInfo
from mecfs_bio.build_system.tasks.simple_tasks import SimpleTasks


def info_purge(info: VerifyingTraceInfo,
               tasks: SimpleTasks
               ) ->VerifyingTraceInfo:
    """
    For each Task stored in Tasks, check whether its dependencies as expressed in info  are a subset of its declared task ependencies
    If this is not the case, purge it from info

    Returns updated info
    """