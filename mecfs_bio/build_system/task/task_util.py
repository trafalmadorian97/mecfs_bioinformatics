from mecfs_bio.build_system.meta.base_meta import FileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec
from mecfs_bio.build_system.task.base_task import Task


def produces_dataframe(task: Task) -> bool:
    """Whether a task's asset can be read as a dataframe, i.e. whether its meta declares a
    DataFrameReadSpec.
    """

    if isinstance(task.meta, FileMeta):
        return isinstance(task.meta.read_spec, DataFrameReadSpec)
    return False
