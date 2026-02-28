from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.verifying_trace_info import VerifyingTraceInfo
from mecfs_bio.build_system.runner.info_purge import info_purge
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.tasks.simple_tasks import SimpleTasks


def test_info_purge():
    task_1 =FakeTask(
                SimpleFileMeta("asset_1"),
                tuple()
            )

    task_2 =FakeTask(
        SimpleFileMeta("asset_2"),
        (task_1,)
    )
    task_3=FakeTask(
        SimpleFileMeta("asset_3"),
        (task_1,task_2)
    )


    tasks = SimpleTasks(
        {
            AssetId("asset_1"): task_1,
            AssetId("asset_2"): task_2,
            AssetId("asset_3"): task_3
        }
    )
    info =VerifyingTraceInfo(
        trace_store={
            AssetId("asset_1"): ("abc",[]),
            AssetId("asset_2"):  ("def",[("asset_1","abc")])   ,
            AssetId("asset_3"): ("hij", [("asset_2","def")] )
        }
    )
    new_info=info_purge(
        info=info,
        tasks=tasks
    )
    assert set(new_info.trace_store.keys()) == {"asset_1", "asset_2"} #asset 3s dependencies are inconsistent between Tasks and info, so it is purged