from pathlib import Path

import pytest

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.metadata_to_path.simple_meta_to_path import (
    SimpleMetaToPath,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.base_tracer import (
    Tracer,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
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
    dependencies_of_targets_from_tasks,
    topological,
)
from mecfs_bio.build_system.task.copy_task import CopyTask
from mecfs_bio.build_system.task.counting_task import CountingTask
from mecfs_bio.build_system.task.discard_deps_task_wrapper import (
    DiscardDepsWrapper,
)
from mecfs_bio.build_system.task.external_file_copy_task import ExternalFileCopyTask
from mecfs_bio.build_system.task.failing_task import FailingTask
from mecfs_bio.build_system.tasks.simple_tasks import find_tasks
from mecfs_bio.build_system.wf.base_wf import SimpleWF

#


@pytest.mark.parametrize(
    argnames="tracer",
    argvalues=[
        SimpleHasher.md5_hasher(),
        ImoHasher.with_xxhash_32(),
        ImoHasher.with_xxhash_128(),
    ],
)
def test_file_copying_task(tmp_path: Path, tracer: Tracer) -> None:
    """
    Test a number of basic properties of the topological scheduler
    """
    incremental_save_path = tmp_path / "incremental_save"
    external_dir = tmp_path / "external"
    external_dir.mkdir(exist_ok=True, parents=True)
    external_file = external_dir / "external_file.txt"
    external_file.write_text("abc123")
    task1 = CountingTask(
        ExternalFileCopyTask(
            meta=SimpleFileMeta(AssetId("file_1")), external_path=external_file
        )
    )

    task2 = CountingTask(
        CopyTask(
            meta=SimpleFileMeta(
                AssetId("file_2"),
            ),
            dep_file_task=task1,
        )
    )

    task3 = CountingTask(
        CopyTask(
            meta=SimpleFileMeta(
                AssetId("file_3"),
            ),
            dep_file_task=task2,
        )
    )

    tasks = find_tasks([task3])

    wf = SimpleWF()
    info: VerifyingTraceInfo = VerifyingTraceInfo.empty()

    asset_dir = tmp_path / "asset_dir"
    asset_dir.mkdir(exist_ok=True, parents=True)
    meta_to_path = SimpleMetaToPath(root=asset_dir)

    rebuilder = VerifyingTraceRebuilder(tracer)

    targets = [task3.meta.asset_id]

    # Verify that all files are created in the correct location
    store, info = topological(
        rebuilder=rebuilder,
        tasks=tasks,
        targets=targets,
        wf=wf,
        info=info,
        meta_to_path=meta_to_path,
        incremental_save_path=incremental_save_path,
    )

    file_1_path = meta_to_path(task1.meta)
    file_2_path = meta_to_path(task2.meta)
    file_3_path = meta_to_path(task3.meta)

    expected_store = {
        task1.meta.asset_id: FileAsset(file_1_path),
        task2.meta.asset_id: FileAsset(file_2_path),
        task3.meta.asset_id: FileAsset(file_3_path),
    }
    assert expected_store == store
    assert task1.run_count == 1
    assert task2.run_count == 1
    assert task3.run_count == 1

    # verify that only the 3rd task if rerun if the 3rd file is deleted

    file_3_path.unlink()
    topological(
        rebuilder=rebuilder,
        tasks=tasks,
        targets=targets,
        wf=wf,
        info=info,
        meta_to_path=meta_to_path,
        incremental_save_path=incremental_save_path,
    )

    assert task1.run_count == 1
    assert task2.run_count == 1
    assert task3.run_count == 2

    # verify early_cutoff: deleting the second file causes only the second task to be rerun

    file_2_path.unlink()
    topological(
        rebuilder=rebuilder,
        tasks=tasks,
        targets=targets,
        wf=wf,
        info=info,
        meta_to_path=meta_to_path,
        incremental_save_path=incremental_save_path,
    )

    assert task1.run_count == 1
    assert task2.run_count == 2
    assert task3.run_count == 2

    # check that verification failure works: if the input data changes, all downstream tasks are run

    file_1_path.unlink()
    external_file.write_text("Modified file")
    topological(
        rebuilder=rebuilder,
        tasks=tasks,
        targets=targets,
        wf=wf,
        info=info,
        meta_to_path=meta_to_path,
        incremental_save_path=incremental_save_path,
    )

    assert task1.run_count == 2
    assert task2.run_count == 3
    assert task3.run_count == 3

    # check that a deserialized info object can be used with the scheduler
    serialization_loc = tmp_path / "serialization_loc" / "info.yaml"
    info.serialize(serialization_loc)
    info_2 = VerifyingTraceInfo.deserialize(serialization_loc)
    topological(
        rebuilder=rebuilder,
        tasks=tasks,
        targets=targets,
        wf=wf,
        info=info_2,
        meta_to_path=meta_to_path,
        incremental_save_path=incremental_save_path,
    )

    assert task1.run_count == 2
    assert task2.run_count == 3
    assert task3.run_count == 3

    # check that we can use DiscardDepsWrapper to materialize the dependencies of
    # task 3 in a temporary directory

    wrapped_task_3 = DiscardDepsWrapper(
        inner=task3,
    )
    assert len(wrapped_task_3.deps) == 0
    tasks_wrapped = find_tasks([wrapped_task_3])

    wrapped_targets = [wrapped_task_3.asset_id]
    topological(
        rebuilder=rebuilder,
        tasks=tasks_wrapped,
        targets=wrapped_targets,
        wf=wf,
        info=info_2,
        meta_to_path=meta_to_path,
        incremental_save_path=incremental_save_path,
    )

    assert task1.run_count == 3
    assert task2.run_count == 4
    assert task3.run_count == 4

    topological(
        rebuilder=rebuilder,
        tasks=tasks_wrapped,
        targets=wrapped_targets,
        wf=wf,
        info=info_2,
        meta_to_path=meta_to_path,
        incremental_save_path=incremental_save_path,
    )

    assert task1.run_count == 3
    assert task2.run_count == 4
    assert task3.run_count == 4


def test_graph_generation(tmp_path: Path):
    """
    verify we can correctly calculate the dependencies and dependees of nodes in a simple
    graph
    """
    external_dir = tmp_path / "external"
    external_dir.mkdir(exist_ok=True, parents=True)
    external_file = external_dir / "external_file.txt"
    external_file.write_text("abc123")
    task1 = ExternalFileCopyTask(
        meta=SimpleFileMeta(AssetId("file_1")), external_path=external_file
    )

    task2 = CopyTask(
        meta=SimpleFileMeta(
            AssetId("file_2"),
        ),
        dep_file_task=task1,
    )

    task3 = CopyTask(
        meta=SimpleFileMeta(
            AssetId("file_3"),
        ),
        dep_file_task=task2,
    )
    tasks = find_tasks([task3])
    graph_1 = dependencies_of_targets_from_tasks(tasks, [task1.asset_id])
    graph_2 = dependencies_of_targets_from_tasks(tasks, [task2.asset_id])
    graph_3 = dependencies_of_targets_from_tasks(tasks, [task3.asset_id])
    assert len(graph_1) == 1
    assert len(graph_2) == 2
    assert len(graph_3) == 3
    graph_a = dependees_of_targets_from_tasks(tasks, [task1.asset_id])
    graph_b = dependees_of_targets_from_tasks(tasks, [task2.asset_id])
    graph_c = dependees_of_targets_from_tasks(tasks, [task3.asset_id])
    assert len(graph_a) == 3
    assert len(graph_b) == 2
    assert len(graph_c) == 1


def test_incremental_save_with_failing_task(tmp_path: Path) -> None:
    """
    If the scheduler fails due to an error, the info store should still be saved
    to the incremental save path
    """
    tracer = SimpleHasher.md5_hasher()
    incremental_save_path = tmp_path / "incremental_save"
    external_dir = tmp_path / "external"
    external_dir.mkdir(exist_ok=True, parents=True)
    external_file = external_dir / "external_file.txt"
    external_file.write_text("abc123")
    task1 = CountingTask(
        ExternalFileCopyTask(
            meta=SimpleFileMeta(AssetId("file_1")), external_path=external_file
        )
    )

    task2 = FailingTask(meta=SimpleFileMeta(AssetId("file_2")), deps=[task1])

    tasks = find_tasks([task2])

    wf = SimpleWF()
    info: VerifyingTraceInfo = VerifyingTraceInfo.empty()

    asset_dir = tmp_path / "asset_dir"
    asset_dir.mkdir(exist_ok=True, parents=True)
    meta_to_path = SimpleMetaToPath(root=asset_dir)

    rebuilder = VerifyingTraceRebuilder(tracer)

    targets = [task2.meta.asset_id]

    assert not incremental_save_path.exists()
    with pytest.raises(ValueError):
        topological(
            rebuilder=rebuilder,
            tasks=tasks,
            targets=targets,
            wf=wf,
            info=info,
            meta_to_path=meta_to_path,
            incremental_save_path=incremental_save_path,
        )
    assert incremental_save_path.exists()
    assert task1.run_count == 1
    info = VerifyingTraceInfo.deserialize(incremental_save_path)
    with pytest.raises(ValueError):
        topological(
            rebuilder=rebuilder,
            tasks=tasks,
            targets=targets,
            wf=wf,
            info=info,
            meta_to_path=meta_to_path,
            incremental_save_path=incremental_save_path,
        )
    assert incremental_save_path.exists()
    assert (
        task1.run_count == 1
    )  # we do not need to rerun task 1, even though task 2 crashed
