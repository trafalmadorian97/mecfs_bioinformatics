from typing import Iterator, Mapping

from attrs import frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.tasks.base_tasks import Tasks


@frozen
class SimpleTasks(Tasks):
    _tasks: Mapping[AssetId, Task]

    def __attrs_post_init__(self):
        for asset_id, task in self._tasks.items():
            assert task.asset_id == asset_id

    def __getitem__(self, asset_id: AssetId) -> Task:
        return self._tasks[asset_id]

    def __len__(self) -> int:
        return len(self._tasks)

    def __iter__(self) -> Iterator[AssetId]:
        yield from self._tasks


def find_tasks(target_tasks: list[Task]) -> SimpleTasks:
    """
    Build a SimpleTasks object by walking the task graph
    """
    _tasks: dict[AssetId, Task] = {}
    visited = set()

    def explore_task(t: Task):
        visited.add(t.asset_id)
        for dep in t.deps:
            if dep.asset_id not in visited:
                explore_task(dep)
        if t.asset_id in _tasks and (t != _tasks[t.asset_id]):
            raise ValueError(
                f"Found two tasks with asset id  {t.asset_id}.  tasks are {t} and {_tasks[t.asset_id]}"
            )
        _tasks[t.asset_id] = t

    for task in target_tasks:
        explore_task(task)
    return SimpleTasks(_tasks)
