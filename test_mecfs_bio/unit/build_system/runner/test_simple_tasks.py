from pathlib import Path

import pytest

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.external_file_copy_task import ExternalFileCopyTask
from mecfs_bio.build_system.tasks.simple_tasks import find_tasks


def test_find_tasks_duplicate_ids(tmp_path: Path):
    """
    find_tasks should reject duplicate IDs
    """
    file_1_path = tmp_path / "file_1.csv"
    file_2_path = tmp_path / "file_2.csv"
    file_1_path.write_text("1,2,3,4,5,6,7,8,9,10")
    file_2_path.write_text("1,2,3,4,5,6,7,8,9,10")
    task_1 = ExternalFileCopyTask(
        meta=SimpleFileMeta(AssetId("my_file")), external_path=file_1_path
    )
    task_2 = ExternalFileCopyTask(
        meta=SimpleFileMeta(AssetId("my_file")), external_path=file_2_path
    )
    with pytest.raises(ValueError):
        find_tasks([task_1, task_2])
