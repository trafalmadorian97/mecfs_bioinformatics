"""Unit tests for the PreformattedSBayesRCDataSource valid-state invariant."""

import pytest

from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.sbayesrc.sbayesrc_data_source import (
    PreformattedSBayesRCDataSource,
)

_FILE_TASK = FakeTask(meta=SimpleFileMeta("ma_file"))
_DIR_TASK = FakeTask(meta=SimpleDirectoryMeta("ma_dir"))


def test_file_task_requires_no_filename():
    PreformattedSBayesRCDataSource(task=_FILE_TASK, filename=None, alias="x")
    with pytest.raises(AssertionError):
        PreformattedSBayesRCDataSource(task=_FILE_TASK, filename="trait.ma", alias="x")


def test_dir_task_requires_filename():
    PreformattedSBayesRCDataSource(task=_DIR_TASK, filename="trait.ma", alias="x")
    with pytest.raises(AssertionError):
        PreformattedSBayesRCDataSource(task=_DIR_TASK, filename=None, alias="x")
