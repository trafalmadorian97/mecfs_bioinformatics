import shutil
from contextlib import contextmanager
from pathlib import Path

from attrs import frozen

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class StampedExternalDirectoryTask(GeneratingTask):
    """
    Copy an external directory into the build system as a DirectoryAsset.  Used in
    system tests to wrap a downloaded LD reference directory.
    """

    meta: SimpleDirectoryMeta
    external_path: Path

    @property
    def deps(self) -> list[Task]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> DirectoryAsset:
        target = scratch_dir / "target"
        shutil.copytree(str(self.external_path), str(target))
        return DirectoryAsset(target)


@frozen
class StampedExternalFileTask(GeneratingTask):
    """
    Copy an external file into the build system, stamping it with an arbitrary
    meta.  Used in system tests to wrap downloaded example files as trait-bearing
    sources (e.g. a GWASSummaryDataFileMeta) for the SBayesRC / polypwas tasks.
    """

    meta: Meta
    external_path: Path

    @property
    def deps(self) -> list[Task]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> FileAsset:
        target = scratch_dir / "target"
        shutil.copy(str(self.external_path), str(target))
        return FileAsset(target)


@contextmanager
def log_on_error(filepath: Path):
    """

    Main use case:
    - A system test fails, and I want to spill the current state of the build system info store for debugging
    """
    try:
        print(f"Debug file path is {filepath}")
        yield
    except Exception:
        print(f"\n--- TRACE: Contents of {filepath} ---")
        try:
            with open(filepath) as f:
                print(f.read())
        except FileNotFoundError:
            print(f"Warning: Debug file '{filepath}' not found.")
        print("--- END OF TRACE ---\n")

        raise
