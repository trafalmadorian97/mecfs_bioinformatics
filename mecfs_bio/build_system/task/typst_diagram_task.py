"""
Task that compiles Typst source into an SVG diagram asset.
"""

from collections.abc import Callable
from pathlib import Path

from attrs import field, frozen

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.diagram_file_meta import (
    DIAGRAM_EXTENSION,
    DiagramFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.task.typst_source import TypstSource
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.subproc.run_command import execute_command


@frozen
class TypstDiagramTask(GeneratingTask):
    """
    Compiles a Typst source (a .typ file or an inline string) into an SVG
    diagram by shelling out to the typst CLI.
    """

    meta: DiagramFileMeta
    source: TypstSource
    executor: Callable[[list[str]], str] = field(default=execute_command)

    @property
    def deps(self) -> list[Task]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> FileAsset:
        src_path = self.source.resolve_to_path(scratch_dir)
        out_path = scratch_dir / (self.meta.asset_id + DIAGRAM_EXTENSION)
        self.executor(["typst", "compile", str(src_path), str(out_path)])
        return FileAsset(out_path)
