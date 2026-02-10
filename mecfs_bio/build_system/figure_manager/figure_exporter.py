from pathlib import Path
from typing import Callable, Mapping, Sequence

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.base_tracer import Tracer
from mecfs_bio.build_system.task.base_task import Task


@frozen
class FigureExporter:
    runner: Callable[[Sequence[Task]],Mapping[AssetId,Asset]]
    tracer: Tracer

    def export(self, to_export: Sequence[Task], fig_dir: Path) -> None:
        pass