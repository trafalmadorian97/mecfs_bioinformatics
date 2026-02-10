from pathlib import Path
from typing import Sequence

from attrs import frozen

from mecfs_bio.build_system.figure_manager.figure_exporter import FigureExporter
from mecfs_bio.build_system.task.base_task import Task


@frozen
class FigureManager:
     fig_dir: Path
     exporter: FigureExporter
     def generate_figs_from_tasks(self,
                                  tasks: Sequence[Task],
                                  )->None: