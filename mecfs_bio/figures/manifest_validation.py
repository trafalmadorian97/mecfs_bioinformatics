"""
Check that the figure manifest is a subset of the figures produced by the
configured task list.

The manifest is the source of truth for which figure files documentation
depends on; the task list is the source of truth for how those files are
built. If the two drift apart --- a manifest path with no task to produce
it --- documentation builds can succeed locally (because the file is
already in the figure directory) but fail for any collaborator who has to
pull and regenerate from scratch. This module catches that drift.
"""

from pathlib import Path
from typing import Sequence

from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.figures.figure_exporter import (
    ValidFigureMeta,
    get_figure_destination,
)
from mecfs_bio.figures.manifest import FigureManifest


class ManifestTaskMismatchError(ValueError):
    """Raised when the manifest references paths no task in the list produces."""


def validate_manifest_subset_of_tasks(
    manifest: FigureManifest,
    tasks: Sequence[Task],
    fig_dir: Path,
) -> None:
    """
    Verify that every path in ``manifest`` is covered by ``tasks``.

    A manifest path is covered when either:
    - it matches exactly the destination of a file-emitting task, or
    - it lies inside the destination directory of a directory-emitting task.

    Raises ManifestTaskMismatchError listing any uncovered paths.
    """
    file_destinations: set[str] = set()
    dir_destinations: list[str] = []
    for task in tasks:
        meta = task.meta
        assert isinstance(meta, ValidFigureMeta)
        dst = get_figure_destination(meta=meta, fig_dir=fig_dir)
        rel = dst.relative_to(fig_dir).as_posix()
        if isinstance(meta, GWASPlotDirectoryMeta):
            dir_destinations.append(rel)
        else:
            file_destinations.add(rel)

    orphan_paths: list[str] = []
    for manifest_path in manifest.figures:
        if manifest_path in file_destinations:
            continue
        if any(
            manifest_path == d or manifest_path.startswith(d + "/")
            for d in dir_destinations
        ):
            continue
        orphan_paths.append(manifest_path)

    if orphan_paths:
        raise ManifestTaskMismatchError(
            "Figure manifest references paths that no task in the supplied "
            "task list produces. Either add a task that generates them or "
            "remove them from the manifest. Offending paths: "
            f"{sorted(orphan_paths)}"
        )
