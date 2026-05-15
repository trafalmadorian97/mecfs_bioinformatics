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
from typing import Iterable, Sequence

from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.figures.figure_exporter import (
    ValidFigureMeta,
    get_figure_destination,
)
from mecfs_bio.figures.manifest import FigureManifest


class ManifestTaskMismatchError(ValueError):
    """Raised when the manifest references paths no task in the list produces."""


def find_orphan_paths(
    paths: Iterable[str],
    tasks: Sequence[Task],
    fig_dir: Path,
) -> list[str]:
    """
    Return the paths in ``paths`` that no task in ``tasks`` produces.

    A path is covered (i.e. *not* an orphan) when either:
    - it matches exactly the destination of a file-emitting task, or
    - it lies inside the destination directory of a directory-emitting task.

    The caller decides which paths to consider --- typically manifest keys
    for validation, or the union of manifest keys and on-disk files for
    pruning. Result is sorted for stable error messages and pruning order.
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
    for path in paths:
        if path in file_destinations:
            continue
        if any(path == d or path.startswith(d + "/") for d in dir_destinations):
            continue
        orphan_paths.append(path)

    return sorted(orphan_paths)


def validate_manifest_subset_of_tasks(
    manifest: FigureManifest,
    tasks: Sequence[Task],
    fig_dir: Path,
) -> None:
    """
    Raises ManifestTaskMismatchError if any manifest path is uncovered by ``tasks``.
    """
    orphan_paths = find_orphan_paths(
        paths=manifest.figures.keys(), tasks=tasks, fig_dir=fig_dir
    )
    if orphan_paths:
        raise ManifestTaskMismatchError(
            "Figure manifest references paths that no task in the supplied "
            "task list produces. Either add a task that generates them or "
            "remove them from the manifest. Offending paths: "
            f"{orphan_paths}"
        )
