"""
Regenerate the DecodeME figures whose content changed after the switch to the gnomAD panel.

A figure task is stale if any file it publishes differs in content (classify_figure_changes)
between docs/_figs and the export made by find_stale_figures. For each stale task, its
manifest entries and local files are deleted, and generate_new_figures then rebuilds just
those tasks into docs/_figs.

The manifest is left without the stale entries; push_figures (or invoke publish-figures)
re-adds them with the new hashes and uploads the blobs.

Run after find_stale_figures.

Usage:
    pixi r python -m experiments.claude.decode_me_panel_switch.refresh_decode_me_figures \
        2>&1 | tee experiments/claude/decode_me_panel_switch/refresh_decode_me_figures.log
"""

import shutil
from pathlib import Path, PurePath

from experiments.claude.decode_me_panel_switch.classify_figure_changes import (
    content_equal,
)
from experiments.claude.decode_me_panel_switch.find_stale_figures import (
    EXPORT_DIR,
    panel_dependent_decode_me_figure_tasks,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.figures.fig_constants import FIGURE_DIRECTORY, FIGURE_MANIFEST_PATH
from mecfs_bio.figures.figure_exporter import ValidFigureMeta, get_figure_destination
from mecfs_bio.figures.key_scripts.generate_figures import FIGURE_EXPORTER
from mecfs_bio.figures.key_scripts.generate_new_figures import generate_new_figures
from mecfs_bio.figures.manifest import FigureManifest


def published_paths(task: Task) -> list[PurePath]:
    """Paths, relative to the figure directory, of the files a figure task publishes."""
    meta = task.meta
    assert isinstance(meta, ValidFigureMeta)
    exported = get_figure_destination(meta=meta, fig_dir=EXPORT_DIR)
    files = sorted(exported.rglob("*")) if exported.is_dir() else [exported]
    return [PurePath(f.relative_to(EXPORT_DIR)) for f in files if f.is_file()]


def is_stale(task: Task) -> bool:
    return any(
        not (FIGURE_DIRECTORY / path).exists()
        or not content_equal(FIGURE_DIRECTORY / path, EXPORT_DIR / path)
        for path in published_paths(task)
    )


def remove_published_figure(task: Task, manifest: FigureManifest) -> FigureManifest:
    for path in published_paths(task):
        manifest = manifest.without_entry(path)
    meta = task.meta
    assert isinstance(meta, ValidFigureMeta)
    destination: Path = get_figure_destination(meta=meta, fig_dir=FIGURE_DIRECTORY)
    if destination.is_dir():
        shutil.rmtree(destination)
    else:
        destination.unlink(missing_ok=True)
    return manifest


def main() -> None:
    stale = [task for task in panel_dependent_decode_me_figure_tasks() if is_stale(task)]
    print(f"{len(stale)} stale figure tasks:", flush=True)
    for task in stale:
        print(f"  {task.asset_id}", flush=True)
    manifest = FigureManifest.load(FIGURE_MANIFEST_PATH)
    for task in stale:
        manifest = remove_published_figure(task, manifest)
    manifest.save(FIGURE_MANIFEST_PATH)
    generate_new_figures(
        all_figure_tasks=stale, fig_dir=FIGURE_DIRECTORY, exporter=FIGURE_EXPORTER
    )
    print("regenerated", flush=True)


if __name__ == "__main__":
    main()
