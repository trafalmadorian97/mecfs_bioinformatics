"""
Point documentation figure references at the vectorized S-LDSC figures.

For every task in ALL_FIGURE_TASKS whose asset id ends in the vectorized S-LDSC suffix (_vec),
the figure it exports is named after that id. The figure exported by the corresponding original
S-LDSC task has the same name without the suffix. This script rewrites each reference to the old
name under docs/ (as docs/_figs/<old> or ../_figs/<old>) to the new name.

Run it after switching a generator to vectorized S-LDSC and exporting the new figures. Pages
whose figures still come from the original S-LDSC tasks are left alone, since their figure tasks
carry no suffix.

A reference matches only when the old name is followed by a path separator, quote, closing
parenthesis, whitespace or end of line, so an old name that is a prefix of a longer name is not
touched.

Usage:
    pixi r python experiments/claude/vectorized_sldsc/rename_sldsc_doc_figure_refs.py           # dry run
    pixi r python experiments/claude/vectorized_sldsc/rename_sldsc_doc_figure_refs.py --apply   # edit files
"""

import argparse
import re
from pathlib import Path
from typing import Mapping, Sequence

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task_generator.sldsc_task_generator import VEC_ID_SUFFIX
from mecfs_bio.figures.fig_constants import DOCS_DIRECTORY, FIGURE_DIRECTORY
from mecfs_bio.figures.figure_exporter import ValidFigureMeta, get_figure_destination
from mecfs_bio.figures.figure_task_list import ALL_FIGURE_TASKS


def vectorized_figure_renames(tasks: Sequence[Task]) -> dict[str, str]:
    """
    Map the old figure name (file or directory name under the figure directory) to the new one,
    for every figure task whose asset id carries the vectorized suffix.
    """
    renames: dict[str, str] = {}
    for task in tasks:
        asset_id = str(task.asset_id)
        if not asset_id.endswith(VEC_ID_SUFFIX):
            continue
        meta = task.meta
        assert isinstance(meta, ValidFigureMeta), f"{asset_id} is not a figure task"
        new_name = get_figure_destination(meta, FIGURE_DIRECTORY).name
        assert new_name.startswith(asset_id), (
            f"figure name {new_name} does not start with asset id {asset_id}"
        )
        old_name = asset_id.removesuffix(VEC_ID_SUFFIX) + new_name[len(asset_id) :]
        renames[old_name] = new_name
    return renames


def rewrite_references(text: str, renames: Mapping[str, str]) -> tuple[str, int]:
    """Return the rewritten text and the number of references replaced."""
    total = 0
    for old_name, new_name in renames.items():
        pattern = re.compile(
            re.escape(f"{FIGURE_DIRECTORY.name}/{old_name}") + r"(?=[/\"')\s]|$)",
            flags=re.MULTILINE,
        )
        text, count = pattern.subn(f"{FIGURE_DIRECTORY.name}/{new_name}", text)
        total += count
    return text, total


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true", help="write the edits")
    args = parser.parse_args()

    renames = vectorized_figure_renames(ALL_FIGURE_TASKS)
    print(f"{len(renames)} vectorized S-LDSC figures in ALL_FIGURE_TASKS")
    for path in sorted(DOCS_DIRECTORY.rglob("*.md")):
        text = path.read_text(encoding="utf-8")
        new_text, count = rewrite_references(text, renames)
        if count == 0:
            continue
        print(f"{path}: {count} reference(s)")
        if args.apply:
            path.write_text(new_text, encoding="utf-8")
    if not args.apply:
        print("Dry run; pass --apply to write the edits.")


if __name__ == "__main__":
    main()
