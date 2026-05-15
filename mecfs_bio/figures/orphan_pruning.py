"""
Remove figure-manifest entries that no figure task produces, deleting their
local files along the way.

This module backs the "I removed a figure task and want the figure system
to catch up" flow. We never prune a path that is still referenced from
documentation --- that would silently break the docs build. Instead,
``prune_orphan_figures`` collects the conflicts and raises, leaving the
manifest and figure directory untouched until the user resolves them
(either by restoring the task to ``ALL_FIGURE_TASKS`` or by removing the
doc references).
"""

import shutil
from pathlib import Path
from typing import Sequence

import structlog

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.figures.manifest import FigureManifest
from mecfs_bio.figures.manifest_validation import find_orphan_paths

logger = structlog.get_logger()

_DOC_EXTENSIONS = (".md", ".mdx")


class OrphanReferencedInDocsError(ValueError):
    """Raised when an orphan manifest path is still referenced by documentation."""


def find_doc_references(rel_path: str, docs_dir: Path) -> list[Path]:
    """
    Return every ``.md`` / ``.mdx`` file under ``docs_dir`` whose contents
    mention ``rel_path`` as a substring.

    The figure paths in the manifest are POSIX-style and relative to the
    figure directory (``docs/_figs``). Every reference style in the repo
    --- ``plotly_embed("docs/_figs/<rel_path>")``,
    ``include_file("docs/_figs/<rel_path>")``,
    ``![alt](../../_figs/<rel_path>)``, and ``<iframe src="..._figs/<rel_path>">``
    --- contains the rel_path verbatim, so a substring search is reliable.
    """
    hits: list[Path] = []
    for path in sorted(docs_dir.rglob("*")):
        if not path.is_file() or path.suffix not in _DOC_EXTENSIONS:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        if rel_path in text:
            hits.append(path)
    return hits


def prune_orphan_figures(
    manifest_path: Path,
    docs_dir: Path,
    fig_dir: Path,
    figure_tasks: Sequence[Task],
) -> None:
    """
    Drop manifest entries that no task in ``figure_tasks`` produces, and
    delete the corresponding files under ``fig_dir``.

    For each orphan, ``docs_dir`` is scanned for references. If any orphan
    is still referenced the function raises ``OrphanReferencedInDocsError``
    without mutating disk; the error distinguishes orphans whose local
    file is also missing, since those require restoring the task before
    the figure can be re-pulled.
    """
    manifest = FigureManifest.load(manifest_path)
    orphans = find_orphan_paths(manifest=manifest, tasks=figure_tasks, fig_dir=fig_dir)
    if not orphans:
        logger.debug("No orphan figure-manifest entries to prune.")
        return

    blockers_present: dict[str, list[Path]] = {}
    blockers_missing: dict[str, list[Path]] = {}
    safe_to_remove: list[str] = []
    for rel_path in orphans:
        refs = find_doc_references(rel_path=rel_path, docs_dir=docs_dir)
        if refs:
            if (fig_dir / rel_path).exists():
                blockers_present[rel_path] = refs
            else:
                blockers_missing[rel_path] = refs
        else:
            safe_to_remove.append(rel_path)

    if blockers_present or blockers_missing:
        raise OrphanReferencedInDocsError(
            _format_blocked_message(
                blockers_present=blockers_present,
                blockers_missing=blockers_missing,
            )
        )

    new_figures = dict(manifest.figures)
    for rel_path in safe_to_remove:
        local = fig_dir / rel_path
        if local.is_file():
            local.unlink()
            logger.debug(f"Deleted orphan figure file {local}.")
        elif local.is_dir():
            shutil.rmtree(local)
            logger.debug(f"Deleted orphan figure directory {local}.")
        new_figures.pop(rel_path, None)

    FigureManifest(figures=new_figures).save(manifest_path)
    logger.info(
        f"Pruned {len(safe_to_remove)} orphan manifest "
        f"{'entry' if len(safe_to_remove) == 1 else 'entries'}: {safe_to_remove}"
    )


def _format_blocked_message(
    blockers_present: dict[str, list[Path]],
    blockers_missing: dict[str, list[Path]],
) -> str:
    lines = [
        "Cannot prune figure manifest: the following entries are not "
        "produced by any task in the supplied task list, but are still "
        "referenced in documentation. Either restore the task to the "
        "list, or remove the doc references.",
        "",
    ]
    for rel_path in sorted(blockers_present):
        lines.append(f"  {rel_path}:")
        lines.append("    referenced in:")
        for ref in blockers_present[rel_path]:
            lines.append(f"      - {ref}")
    for rel_path in sorted(blockers_missing):
        lines.append(f"  {rel_path}:")
        lines.append("    referenced in:")
        for ref in blockers_missing[rel_path]:
            lines.append(f"      - {ref}")
        lines.append(
            "    (local file is also missing; restoring the task is "
            "required to be able to re-pull this figure)"
        )
    return "\n".join(lines)
