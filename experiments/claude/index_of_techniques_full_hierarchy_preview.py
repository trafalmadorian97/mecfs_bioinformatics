"""Preview the Index of Techniques before/after showing the full folder hierarchy.

The material/tags plugin builds docs/Analysis/Index_of_Techniques.md by grouping
every page that carries a technique tag under that tag, using the page's
(disambiguated) title as the link text. docs_hooks.py derives that title from the
page's folder path.

This script reproduces that grouping offline so we can compare the link text
produced by the OLD two-level logic (Study, Trait) against the NEW full-hierarchy
logic (deepest-first folder chain, ending with the trait) that docs_hooks now
implements -- without running a full mkdocs build.

Run: pixi r python experiments/claude/index_of_techniques_full_hierarchy_preview.py
"""

from __future__ import annotations

import re
from pathlib import Path, PurePath

import yaml

import docs_hooks

DOCS = Path("docs")
TAGS_ALLOWED = [
    "LDSC",
    "S-LDSC",
    "CT-LDSC",
    "MAGMA",
    "H-MAGMA",
    "Mendelian randomization",
    "MiXeR",
    "SuSiE",
    "LCV",
]
H1 = re.compile(r"^\s*#\s+(.+?)\s*$", re.MULTILINE)
FRONTMATTER = re.compile(r"^---\n(.*?)\n---\n", re.DOTALL)


def _old_context(src_uri: str) -> str | None:
    """The previous two-level behaviour: only Study and Trait."""
    parts = PurePath(src_uri).parts
    if len(parts) < docs_hooks.MIN_ANALYSIS_PARTS or parts[0] != docs_hooks.ANALYSIS_ROOT:
        return None
    trait = docs_hooks._trait_display(parts[1])
    study = docs_hooks._humanize(parts[2])
    return f"{study}, {trait}"


def _page_short_title(text: str) -> str:
    match = H1.search(text)
    return match.group(1) if match else "(no H1)"


def _page_tags(text: str) -> list[str]:
    match = FRONTMATTER.search(text)
    if not match:
        return []
    meta = yaml.safe_load(match.group(1)) or {}
    return list(meta.get("tags", []) or [])


def _render(title_fn) -> str:
    by_tag: dict[str, list[str]] = {tag: [] for tag in TAGS_ALLOWED}
    for md in sorted(DOCS.rglob("*.md")):
        text = md.read_text(encoding="utf-8")
        src_uri = md.relative_to(DOCS).as_posix()
        context = title_fn(src_uri)
        if context is None:
            continue
        tags = _page_tags(text)
        short = _page_short_title(text)
        title = f"{short} ({context})"
        for tag in tags:
            if tag in by_tag:
                by_tag[tag].append(title)

    lines: list[str] = []
    for tag in TAGS_ALLOWED:
        titles = sorted(by_tag[tag])
        if not titles:
            continue
        lines.append(f"## {tag}")
        lines.extend(f"- {t}" for t in titles)
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    print("=" * 80)
    print("CURRENT (two-level: Study, Trait)")
    print("=" * 80)
    print(_render(_old_context))
    print("=" * 80)
    print("PROPOSED (full hierarchy, deepest-first, ending with Trait)")
    print("=" * 80)
    print(_render(docs_hooks._analysis_context))


if __name__ == "__main__":
    main()
