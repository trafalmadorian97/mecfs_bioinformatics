"""
Describe how each regenerated DecodeME figure differs from its published version, to find the
statements in docs/ that the new figures invalidate.

The published versions are the manifest entries that refresh_decode_me_figures removed (the
manifest at git HEAD), downloaded from the figure release into PUBLISHED_DIR. The new versions
are in docs/_figs. Data-bearing files are summarized:

- markdown tables (.mdx, other than SuSiE variant lists): a line diff;
- SuSiE result tables: per credible set, its size and top variants, and how far the old and
  new sets overlap;
- PolyFun explainability and prior-comparison tables: per credible-set column, the same
  summary, plus the top variants by each PIP column;
- plotly HTML: the top points of each trace by y value.

Images (.png, .svg) are drawn from the same data and are not summarized.

Run after refresh_decode_me_figures.

Usage:
    pixi r python -m experiments.claude.decode_me_panel_switch.describe_figure_changes \
        2>&1 | tee experiments/claude/decode_me_panel_switch/describe_figure_changes.log
"""

import base64
import difflib
import json
import re
import subprocess
from pathlib import Path, PurePath

import numpy as np
import polars as pl

from mecfs_bio.constants.gh_constants import GH_REPO_NAME
from mecfs_bio.figures.fig_constants import (
    FIGURE_DIRECTORY,
    FIGURE_GITHUB_RELEASE_TAG,
    FIGURE_MANIFEST_PATH,
)
from mecfs_bio.figures.manifest import FigureManifest
from mecfs_bio.util.github_commands.upload_download import download_release_asset

PUBLISHED_DIR = Path("experiments/claude/decode_me_panel_switch/published_figs")
TOP_N = 5
SUSIE_RESULT_SUFFIX = "_parquet_result_table.parquet"
THRESHOLD_LINE = re.compile(r'"shapes":\[\{[^}]*\},"type":"line"[^}]*"y0":([0-9.]+)')
PLOTLY_CALL = re.compile(r"Plotly\.newPlot\(\s*\"[^\"]+\",\s*(\[.*?\]),\s*\{", re.S)


def removed_entries() -> dict[PurePath, str]:
    head = subprocess.run(
        ["git", "show", f"HEAD:{FIGURE_MANIFEST_PATH.as_posix()}"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    old = {PurePath(k): v for k, v in json.loads(head)["figures"].items()}
    current = FigureManifest.load(FIGURE_MANIFEST_PATH).figures
    return {path: sha for path, sha in old.items() if path not in current}


def fetch_published(entries: dict[PurePath, str]) -> None:
    for path, sha in entries.items():
        dest = PUBLISHED_DIR / path
        if not dest.exists():
            download_release_asset(
                release_tag=FIGURE_GITHUB_RELEASE_TAG,
                repo_name=GH_REPO_NAME,
                asset_name=sha,
                dest=dest,
            )


def variant_id(row: dict) -> str:
    chrom = row.get("CHR", row.get("chr"))
    pos = row.get("POS", row.get("pos"))
    ea = row.get("EA", row.get("ea"))
    nea = row.get("NEA", row.get("nea"))
    return f"{chrom}:{pos}:{nea}:{ea}"


def credible_sets(
    df: pl.DataFrame, cs_col: str, pip_col: str
) -> dict[str, list[tuple[str, float]]]:
    sets: dict[str, list[tuple[str, float]]] = {}
    for row in df.filter(pl.col(cs_col).is_not_null()).iter_rows(named=True):
        if row[cs_col] in (-1, "", None):
            continue
        sets.setdefault(str(row[cs_col]), []).append((variant_id(row), row[pip_col]))
    return {
        cs: sorted(members, key=lambda m: -m[1]) for cs, members in sorted(sets.items())
    }


def describe_sets(label: str, sets: dict[str, list[tuple[str, float]]]) -> None:
    print(f"    {label}: {len(sets)} credible sets")
    for cs, members in sets.items():
        top = ", ".join(f"{v} ({p:.3f})" for v, p in members[:TOP_N])
        print(f"      {cs}: {len(members)} variants; top: {top}")


def describe_overlap(
    old: dict[str, list[tuple[str, float]]], new: dict[str, list[tuple[str, float]]]
) -> None:
    for cs, members in old.items():
        old_ids = {v for v, _ in members}
        best = max(
            ((len(old_ids & {v for v, _ in m}), c) for c, m in new.items()),
            default=(0, None),
        )
        print(
            f"      old {cs} ({len(old_ids)}) best matches new {best[1]}: "
            f"{best[0]} shared"
        )


def compare_credible_set_columns(old: pl.DataFrame, new: pl.DataFrame) -> None:
    cs_cols = [c for c in new.columns if c == "cs" or c.startswith("cs_")]
    for cs_col in cs_cols:
        pip_col = "PIP" if cs_col == "cs" else cs_col.replace("cs_", "pip_")
        if pip_col not in new.columns or cs_col not in old.columns:
            continue
        print(f"  [{cs_col} / {pip_col}]")
        old_sets = credible_sets(old, cs_col, pip_col)
        new_sets = credible_sets(new, cs_col, pip_col)
        describe_sets("published", old_sets)
        describe_sets("new", new_sets)
        describe_overlap(old_sets, new_sets)
    for pip_col in [c for c in new.columns if c.startswith("pip_")]:
        for label, df in (("published", old), ("new", new)):
            top = df.sort(pip_col, descending=True).head(TOP_N)
            listed = ", ".join(
                f"{variant_id(r)} ({r[pip_col]:.3f})" for r in top.iter_rows(named=True)
            )
            print(f"  top by {pip_col} ({label}): {listed}")


def compare_annotation_table(old: pl.DataFrame, new: pl.DataFrame) -> None:
    def variants(df: pl.DataFrame) -> list[str]:
        return [c for c in df.columns if c.count(":") == 3]

    print(f"  variant columns (published): {variants(old)}")
    print(f"  variant columns (new):       {variants(new)}")


def plotly_traces(path: Path) -> list[dict]:
    match = PLOTLY_CALL.search(path.read_text())
    if match is None:
        return []
    return json.loads(match.group(1))


def decoded(values: object) -> list:
    """Plotly writes numeric arrays either as lists or as base64 typed arrays."""
    if isinstance(values, dict) and "bdata" in values:
        return np.frombuffer(
            base64.b64decode(values["bdata"]), dtype=np.dtype(values["dtype"])
        ).tolist()
    return values if isinstance(values, list) else []


def point_labels(trace: dict) -> list:
    if isinstance(trace.get("hovertext"), list):
        return trace["hovertext"]
    if isinstance(trace.get("customdata"), list):
        return [row[0] if isinstance(row, list) else row for row in trace["customdata"]]
    return decoded(trace.get("x"))


def top_points(path: Path) -> list[tuple[float, str]]:
    points = [
        (float(y), str(label))
        for trace in plotly_traces(path)
        for y, label in zip(decoded(trace.get("y")), point_labels(trace))
        if isinstance(y, (int, float))
    ]
    return sorted(points, key=lambda p: -p[0])[: 2 * TOP_N]


def bonferroni_genes(path: Path) -> set[str]:
    """Genes of a gene Manhattan plot above the plot's own significance line.

    The newer plots draw only genes with -log10(p) above 2, so the threshold cannot be
    recomputed from the number of plotted points.
    """
    match = THRESHOLD_LINE.search(path.read_text())
    assert match is not None, f"no significance line in {path}"
    threshold = float(match.group(1))
    return {
        str(label)
        for trace in plotly_traces(path)
        for y, label in zip(decoded(trace.get("y")), point_labels(trace))
        if y > threshold
    }


def compare_plotly(old: Path, new: Path) -> None:
    for label, path in (("published", old), ("new", new)):
        listed = "; ".join(f"{name[:40]} = {y:.4g}" for y, name in top_points(path))
        print(f"  top points ({label}): {listed}")
    if "gene_manhattan" in old.name:
        old_genes = bonferroni_genes(old)
        new_genes = bonferroni_genes(new)
        print(f"  Bonferroni-significant genes: {len(old_genes)} -> {len(new_genes)}")
        print(f"    lost: {sorted(old_genes - new_genes)}")
        print(f"    gained: {sorted(new_genes - old_genes)}")


def describe(path: PurePath) -> None:
    old = PUBLISHED_DIR / path
    new = FIGURE_DIRECTORY / path
    print(f"\n=== {path}")
    name = path.name
    if name.endswith(".mdx") and "convert_cs_to_markdown" not in name:
        diff = difflib.unified_diff(
            old.read_text().splitlines(), new.read_text().splitlines(), lineterm="", n=0
        )
        for line in list(diff)[2:]:
            print(f"  {line}")
    elif name.endswith(".parquet") and "per_variant_annotation" in name:
        compare_annotation_table(pl.read_parquet(old), pl.read_parquet(new))
    elif name.endswith(".parquet"):
        compare_credible_set_columns(pl.read_parquet(old), pl.read_parquet(new))
    elif name.endswith(".html"):
        compare_plotly(old, new)
    else:
        print("  (image or SuSiE variant list; see the matching table)")


def main() -> None:
    entries = removed_entries()
    print(f"{len(entries)} figures were regenerated")
    fetch_published(entries)
    for path in sorted(entries):
        describe(path)


if __name__ == "__main__":
    main()
