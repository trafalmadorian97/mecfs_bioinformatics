"""
Classify each DecodeME figure that find_stale_figures reported as byte-changed: did its
content change, or only incidental bytes (plotly div ids, PNG metadata, SVG clip-path ids)?

Compares the published figures in docs/_figs with the export in
find_stale_figures.EXPORT_DIR:

- .mdx: text equality;
- .parquet: polars frame equality;
- .png: pixel equality;
- .html and .svg: text equality after blanking random ids.

Usage:
    pixi r python -m experiments.claude.decode_me_panel_switch.classify_figure_changes \
        2>&1 | tee experiments/claude/decode_me_panel_switch/classify_figure_changes.log
"""

import re
from pathlib import Path, PurePath

import numpy as np
import polars as pl
from PIL import Image

from experiments.claude.decode_me_panel_switch.find_stale_figures import EXPORT_DIR
from mecfs_bio.figures.fig_constants import FIGURE_DIRECTORY
from mecfs_bio.figures.manifest import scan_figure_dir, sha256_of_file

_UUID = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}")
_SVG_ID = re.compile(r'(id="|url\(#|href="#)[^")]+')


def _normalized_text(path: Path) -> str:
    return _SVG_ID.sub(r"\1X", _UUID.sub("X", path.read_text()))


def content_equal(old: Path, new: Path) -> bool:
    suffix = old.suffix
    if suffix == ".parquet":
        return pl.read_parquet(old).equals(pl.read_parquet(new))
    if suffix == ".png":
        return np.array_equal(
            np.asarray(Image.open(old).convert("RGBA")),
            np.asarray(Image.open(new).convert("RGBA")),
        )
    if suffix in (".html", ".svg"):
        return _normalized_text(old) == _normalized_text(new)
    return old.read_text() == new.read_text()


def main() -> None:
    exported = scan_figure_dir(EXPORT_DIR).figures
    byte_changed: list[PurePath] = [
        path
        for path, sha in sorted(exported.items())
        if sha256_of_file(FIGURE_DIRECTORY / path) != sha
    ]
    content_changed = [
        path
        for path in byte_changed
        if not content_equal(FIGURE_DIRECTORY / path, EXPORT_DIR / path)
    ]
    bytes_only = [path for path in byte_changed if path not in content_changed]
    print(f"{len(bytes_only)} changed in bytes only:")
    for path in bytes_only:
        print(f"  {path}")
    print(f"\n{len(content_changed)} changed in content:")
    for path in content_changed:
        print(f"  {path}")


if __name__ == "__main__":
    main()
