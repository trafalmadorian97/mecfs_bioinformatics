from pathlib import Path, PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_manhattan_plot_meta import (
    GWASLabManhattanQQPlotMeta,
)
from mecfs_bio.build_system.meta.markdown_file_meta import MarkdownFileMeta
from mecfs_bio.build_system.meta.plot_file_meta import GWASPlotFileMeta
from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.figures.figure_exporter import get_figure_destination


def test_get_figure_destination_plot_file_meta(tmp_path: Path):
    meta = GWASPlotFileMeta(
        trait="t",
        project="p",
        extension=".html",
        id=AssetId("my_plot"),
    )
    assert (
        get_figure_destination(meta=meta, fig_dir=tmp_path) == tmp_path / "my_plot.html"
    )


def test_get_figure_destination_manhattan_meta(tmp_path: Path):
    meta = GWASLabManhattanQQPlotMeta(
        trait="t",
        project="p",
        id=AssetId("my_manhattan"),
    )
    assert (
        get_figure_destination(meta=meta, fig_dir=tmp_path)
        == tmp_path / "my_manhattan.png"
    )


def test_get_figure_destination_plot_directory_meta(tmp_path: Path):
    meta = GWASPlotDirectoryMeta(
        trait="t",
        project="p",
        id=AssetId("my_dir"),
    )
    assert get_figure_destination(meta=meta, fig_dir=tmp_path) == tmp_path / "my_dir"


def test_get_figure_destination_markdown_meta(tmp_path: Path):
    meta = MarkdownFileMeta(
        id=AssetId("my_md"),
        trait="t",
        project="p",
        sub_dir=PurePath("analysis/markdown"),
    )
    assert get_figure_destination(meta=meta, fig_dir=tmp_path) == tmp_path / "my_md.mdx"
