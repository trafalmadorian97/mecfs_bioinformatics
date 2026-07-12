from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import pytest

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.gene_manhattan_plot_task import (
    GenePValueTableSource,
    MagmaGeneSource,
    build_manhattan_plot,
)
from mecfs_bio.build_system.task.magma.magma_gene_analysis_task import (
    GENE_ANALYSIS_OUTPUT_STEM_NAME,
)
from mecfs_bio.constants.genomic_coordinate_constants import EXTENDED_MHC_BUILD_37


def _synthetic_df() -> pd.DataFrame:
    """Three chromosomes, four genes each, with one strongly-significant gene per chrom."""
    rows: list[dict] = []
    for chrom in ["1", "2", "X"]:
        for offset, p in enumerate([0.5, 0.1, 1e-3, 1e-9]):
            pos = 100_000 + offset * 50_000
            ensembl_id = f"ENSG{chrom}_{offset:03d}"
            rows.append(
                {
                    "chrom": chrom,
                    "pos": float(pos),
                    "ensembl_id": ensembl_id,
                    "gene_name": f"GENE_{chrom}_{offset}",
                    "p_value": p,
                }
            )
    return pd.DataFrame(rows)


def test_build_manhattan_plot_basic_shape():
    df = _synthetic_df()
    fig = build_manhattan_plot(
        df=df,
        sig_threshold=5e-8,
        point_size=5,
        colors=("#1f77b4", "#ff7f0e"),
        sig_line_color="red",
        title="test",
        genome_build="19",
    )

    assert isinstance(fig, go.Figure)
    # One Scattergl trace per chromosome.
    chroms = sorted(df["chrom"].unique().tolist())
    assert len(fig.data) == len(chroms)
    for trace in fig.data:
        assert isinstance(trace, go.Scattergl)

    # x-axis labels match the chromosome names in canonical order (1, 2, X).
    xaxis = fig.layout.xaxis
    assert tuple(xaxis.ticktext) == ("1", "2", "X")
    # Tick positions are monotonically increasing as we walk along chromosomes.
    tickvals = list(xaxis.tickvals)
    assert tickvals == sorted(tickvals)


def test_build_manhattan_plot_alternates_colors_across_chromosomes():
    df = _synthetic_df()
    blue, orange = "#1f77b4", "#ff7f0e"
    fig = build_manhattan_plot(
        df=df,
        sig_threshold=5e-8,
        point_size=5,
        colors=(blue, orange),
        sig_line_color="red",
        title=None,
        genome_build="19",
    )
    chr_to_color = {t.name: t.marker.color for t in fig.data}
    # Canonical order: chr1 first => blue, chr2 second => orange, chrX third => blue.
    assert chr_to_color["chr1"] == blue
    assert chr_to_color["chr2"] == orange
    assert chr_to_color["chrX"] == blue


def test_build_manhattan_plot_sig_line_present():
    df = _synthetic_df()
    sig_threshold = 5e-8
    fig = build_manhattan_plot(
        df=df,
        sig_threshold=sig_threshold,
        point_size=5,
        colors=("#1f77b4", "#ff7f0e"),
        sig_line_color="red",
        title=None,
        genome_build="19",
    )

    # add_hline produces a shape on the layout at y = -log10(threshold).
    expected_y = float(-np.log10(sig_threshold))
    shapes = fig.layout.shapes
    assert len(shapes) >= 1
    matching = [s for s in shapes if s.y0 == expected_y and s.y1 == expected_y]
    assert len(matching) == 1


def test_build_manhattan_plot_default_sig_threshold_is_bonferroni():
    df = _synthetic_df()
    fig = build_manhattan_plot(
        df=df,
        sig_threshold=None,
        point_size=5,
        colors=("#1f77b4", "#ff7f0e"),
        sig_line_color="red",
        title=None,
        genome_build="19",
    )

    expected_y = float(-np.log10(0.05 / len(df)))
    shapes = fig.layout.shapes
    matching = [
        s
        for s in shapes
        if np.isclose(float(s.y0), expected_y) and np.isclose(float(s.y1), expected_y)
    ]
    assert len(matching) == 1


def test_build_manhattan_plot_bonferroni_uses_correction_count():
    df = _synthetic_df()
    # Supplying a correction count larger than the plotted-row count (e.g. because
    # high-p-value genes were filtered out before plotting) must drive the
    # threshold, keeping the significance line invariant to that filtering.
    num_genes_for_correction = 5_000
    fig = build_manhattan_plot(
        df=df,
        sig_threshold=None,
        point_size=5,
        colors=("#1f77b4", "#ff7f0e"),
        sig_line_color="red",
        title=None,
        genome_build="19",
        num_genes_for_correction=num_genes_for_correction,
    )

    expected_y = float(-np.log10(0.05 / num_genes_for_correction))
    shapes = fig.layout.shapes
    matching = [
        s
        for s in shapes
        if np.isclose(float(s.y0), expected_y) and np.isclose(float(s.y1), expected_y)
    ]
    assert len(matching) == 1


def test_build_manhattan_plot_drops_nonpositive_and_null_pvalues():
    df = _synthetic_df()
    # Inject a row with p == 0 and a row with p == NaN. Both should be dropped
    # rather than cause a -log10 blowup.
    df = pd.concat(
        [
            df,
            pd.DataFrame(
                [
                    {
                        "chrom": "1",
                        "pos": 999_999.0,
                        "ensembl_id": "ENSG_BAD_ZERO",
                        "gene_name": "BAD_ZERO",
                        "p_value": 0.0,
                    },
                    {
                        "chrom": "1",
                        "pos": 1_000_001.0,
                        "ensembl_id": "ENSG_BAD_NAN",
                        "gene_name": "BAD_NAN",
                        "p_value": float("nan"),
                    },
                ]
            ),
        ],
        ignore_index=True,
    )

    fig = build_manhattan_plot(
        df=df,
        sig_threshold=5e-8,
        point_size=5,
        colors=("#1f77b4", "#ff7f0e"),
        sig_line_color="red",
        title=None,
        genome_build="19",
    )

    chr1_trace = next(t for t in fig.data if t.name == "chr1")
    # Original 4 chr1 points, no extras for the bad rows.
    assert len(chr1_trace.y) == 4


def _write_gene_pvalue_source(
    tmp_path: Path,
    max_p_value: float | None,
) -> tuple[GenePValueTableSource, Fetch]:
    """Build a GenePValueTableSource over four genes with p-values 0.5, 0.02, 1e-3, 1e-9.

    Returns the source plus a fetch callable that resolves its two task assets.
    """
    read_spec = DataFrameReadSpec(DataFrameParquetFormat())

    ensembl_ids = ["ENSG000", "ENSG001", "ENSG002", "ENSG003"]
    table_df = pd.DataFrame({"gene": ensembl_ids, "pval": [0.5, 0.02, 1e-3, 1e-9]})
    table_path = tmp_path / "table.parquet"
    table_df.to_parquet(table_path)

    loc_df = pd.DataFrame(
        {
            "ensembl_name": ensembl_ids,
            "chrom": ["1", "1", "2", "X"],
            "gene_start": [100_000.0, 200_000.0, 100_000.0, 100_000.0],
            "gene_end": [150_000.0, 250_000.0, 150_000.0, 150_000.0],
            "gene_name": ["GENE0", "GENE1", "GENE2", "GENE3"],
        }
    )
    loc_path = tmp_path / "loc.parquet"
    loc_df.to_parquet(loc_path)

    source = GenePValueTableSource(
        table_task=FakeTask(SimpleFileMeta("table", read_spec=read_spec)),
        gene_locations_task=FakeTask(SimpleFileMeta("loc", read_spec=read_spec)),
        gene_col="gene",
        p_col="pval",
        genome_build="19",
        max_p_value=max_p_value,
    )

    def fetch(asset_id: AssetId):
        if asset_id == "table":
            return FileAsset(table_path)
        if asset_id == "loc":
            return FileAsset(loc_path)
        raise ValueError(f"Unknown asset id {asset_id}")

    return source, fetch


def test_gene_pvalue_source_filters_by_max_p_value(tmp_path: Path):
    source, fetch = _write_gene_pvalue_source(tmp_path, max_p_value=0.01)
    data = source.load_df(fetch=fetch)
    # 0.5 and 0.02 are at or above 0.01 and dropped; 1e-3 and 1e-9 remain.
    assert sorted(data.df["p_value"].tolist()) == [1e-9, 1e-3]


def test_gene_pvalue_source_no_filter_when_max_p_value_none(tmp_path: Path):
    source, fetch = _write_gene_pvalue_source(tmp_path, max_p_value=None)
    data = source.load_df(fetch=fetch)
    assert len(data.df) == 4


def test_gene_pvalue_source_correction_count_is_pre_filter(tmp_path: Path):
    # All four genes have valid positive p-values, so the multiple-testing count
    # must be 4 regardless of how many rows survive the max_p_value filter.
    source, fetch = _write_gene_pvalue_source(tmp_path, max_p_value=0.01)
    data = source.load_df(fetch=fetch)
    assert len(data.df) == 2
    assert data.num_genes_for_correction == 4


def _write_magma_source(
    tmp_path: Path,
    max_p_value: float | None,
) -> tuple[MagmaGeneSource, Fetch]:
    """Build a MagmaGeneSource over four genes with p-values 0.5, 0.02, 1e-3, 1e-9.

    Returns the source plus a fetch callable that resolves its two task assets.
    """
    magma_dir = tmp_path / "magma"
    magma_dir.mkdir()
    genes_out = magma_dir / f"{GENE_ANALYSIS_OUTPUT_STEM_NAME}.genes.out"
    genes_out.write_text(
        "GENE CHR START STOP P\n"
        "ENSG000 1 100000 150000 0.5\n"
        "ENSG001 1 200000 250000 0.02\n"
        "ENSG002 2 100000 150000 1e-3\n"
        "ENSG003 X 100000 150000 1e-9\n"
    )

    thesaurus_df = pd.DataFrame(
        {
            "Gene stable ID": ["ENSG000", "ENSG001", "ENSG002", "ENSG003"],
            "Gene name": ["GENE0", "GENE1", "GENE2", "GENE3"],
        }
    )
    thesaurus_path = tmp_path / "thesaurus.parquet"
    thesaurus_df.to_parquet(thesaurus_path)

    source = MagmaGeneSource(
        magma_task=FakeTask(SimpleDirectoryMeta("magma")),
        gene_thesaurus_task=FakeTask(
            SimpleFileMeta(
                "thesaurus", read_spec=DataFrameReadSpec(DataFrameParquetFormat())
            )
        ),
        genome_build="19",
        max_p_value=max_p_value,
    )

    def fetch(asset_id: AssetId):
        if asset_id == "magma":
            return DirectoryAsset(magma_dir)
        if asset_id == "thesaurus":
            return FileAsset(thesaurus_path)
        raise ValueError(f"Unknown asset id {asset_id}")

    return source, fetch


def test_magma_source_filters_by_max_p_value(tmp_path: Path):
    source, fetch = _write_magma_source(tmp_path, max_p_value=0.01)
    data = source.load_df(fetch=fetch)
    # 0.5 and 0.02 are at or above 0.01 and dropped; 1e-3 and 1e-9 remain. The
    # correction count still reflects all four tested genes.
    assert sorted(data.df["p_value"].tolist()) == [1e-9, 1e-3]
    assert data.num_genes_for_correction == 4


def test_magma_source_no_filter_when_max_p_value_none(tmp_path: Path):
    source, fetch = _write_magma_source(tmp_path, max_p_value=None)
    data = source.load_df(fetch=fetch)
    assert len(data.df) == 4
    assert data.num_genes_for_correction == 4


def test_build_manhattan_plot_y_axis_start_sets_lower_bound():
    df = _synthetic_df()
    y_axis_start = float(-np.log10(0.1))
    fig = build_manhattan_plot(
        df=df,
        sig_threshold=5e-8,
        point_size=5,
        colors=("#1f77b4", "#ff7f0e"),
        sig_line_color="red",
        title=None,
        genome_build="19",
        y_axis_start=y_axis_start,
    )
    y_range = fig.layout.yaxis.range
    assert y_range is not None
    assert float(y_range[0]) == pytest.approx(y_axis_start)
    # The upper bound stays above the most significant point (-log10(1e-9)).
    assert float(y_range[1]) > 9.0


def test_build_manhattan_plot_no_y_range_without_start():
    df = _synthetic_df()
    fig = build_manhattan_plot(
        df=df,
        sig_threshold=5e-8,
        point_size=5,
        colors=("#1f77b4", "#ff7f0e"),
        sig_line_color="red",
        title=None,
        genome_build="19",
    )
    # Without y_axis_start the axis is left to auto-range.
    assert fig.layout.yaxis.range is None


def _synthetic_df_with_hla() -> pd.DataFrame:
    """One chr6 gene inside the extended MHC region and one outside it."""
    in_pos = float((EXTENDED_MHC_BUILD_37.start + EXTENDED_MHC_BUILD_37.end) / 2)
    out_pos = float(EXTENDED_MHC_BUILD_37.end + 1_000_000)
    return pd.DataFrame(
        [
            {
                "chrom": "6",
                "pos": in_pos,
                "ensembl_id": "ENSG_HLA",
                "gene_name": "HLA_GENE",
                "p_value": 1e-6,
            },
            {
                "chrom": "6",
                "pos": out_pos,
                "ensembl_id": "ENSG_OUT",
                "gene_name": "OUT_GENE",
                "p_value": 1e-4,
            },
        ]
    )


def test_build_manhattan_plot_marks_hla_region_genes():
    df = _synthetic_df_with_hla()
    fig = build_manhattan_plot(
        df=df,
        sig_threshold=5e-8,
        point_size=5,
        colors=("#1f77b4", "#ff7f0e"),
        sig_line_color="red",
        title=None,
        genome_build="19",
        hla_interval=EXTENDED_MHC_BUILD_37,
        hla_marker_symbol="diamond",
    )
    chr6_trace = next(t for t in fig.data if t.name == "chr6")
    # The in-region gene sorts before the out-of-region gene by position, so the
    # symbol array is (diamond, circle).
    assert tuple(chr6_trace.marker.symbol) == ("diamond", "circle")
    # A legend entry describing the HLA marker is added.
    legend_names = [t.name for t in fig.data if t.showlegend]
    assert "Extended HLA/MHC region" in legend_names


def test_build_manhattan_plot_no_hla_marking_by_default():
    df = _synthetic_df_with_hla()
    fig = build_manhattan_plot(
        df=df,
        sig_threshold=5e-8,
        point_size=5,
        colors=("#1f77b4", "#ff7f0e"),
        sig_line_color="red",
        title=None,
        genome_build="19",
    )
    chr6_trace = next(t for t in fig.data if t.name == "chr6")
    # No hla_interval => plain circle marker (scalar symbol, not an array).
    assert chr6_trace.marker.symbol is None or chr6_trace.marker.symbol == "circle"


def test_build_manhattan_plot_empty_df_raises():
    df = pd.DataFrame(
        {
            "chrom": pd.Series(dtype=str),
            "pos": pd.Series(dtype=float),
            "ensembl_id": pd.Series(dtype=str),
            "gene_name": pd.Series(dtype=str),
            "p_value": pd.Series(dtype=float),
        }
    )
    with pytest.raises(ValueError):
        build_manhattan_plot(
            df=df,
            sig_threshold=5e-8,
            point_size=5,
            colors=("#1f77b4", "#ff7f0e"),
            sig_line_color="red",
            title=None,
            genome_build="19",
        )
