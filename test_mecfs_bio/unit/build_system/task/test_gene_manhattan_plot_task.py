import numpy as np
import pandas as pd
import plotly.graph_objs as go
import pytest

from mecfs_bio.build_system.task.gene_manhattan_plot_task import (
    build_manhattan_plot,
)


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


def test_build_manhattan_plot_hover_and_yvalues():
    df = _synthetic_df()
    fig = build_manhattan_plot(
        df=df,
        sig_threshold=5e-8,
        point_size=5,
        colors=("#1f77b4", "#ff7f0e"),
        sig_line_color="red",
        title=None,
    )

    chr1_trace = next(t for t in fig.data if t.name == "chr1")
    # Hover template surfaces gene name, Ensembl ID, and -log10(p).
    assert "%{customdata[0]}" in chr1_trace.hovertemplate
    assert "Ensembl" in chr1_trace.hovertemplate
    assert "log" in chr1_trace.hovertemplate

    # The y-values are -log10(p) for chr1's four genes.
    expected_y = -np.log10(np.array([0.5, 0.1, 1e-3, 1e-9]))
    np.testing.assert_allclose(sorted(chr1_trace.y), sorted(expected_y))


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
    )

    expected_y = float(-np.log10(0.05 / len(df)))
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
    )

    chr1_trace = next(t for t in fig.data if t.name == "chr1")
    # Original 4 chr1 points, no extras for the bad rows.
    assert len(chr1_trace.y) == 4


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
        )
