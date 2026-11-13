"""
Matplotlib helpers shared by the locus plots: Manhattan panels, recombination-rate tracks from a
genetic map, and gene-track panels.

-log10 p is the negative base-10 logarithm of a variant's association p-value. The recombination
rate is in cM/Mb (centimorgans per megabase) and comes from a genetic map table with the
GMAP_* columns, in whatever genome build the map was made for.
"""

import narwhals as nw
import numpy as np
import polars as pl
from matplotlib.axes import Axes

from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
    GENE_INFO_STRAND_COL,
    GWAS_SIGNIFICANCE_MLOG10P,
    plot_gene_tracks,
)
from mecfs_bio.constants.genetic_map_constants import GMAP_POS_COL, GMAP_RATE_COL
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_CHROM_NAME_FOR_CODE,
)


def load_recomb(
    fetch: Fetch, task: Task, chrom: int, bp_min: int, bp_max: int
) -> pl.DataFrame:
    """Locus-windowed recombination rate (cM/Mb) from a genetic map task."""
    return (
        scan_dataframe_asset(fetch(task.asset_id), task.meta)
        .filter(
            (nw.col(GWASLAB_CHROM_COL) == chrom)
            & (nw.col(GMAP_POS_COL) >= bp_min)
            & (nw.col(GMAP_POS_COL) <= bp_max)
        )
        .select(GWASLAB_CHROM_COL, GMAP_POS_COL, GMAP_RATE_COL)
        .collect()
        .to_polars()
        .sort(GMAP_POS_COL)
    )


def plot_recomb(ax: Axes, recomb: pl.DataFrame) -> None:
    bp = recomb[GMAP_POS_COL].to_numpy().astype(float)
    rate = recomb[GMAP_RATE_COL].to_numpy().astype(float)
    if len(bp) >= 1:
        (line,) = ax.plot(bp, rate, color="tab:red", alpha=0.4, linewidth=1)
        line.set_rasterized(True)
    # The recomb line is red, but its axis label/ticks stay black to match the
    # other panels' axis labels (the color only carries meaning on the line).
    ax.set_ylabel("cM/Mb")


def plot_association(
    ax: Axes,
    pos: np.ndarray,
    mlog10p: np.ndarray,
    label: str,
    significance_line: bool,
) -> int:
    """Scatter -log10 p against position, mark the lead variant, and return its index."""
    ax.scatter(pos, mlog10p, s=8, color="tab:blue").set_rasterized(True)
    lead = int(np.nanargmax(mlog10p))
    ax.scatter(pos[lead], mlog10p[lead], s=40, marker="^", color="black")
    if significance_line:
        ax.axhline(
            y=GWAS_SIGNIFICANCE_MLOG10P, color="grey", linestyle="--", linewidth=1
        )
    ax.set_ylabel(f"{label}\n-log10 p")
    return lead


def plot_gene_track_panel(
    ax: Axes, genes: pl.DataFrame, chrom: int, start: int, end: int
) -> None:
    """Gene tracks of the genes on chrom within [start, end], with a GRCh38 x-axis label."""
    chrom_name = GWASLAB_CHROM_NAME_FOR_CODE.get(chrom, None) or str(chrom)
    plot_gene_tracks(
        ax=ax,
        gene_df=genes.filter(pl.col(GENE_INFO_CHROM_COL).cast(pl.String) == chrom_name),
        start_bp=start,
        end_bp=end,
        gene_start_col=GENE_INFO_START_COL,
        gene_end_col=GENE_INFO_END_COL,
        gene_name_col=GENE_INFO_NAME_COL,
        gene_strand_col=GENE_INFO_STRAND_COL,
    )
    ax.set_ylabel("genes")
    ax.set_xlabel(f"GRCh38 chr{chrom_name} position (bp)")
