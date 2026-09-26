"""Stacked explainability figure for a uniform SUSIE run against one or more
prior-weighted SUSIE runs at the same locus.

Panels, top to bottom, sharing the genomic-position x-axis:
  1. Manhattan (-log10 p), points colored by LD r^2 with the min-p lead
     variant, with local recombination rate (cM/Mb from the hg19 genetic map)
     on a secondary axis. Drawn from the uniform run, whose filtered variants
     and LD every prior-weighted run shares (a prior run requires its prior to
     cover every variant).
  2. PIP, uniform run (credible-set variants only).
  3. One PIP row per prior-weighted run, in the order given. All PIP rows share
     one y-scale so they are directly comparable, and each prior row carries a
     callout on each credible set's prior-boosted variant naming its key
     annotation families (from that run's contrast task callouts.parquet).
  4. Genes.

With a single prior run this is the polyfun-vs-uniform explainability figure;
with two it compares two priors (e.g. the precomputed PolyFun prior and a
trait-specific one) against the same uniform run.

Writes both explain_plot.png and explain_plot.svg. Inspired by
SusieStackPlotTask but independent of it.

Figures are built through matplotlib's object-oriented API (a directly
constructed Figure), not pyplot. That keeps rendering off any global backend
state: the process-wide interactive backend is never selected or mutated, the
output format is chosen per file extension by savefig (so the .svg is a true
vector file and the .png raster), and no figure is registered in pyplot's
global manager, so there is nothing to close to reclaim memory.
"""

from pathlib import Path, PurePath

import narwhals as nw
import numpy as np
import polars as pl
import textalloc as ta
from attrs import frozen
from matplotlib.figure import Figure

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (
    scan_dataframe_asset,
)
from mecfs_bio.build_system.meta.result_directory_meta import (
    ResultDirectoryMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import (
    DataProcessingPipe,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    CALLOUT_LABEL_COL,
    CALLOUT_PIP_PF_COL,
    CALLOUTS_FILENAME,
)
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    COMBINED_CS_FILENAME,
    FILTERED_GWAS_FILENAME,
    FILTERED_LD_FILENAME,
    PIP_COLUMN,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
    GENE_INFO_STRAND_COL,
    GWAS_SIGNIFICANCE_MLOG10P,
    plot_gene_tracks,
    plot_susie_track,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.genetic_map_constants import (
    GMAP_POS_COL,
    GMAP_RATE_COL,
)
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)

PLOT_PNG_FILENAME = "explain_plot.png"
PLOT_SVG_FILENAME = "explain_plot.svg"
# PIP units of empty space kept above the tallest stem for the callout labels.
_PIP_LABEL_HEADROOM = 0.4
# Callout font sizes as (single-line, wrapped-multi-line) pairs. The larger pair
# is used by default; once the panel carries _CALLOUT_CROWDED_COUNT or more
# callouts the smaller pair keeps a dense locus from turning into a wall of
# overlapping labels.
_CALLOUT_FONT_DEFAULT = (10.0, 9.5)
_CALLOUT_FONT_CROWDED = (8.5, 8.0)
_CALLOUT_CROWDED_COUNT = 3
# Figure sizing. The panels stack vertically sharing the x-axis; the figure is
# sized well beyond matplotlib's defaults so the Manhattan points, PIP stems, and
# gene labels stay legible when the figure fills a laptop screen. Width is fixed;
# height scales with the panel count.
_FIGURE_WIDTH_IN = 18.0
_PANEL_HEIGHT_IN = 2.9
# Raster resolution
_PNG_DPI = 200


UNIFORM_ROW_LABEL = "uniform"


@frozen(slots=True)
class PriorRun:
    """One prior-weighted SUSIE run, drawn as its own PIP row.

    label names the row (its y-axis reads "PIP (label)"). contrast_task is the
    PolyfunExplainContrastTask explaining susie_task against the uniform run; its
    callouts.parquet supplies the row's callouts.
    """

    label: str
    susie_task: Task
    contrast_task: Task


@frozen(slots=True)
class PolyfunExplainPlotTask(Task):
    """Render the uniform-vs-prior explainability figure, one PIP row per prior
    run."""

    meta: Meta
    susie_uniform_task: Task
    prior_runs: tuple[PriorRun, ...]
    gene_info_task: Task
    genetic_map_task: Task
    genome_build: GenomeBuild = "19"
    gene_info_pipe: DataProcessingPipe = IdentityPipe()

    def __attrs_post_init__(self) -> None:
        assert len(self.prior_runs) >= 1, "Need at least one prior run to plot"
        labels = [UNIFORM_ROW_LABEL] + [run.label for run in self.prior_runs]
        assert len(set(labels)) == len(labels), (
            f"PIP row labels must be unique: {labels}"
        )

    @property
    def deps(self) -> list["Task"]:
        prior_deps = [
            task
            for run in self.prior_runs
            for task in (run.susie_task, run.contrast_task)
        ]
        return [
            self.susie_uniform_task,
            *prior_deps,
            self.gene_info_task,
            self.genetic_map_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        uni_dir = _dir(fetch, self.susie_uniform_task)
        locus = _load_locus(uni_dir)
        pip_rows = [_uniform_pip_row(uni_dir)] + [
            _prior_pip_row(fetch, run) for run in self.prior_runs
        ]
        genes = (
            self.gene_info_pipe.process(
                scan_dataframe_asset(
                    fetch(self.gene_info_task.asset_id), self.gene_info_task.meta
                )
            )
            .collect()
            .to_polars()
        )
        recomb = _load_recomb(
            fetch, self.genetic_map_task, locus.chrom, locus.bp_min, locus.bp_max
        )
        _render(
            scratch_dir=scratch_dir,
            locus=locus,
            pip_rows=pip_rows,
            genes=genes,
            recomb=recomb,
            genome_build=self.genome_build,
        )
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        susie_uniform_task: Task,
        prior_runs: tuple[PriorRun, ...],
        gene_info_task: Task,
        genetic_map_task: Task,
        genome_build: GenomeBuild = "19",
        gene_info_pipe: DataProcessingPipe = IdentityPipe(),
    ) -> "PolyfunExplainPlotTask":
        source_meta = susie_uniform_task.meta
        if not isinstance(source_meta, ResultDirectoryMeta):
            raise ValueError(f"Unknown meta for uniform susie task: {source_meta}")
        meta = ResultDirectoryMeta(
            id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=PurePath("analysis"),
        )
        return cls(
            meta=meta,
            susie_uniform_task=susie_uniform_task,
            prior_runs=prior_runs,
            gene_info_task=gene_info_task,
            genetic_map_task=genetic_map_task,
            genome_build=genome_build,
            gene_info_pipe=gene_info_pipe,
        )


@frozen(slots=True)
class _Locus:
    """The uniform run's filtered variants and LD, plus the locus window they
    span."""

    gwas: pl.DataFrame
    ld: np.ndarray
    chrom: int
    bp_min: int
    bp_max: int


@frozen(slots=True)
class _PipRow:
    """One PIP row: its credible-set variants and, for a prior run, the callouts
    placed above them."""

    label: str
    cs: pl.DataFrame
    callouts: pl.DataFrame | None


def _dir(fetch: Fetch, task: Task) -> Path:
    asset = fetch(task.asset_id)
    assert isinstance(asset, DirectoryAsset)
    return asset.path


def _load_locus(uni_dir: Path) -> _Locus:
    gwas = pl.read_parquet(uni_dir / FILTERED_GWAS_FILENAME)
    pos = gwas[GWASLAB_POS_COL].to_numpy()
    return _Locus(
        gwas=gwas,
        ld=np.load(uni_dir / FILTERED_LD_FILENAME),
        chrom=int(gwas[GWASLAB_CHROM_COL][0]),
        bp_min=int(pos.min()),
        bp_max=int(pos.max()),
    )


def _uniform_pip_row(uni_dir: Path) -> _PipRow:
    # Credible-set membership: the PIP rows plot only these variants, colored and
    # legended by credible set (empty frame if the run found none).
    return _PipRow(
        label=UNIFORM_ROW_LABEL,
        cs=pl.read_parquet(uni_dir / COMBINED_CS_FILENAME),
        callouts=None,
    )


def _prior_pip_row(fetch: Fetch, run: PriorRun) -> _PipRow:
    return _PipRow(
        label=run.label,
        cs=pl.read_parquet(_dir(fetch, run.susie_task) / COMBINED_CS_FILENAME),
        callouts=pl.read_parquet(_dir(fetch, run.contrast_task) / CALLOUTS_FILENAME),
    )


def _norm_sf(z: np.ndarray) -> np.ndarray:
    from scipy.stats import norm

    return norm.sf(z)


def _load_recomb(
    fetch: Fetch, task: Task, chrom: int, bp_min: int, bp_max: int
) -> pl.DataFrame:
    """Locus-windowed recombination rate (cM/Mb) from the hg19 genetic map."""
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


def _plot_recomb(ax, recomb: pl.DataFrame) -> None:
    bp = recomb[GMAP_POS_COL].to_numpy().astype(float)
    rate = recomb[GMAP_RATE_COL].to_numpy().astype(float)
    if len(bp) >= 1:
        (line,) = ax.plot(bp, rate, color="tab:red", alpha=0.4, linewidth=1)
        line.set_rasterized(True)
    # The recomb line is red, but its axis label/ticks stay black to match the
    # other panels' axis labels (the color only carries meaning on the line).
    ax.set_ylabel("cM/Mb")


def _render(
    scratch_dir: Path,
    locus: _Locus,
    pip_rows: list[_PipRow],
    genes: pl.DataFrame,
    recomb: pl.DataFrame,
    genome_build: GenomeBuild,
) -> None:
    # manhattan + one row per PIP run + genes, one panel each.
    n_panels = 1 + len(pip_rows) + 1
    fig = Figure(figsize=(_FIGURE_WIDTH_IN, _PANEL_HEIGHT_IN * n_panels))
    # Left column holds the tracks; the narrow right column is reserved for
    # legends/colorbars so nothing overlaps the data (mirrors SusieStackPlotTask).
    gs = fig.add_gridspec(
        nrows=n_panels,
        ncols=2,
        width_ratios=[1.0, 0.15],
        hspace=0.12,
        wspace=0.02,
    )
    ax0 = fig.add_subplot(gs[0, 0])
    axes = [ax0] + [fig.add_subplot(gs[i, 0], sharex=ax0) for i in range(1, n_panels)]

    _plot_manhattan(fig, gs, ax0, locus, recomb)
    _plot_pip_rows(
        fig,
        gs,
        axes[1:-1],
        pip_rows,
        xlims=(float(locus.bp_min), float(locus.bp_max)),
    )
    _plot_genes(
        axes[-1],
        genes,
        chrom=locus.chrom,
        bp_min=locus.bp_min,
        bp_max=locus.bp_max,
        genome_build=genome_build,
    )
    _tidy_shared_x(axes, bp_min=locus.bp_min, bp_max=locus.bp_max)

    fig.savefig(scratch_dir / PLOT_PNG_FILENAME, dpi=_PNG_DPI, bbox_inches="tight")
    fig.savefig(scratch_dir / PLOT_SVG_FILENAME, bbox_inches="tight")


def _plot_manhattan(fig: Figure, gs, ax0, locus: _Locus, recomb: pl.DataFrame) -> None:
    """Manhattan colored by LD with the lead (min-p) variant, with the
    recombination rate on a twin axis and their keys in the right column."""
    x = locus.gwas[GWASLAB_POS_COL].to_numpy()

    z = (locus.gwas[GWASLAB_BETA_COL] / locus.gwas[GWASLAB_SE_COL]).to_numpy()
    neglogp = -np.log10(2.0 * _norm_sf(np.abs(z)))
    lead = int(np.argmax(np.abs(z)))
    r2 = locus.ld[lead, :] ** 2
    sc = ax0.scatter(x, neglogp, c=r2, cmap="viridis", vmin=0, vmax=1, s=10)
    sc.set_rasterized(True)
    # Mark the lead (min-p) variant with a black triangle, matching
    # SusieStackPlotTask's Manhattan panel.
    ax0.scatter(x[lead], neglogp[lead], s=35, marker="^", c="black")
    # Dashed grey reference line at the default GWAS genome-wide significance
    # threshold (matches SusieStackPlotTask's Manhattan panel).
    ax0.axhline(
        y=GWAS_SIGNIFICANCE_MLOG10P,
        color="grey",
        linestyle="--",
        linewidth=1.5,
    )
    ax0.set_ylabel("-log10 p")
    # Right-column cell for panel 1: a shortened colorbar in the upper portion
    # (anchored right, clear of the recomb axis's ticks/label), and below it a
    # small legend telling the reader the red twin-axis line is recombination
    # rate.
    cbar_cell = fig.add_subplot(gs[0, 1])
    cbar_cell.axis("off")
    # Colorbar and recomb key share one aligned column in the RIGHT half of the
    # cell (the left half holds the cM/Mb twin-axis ticks and label): the colorbar
    # bar occupies [cbar_x0, cbar_x0+cbar_w] and the recomb line segment sits
    # directly beneath it over the same x-extent, so the two keys line up.
    cbar_x0, cbar_w = 0.52, 0.20
    cax = cbar_cell.inset_axes((cbar_x0, 0.42, cbar_w, 0.52))
    fig.colorbar(sc, cax=cax, label="r$^2$ w/ lead")

    ax0b = ax0.twinx()
    _plot_recomb(ax0b, recomb)
    cbar_cell.plot(
        [cbar_x0, cbar_x0 + cbar_w],
        [0.24, 0.24],
        transform=cbar_cell.transAxes,
        color="tab:red",
        alpha=0.6,
        linewidth=1.2,
    )
    cbar_cell.text(
        cbar_x0 + cbar_w / 2.0,
        0.14,
        "recomb rate",
        transform=cbar_cell.transAxes,
        ha="center",
        va="top",
        fontsize=7,
    )


def _plot_pip_rows(
    fig: Figure,
    gs,
    axes: list,
    pip_rows: list[_PipRow],
    xlims: tuple[float, float],
) -> None:
    """PIP rows as vertical stems restricted to credible-set variants, colored
    per credible set with a legend in the right column (reuses the stackplot's
    susie track), and callouts above each prior row's stems."""
    assert len(axes) == len(pip_rows)
    for i, row in enumerate(pip_rows):
        _plot_pip_panel(fig, gs, i + 1, axes[i], row.cs, f"PIP ({row.label})")
    # Share one y-scale across every PIP row so their stem heights are directly
    # comparable (PIP in [0, 1]; scale to the tallest stem, else full range).
    # Reserve a fixed band above the tallest stem for the callout labels; raise
    # every row equally so their data scale stays shared. A fixed (not
    # stem-proportional) headroom keeps the band just large enough for labels
    # even when a stem already reaches PIP ~1.
    label_top = _shared_pip_top([row.cs for row in pip_rows]) + _PIP_LABEL_HEADROOM
    for ax, row in zip(axes, pip_rows):
        ax.set_ylim(0.0, label_top)
        # PIP is a probability, so never draw ticks/labels above 1.0 even though
        # the panel extends higher to fit the callouts.
        _cap_pip_ticks(ax)
        # Place the callout labels in the headroom above the stems, angled off to
        # the side so the leader line is clearly distinct from a vertical stem.
        if row.callouts is not None:
            _place_callouts(
                ax,
                row.callouts,
                row.cs,
                xlims=xlims,
                ylims=(0.0, label_top),
            )


def _plot_genes(
    ax_gene,
    genes: pl.DataFrame,
    chrom: int,
    bp_min: int,
    bp_max: int,
    genome_build: GenomeBuild,
) -> None:
    # Genes: reuse the stackplot's lane-packed gene track. Filter to this
    # chromosome first (the reference lists every chromosome; the helper windows
    # only by position).
    genes_chrom = genes.filter(
        pl.col(GENE_INFO_CHROM_COL).cast(pl.String) == str(chrom)
    )
    plot_gene_tracks(
        ax=ax_gene,
        gene_df=genes_chrom,
        start_bp=bp_min,
        end_bp=bp_max,
        gene_start_col=GENE_INFO_START_COL,
        gene_end_col=GENE_INFO_END_COL,
        gene_name_col=GENE_INFO_NAME_COL,
        gene_strand_col=GENE_INFO_STRAND_COL,
    )
    ax_gene.set_ylabel("genes")
    ax_gene.set_xlabel(f"hg{genome_build} chr{chrom} position (bp)")


def _tidy_shared_x(axes: list, bp_min: int, bp_max: int) -> None:
    # Lock every panel to the locus window and tidy the shared x-axis: only the
    # bottom (genes) panel keeps tick labels; drop top+right spines throughout.
    # The gene panel additionally drops its left spine so it has no vertical
    # frame lines at all (matching the stackplot).
    axes[0].set_xlim(bp_min, bp_max)
    for ax in axes[:-1]:
        ax.tick_params(axis="x", which="both", labelbottom=False, bottom=False)
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[-1].spines["left"].set_visible(False)


def _shared_pip_top(cs_frames: list[pl.DataFrame]) -> float:
    """Common y-axis top for every PIP row: a little above the tallest stem
    across the runs, or the full [0, 1] range when none has a credible set."""
    tops = [
        float(df[PIP_COLUMN].to_numpy().max())
        for df in cs_frames
        if df.height > 0 and PIP_COLUMN in df.columns
    ]
    if not tops:
        return 1.0
    return min(1.0, max(tops) * 1.05)


def _cap_pip_ticks(ax) -> None:
    """Keep the PIP axis's ticks/labels within [0, 1]: the panel extends above 1
    to make room for callout labels, but PIP is a probability so a tick above 1
    would misread. Drops any auto tick outside [0, 1] (a no-op for low-PIP panels
    whose ticks already stay under 1)."""
    ax.set_yticks([t for t in ax.get_yticks() if -1e-9 <= t <= 1.0 + 1e-9])


def _plot_pip_panel(
    fig: Figure,
    gs,
    row: int,
    ax_pip,
    cs_df: pl.DataFrame,
    label: str,
) -> None:
    """Draw one PIP panel: credible-set-colored stems on the left-column axis and
    a per-credible-set legend in the reserved right-column cell. The legend is
    numbered by SUSIE's own credible-set index (L2 -> "CS 2") so it matches the
    contrast task's detailed table (cs_pf/cs_u), rather than positionally."""
    legend_ax = fig.add_subplot(gs[row, 1])
    legend_ax.axis("off")
    plot_susie_track(
        susie_cs_df=cs_df,
        ax_pip=ax_pip,
        pip_legend_ax=legend_ax,
        cs_label_mode="susie_index",
    )
    ax_pip.set_ylabel(label)


def _place_callouts(
    ax_pf,
    callouts: pl.DataFrame,
    stem_df: pl.DataFrame,
    xlims: tuple[float, float],
    ylims: tuple[float, float],
) -> None:
    """Annotate a prior run's PIP panel: one text label per callout row, anchored at
    (POS, pip_pf). textalloc treats every PIP stem (0 -> pip) as an obstacle line
    and places each label in the free space above them, angled north-east so the
    leader line meets the box corner and reads distinctly from the vertical stems.
    xlims/ylims are the panel's true axis limits (textalloc normalizes distances
    against them). Empty frame -> no-op."""
    if callouts.height == 0:
        return
    xs = callouts[GWASLAB_POS_COL].to_numpy().astype(float).tolist()
    ys = callouts[CALLOUT_PIP_PF_COL].to_numpy().astype(float).tolist()
    # A long (many-family) callout is wrapped one family per line, keeping its box
    # narrow so textalloc can seat it without crossing a neighbouring stem; short
    # callouts render single-line. The whole panel drops to the smaller font pair
    # once it is crowded (>= _CALLOUT_CROWDED_COUNT callouts).
    fonts = (
        _CALLOUT_FONT_CROWDED
        if callouts.height >= _CALLOUT_CROWDED_COUNT
        else _CALLOUT_FONT_DEFAULT
    )
    wrapped = [
        _wrap_callout_label(t, fonts) for t in callouts[CALLOUT_LABEL_COL].to_list()
    ]
    texts = [text for text, _ in wrapped]
    sizes = [size for _, size in wrapped]
    # Every stem (vertical line from 0 to its PIP) is an obstacle to route labels
    # and leader lines around.
    stem_x = stem_df[GWASLAB_POS_COL].to_numpy().astype(float)
    stem_pip = stem_df[PIP_COLUMN].to_numpy().astype(float)
    x_lines: list[np.ndarray | list[float]] = [[float(px), float(px)] for px in stem_x]
    y_lines: list[np.ndarray | list[float]] = [[0.0, float(pp)] for pp in stem_pip]
    # Seed so textalloc's candidate search is reproducible across builds (keeps
    # the committed SVG stable); textalloc draws from numpy's global RNG.
    np.random.seed(0)
    ta.allocate(
        ax_pf,
        xs,
        ys,
        texts,
        x_scatter=stem_x.tolist(),
        y_scatter=stem_pip.tolist(),
        x_lines=x_lines,
        y_lines=y_lines,
        textsize=sizes,
        linecolor="black",
        linewidth=0.6,
        direction=["northeast", "northwest", "north", "east"],
        min_distance=0.03,
        max_distance=0.8,
        xlims=xlims,
        ylims=ylims,
        avoid_label_lines_overlap=True,
        avoid_crossing_label_lines=True,
        nbr_candidates=5000,
    )


def _wrap_callout_label(
    label: str, fonts: tuple[float, float], max_chars: int = 45
) -> tuple[str, float]:
    """Lay out one callout label. A label longer than max_chars (a variant with
    several or verbose annotation families) is wrapped to one family per line, so
    its box is only as wide as the longest single family and can be placed clear of
    neighbouring stems. Shorter labels are returned unchanged. fonts is the
    (single-line, wrapped) size pair the caller chose for the panel's crowding.

    Parses the "pos:nea:ea (fam ++, fam +, ...)" form produced by the contrast
    task; anything not matching that shape is returned as-is."""
    single_line_size, wrapped_size = fonts
    if len(label) <= max_chars or " (" not in label or not label.endswith(")"):
        return label, single_line_size
    head, inner = label.split(" (", 1)
    families = inner[:-1].split(", ")
    return "\n".join([head, *families]), wrapped_size
