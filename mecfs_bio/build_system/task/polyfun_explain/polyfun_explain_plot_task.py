"""Stacked explainability figure for a polyfun-vs-uniform SUSIE result.

Panels, top to bottom, sharing the genomic-position x-axis:
  1. Manhattan (-log10 p), points colored by LD r^2 with the min-p lead
     variant, with local recombination rate (cM/Mb from the hg19 genetic map)
     on a secondary axis.
  2..(1+x). One panel per important family: sum_c gamma_raw_c * a_ic across
     the locus (raw scaled value; the contrast is the profile, not plotted
     here).
  (2+x). Prior fold m*pi_i (log scale).
  (3+x). PIP, uniform run.
  (4+x). PIP, polyfun run.
  (5+x). Genes.

Writes both explain_plot.png and explain_plot.svg. Inspired by
SusieStackPlotTask but independent of it.

Figures are built through matplotlib's object-oriented API (a directly
constructed Figure), not pyplot. That keeps rendering off any global backend
state: the process-wide interactive backend is never selected or mutated, the
output format is chosen per file extension by savefig (so the .svg is a true
vector file and the .png raster), and no figure is registered in pyplot's
global manager, so there is nothing to close to reclaim memory.

The important families and the focal variant come from the contrast task's
selection.json, so the figure agrees with the contrast task's tables. The
per-variant family-scaled tracks are recomputed here (rather than read back
from the contrast task's output) directly from the ridge weights and the
annotation parquet, using the same allele-aware join and the same
_family_scaled helper the contrast task itself uses, so the values plotted
match the contrast task's family_scaled table exactly.
"""

import json
from pathlib import Path, PurePath

import narwhals as nw
import numpy as np
import polars as pl
from attrs import frozen
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    ANNOTATION_COL,
    GAMMA_RAW_COL,
)
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    FAMILY_COL as WEIGHTS_FAMILY_COL,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genetic_map.parse_genetic_map_task import (
    GMAP_POS_COL,
    GMAP_RATE_COL,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import (
    DataProcessingPipe,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    _ANNOT_KEY,
    FAMILY_COL,
    FAMILY_SCALED_COL,
    SELECTION_IMPORTANT_FAMILIES_KEY,
    SELECTION_JSON_FILENAME,
    _family_scaled,
    _load_annotations,
    _load_run_variants,
    _load_weights,
)
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    COMBINED_CS_FILENAME,
    FILTERED_GWAS_FILENAME,
    FILTERED_LD_FILENAME,
    PRIOR_FILENAME,
    PRIOR_WEIGHT_COLUMN,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
    GENE_INFO_STRAND_COL,
    plot_gene_tracks,
    plot_susie_track,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)

PLOT_PNG_FILENAME = "explain_plot.png"
PLOT_SVG_FILENAME = "explain_plot.svg"
# _ANNOT_KEY (imported from the contrast task) is the allele-aware join key
# (CHR, POS, unordered-allele-key); _load_annotations and _load_run_variants both
# attach that key, so the family-scaled values plotted here match the tables.


@frozen
class PolyfunExplainPlotTask(Task):
    """Render the 8-panel polyfun-vs-uniform explainability figure."""

    meta: Meta
    susie_uniform_task: Task
    susie_polyfun_task: Task
    contrast_task: Task
    annotation_parquet_task: Task
    gene_info_task: Task
    ridge_weights_task: Task
    genetic_map_task: Task
    gene_info_pipe: DataProcessingPipe = IdentityPipe()
    n_family_panels: int = 3

    @property
    def deps(self) -> list["Task"]:
        return [
            self.susie_uniform_task,
            self.susie_polyfun_task,
            self.contrast_task,
            self.ridge_weights_task,
            self.annotation_parquet_task,
            self.gene_info_task,
            self.genetic_map_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        uni_dir = _dir(fetch, self.susie_uniform_task)
        pf_dir = _dir(fetch, self.susie_polyfun_task)
        contrast_dir = _dir(fetch, self.contrast_task)
        selection = json.loads((contrast_dir / SELECTION_JSON_FILENAME).read_text())
        families = selection[SELECTION_IMPORTANT_FAMILIES_KEY][: self.n_family_panels]

        pf_variants = _load_run_variants(pf_dir).sort(GWASLAB_POS_COL)
        uni_variants = _load_run_variants(uni_dir)
        # Credible-set membership per run: the PIP panels plot only these variants,
        # colored/legended by credible set (empty frame if the run found none).
        uni_cs = pl.read_parquet(uni_dir / COMBINED_CS_FILENAME)
        pf_cs = pl.read_parquet(pf_dir / COMBINED_CS_FILENAME)
        pf_full = pl.read_parquet(pf_dir / FILTERED_GWAS_FILENAME)
        ld = np.load(pf_dir / FILTERED_LD_FILENAME)
        prior_w = pl.read_parquet(pf_dir / PRIOR_FILENAME)[
            PRIOR_WEIGHT_COLUMN
        ].to_numpy()

        chrom = int(pf_variants[GWASLAB_CHROM_COL][0])
        bp_min = int(pf_variants[GWASLAB_POS_COL].to_numpy().min())
        bp_max = int(pf_variants[GWASLAB_POS_COL].to_numpy().max())

        weights = _load_weights(fetch, self.ridge_weights_task)
        annot_cols = weights[ANNOTATION_COL].to_list()
        annot = _load_annotations(
            fetch, self.annotation_parquet_task, chrom, bp_min, bp_max, annot_cols
        )
        # Allele-aware join on (CHR, POS, allele-key) -- mirrors the contrast
        # task's own join exactly, so the family-scaled values plotted here match
        # its per_family_contrast/family_scaled tables.
        pf_annot = pf_variants.join(annot, on=_ANNOT_KEY, how="inner")
        gamma = dict(zip(weights[ANNOTATION_COL], weights[GAMMA_RAW_COL]))
        family = dict(zip(weights[ANNOTATION_COL], weights[WEIGHTS_FAMILY_COL]))
        scaled = _family_scaled(pf_annot, annot_cols, gamma, family)

        genes = (
            self.gene_info_pipe.process(
                scan_dataframe_asset(
                    fetch(self.gene_info_task.asset_id), self.gene_info_task.meta
                )
            )
            .collect()
            .to_polars()
        )

        _render(
            scratch_dir=scratch_dir,
            pf_variants=pf_variants,
            uni_variants=uni_variants,
            uni_cs=uni_cs,
            pf_cs=pf_cs,
            pf_full=pf_full,
            ld=ld,
            prior_w=prior_w,
            families=families,
            scaled=scaled,
            genes=genes,
            fetch=fetch,
            genetic_map_task=self.genetic_map_task,
            chrom=chrom,
            bp_min=bp_min,
            bp_max=bp_max,
        )
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        susie_uniform_task: Task,
        susie_polyfun_task: Task,
        contrast_task: Task,
        annotation_parquet_task: Task,
        gene_info_task: Task,
        ridge_weights_task: Task,
        genetic_map_task: Task,
        gene_info_pipe: DataProcessingPipe = IdentityPipe(),
        n_family_panels: int = 3,
    ) -> "PolyfunExplainPlotTask":
        source_meta = susie_polyfun_task.meta
        if not isinstance(source_meta, ResultDirectoryMeta):
            raise ValueError(f"Unknown meta for polyfun susie task: {source_meta}")
        meta = ResultDirectoryMeta(
            id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=PurePath("analysis"),
        )
        return cls(
            meta=meta,
            susie_uniform_task=susie_uniform_task,
            susie_polyfun_task=susie_polyfun_task,
            contrast_task=contrast_task,
            annotation_parquet_task=annotation_parquet_task,
            gene_info_task=gene_info_task,
            ridge_weights_task=ridge_weights_task,
            genetic_map_task=genetic_map_task,
            gene_info_pipe=gene_info_pipe,
            n_family_panels=n_family_panels,
        )


def _dir(fetch: Fetch, task: Task) -> Path:
    asset = fetch(task.asset_id)
    assert isinstance(asset, DirectoryAsset)
    return asset.path


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
    ax.set_ylabel("cM/Mb", color="tab:red")
    ax.tick_params(axis="y", labelcolor="tab:red")


def _render(
    scratch_dir: Path,
    pf_variants: pl.DataFrame,
    uni_variants: pl.DataFrame,
    uni_cs: pl.DataFrame,
    pf_cs: pl.DataFrame,
    pf_full: pl.DataFrame,
    ld: np.ndarray,
    prior_w: np.ndarray,
    families: list[str],
    scaled: pl.DataFrame,
    genes: pl.DataFrame,
    fetch: Fetch,
    genetic_map_task: Task,
    chrom: int,
    bp_min: int,
    bp_max: int,
) -> None:
    # manhattan + families + prior fold + 2 pip + genes, one panel each.
    n_panels = 1 + len(families) + 3 + 1
    fig = Figure(figsize=(12, 1.7 * n_panels))
    # Left column holds the tracks; the narrow right column is reserved for
    # legends/colorbars so nothing overlaps the data (mirrors SusieStackPlotTask).
    gs = fig.add_gridspec(
        nrows=n_panels,
        ncols=2,
        width_ratios=[1.0, 0.10],
        hspace=0.12,
        wspace=0.02,
    )
    ax0 = fig.add_subplot(gs[0, 0])
    axes = [ax0] + [fig.add_subplot(gs[i, 0], sharex=ax0) for i in range(1, n_panels)]
    x = pf_full[GWASLAB_POS_COL].to_numpy()

    # Panel 1: Manhattan colored by LD with lead (min-p) variant + recomb rate.
    z = (pf_full[GWASLAB_BETA_COL] / pf_full[GWASLAB_SE_COL]).to_numpy()
    neglogp = -np.log10(2.0 * _norm_sf(np.abs(z)))
    lead = int(np.argmax(np.abs(z)))
    r2 = ld[lead, :] ** 2
    sc = ax0.scatter(x, neglogp, c=r2, cmap="viridis", vmin=0, vmax=1, s=10)
    sc.set_rasterized(True)
    ax0.set_ylabel("-log10 p")
    # Colorbar sits at the RIGHT edge of the reserved right-column cell (a thin
    # inset), leaving the cell's left half for the recomb axis's ticks/label so
    # the two do not collide.
    cbar_cell = fig.add_subplot(gs[0, 1])
    cbar_cell.axis("off")
    cax = inset_axes(cbar_cell, width="30%", height="100%", loc="center right")
    fig.colorbar(sc, cax=cax, label="r$^2$ w/ lead")

    recomb = _load_recomb(fetch, genetic_map_task, chrom, bp_min, bp_max)
    ax0b = ax0.twinx()
    _plot_recomb(ax0b, recomb)

    # Family panels: sum_c gamma_raw_c * a_ic across the locus.
    for k, fam in enumerate(families):
        axk = axes[1 + k]
        fam_series = scaled.filter(pl.col(FAMILY_COL) == fam).sort(GWASLAB_POS_COL)
        (line,) = axk.plot(
            fam_series[GWASLAB_POS_COL].to_numpy(),
            fam_series[FAMILY_SCALED_COL].to_numpy(),
            linewidth=0.8,
        )
        line.set_rasterized(True)
        axk.set_ylabel(fam, fontsize=8)

    # Prior fold (log): the polyfun prior's multiplicative lift over uniform.
    lift = pf_variants.height * prior_w / prior_w.sum()
    ax_lift = axes[1 + len(families)]
    lift_sc = ax_lift.scatter(pf_full[GWASLAB_POS_COL].to_numpy(), lift, s=6)
    lift_sc.set_rasterized(True)
    ax_lift.set_yscale("log")
    ax_lift.set_ylabel("prior fold")

    # PIP uniform / polyfun as vertical stems restricted to credible-set variants,
    # colored per credible set with a legend in the right column (reuses the
    # stackplot's susie track). Empty frame -> no stems, no legend.
    _plot_pip_panel(fig, gs, 2 + len(families), axes, uni_cs, "PIP uniform")
    _plot_pip_panel(fig, gs, 3 + len(families), axes, pf_cs, "PIP polyfun")

    # Genes: reuse the stackplot's lane-packed gene track. Filter to this
    # chromosome first (the reference lists every chromosome; the helper windows
    # only by position).
    ax_gene = axes[-1]
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
    ax_gene.set_xlabel(f"chr{chrom} position (bp)")

    # Lock every panel to the locus window and tidy the shared x-axis: only the
    # bottom (genes) panel keeps tick labels; drop top+right spines throughout.
    # The gene panel additionally drops its left spine so it has no vertical
    # frame lines at all (matching the stackplot).
    ax0.set_xlim(bp_min, bp_max)
    for ax in axes[:-1]:
        ax.tick_params(axis="x", which="both", labelbottom=False, bottom=False)
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    ax_gene.spines["left"].set_visible(False)

    fig.savefig(scratch_dir / PLOT_PNG_FILENAME, dpi=150, bbox_inches="tight")
    fig.savefig(scratch_dir / PLOT_SVG_FILENAME, bbox_inches="tight")


def _plot_pip_panel(
    fig: Figure,
    gs,
    row: int,
    axes: list,
    cs_df: pl.DataFrame,
    label: str,
) -> None:
    """Draw one PIP panel: credible-set-colored stems on the left-column axis and
    a per-credible-set legend in the reserved right-column cell."""
    ax_pip = axes[row]
    legend_ax = fig.add_subplot(gs[row, 1])
    legend_ax.axis("off")
    plot_susie_track(susie_cs_df=cs_df, ax_pip=ax_pip, pip_legend_ax=legend_ax)
    ax_pip.set_ylim(bottom=0.0)
    ax_pip.set_ylabel(label)
