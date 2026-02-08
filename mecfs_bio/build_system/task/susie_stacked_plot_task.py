import matplotlib.patheffects as pe
import xarray as xr
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from attrs import frozen
from h5py._hl import dims
# import pandas as pd
from matplotlib import gridspec

import polars as pl
from matplotlib.figure import Figure

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_mlog10p_pipe import ComputeMlog10pIfNeededPipe
from mecfs_bio.build_system.task.pipes.compute_p_from_beta_se import ComputePFromBetaSEPipeIfNeeded
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import FILTERED_GWAS_FILENAME, COMBINED_CS_FILENAME, \
    FILTERED_LD_FILENAME, PIP_COLUMN, CS_COLUMN
from mecfs_bio.build_system.task.susie_trackplot_task import EnsemblGeneInfoSource, RegionSelect, get_region
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import GWASLAB_POS_COL, GWASLAB_MLOG10P_COL
from mecfs_bio.util.plotting.save_fig import write_plots_to_dir

GENE_INFO_START_COL = "gene_start"
GENE_INFO_END_COL = "gene_end"
GENE_INFO_NAME_COL ="gene_name"
GENE_INFO_STRAND_COL = "strand"

_gwas_pipe = CompositePipe(
    [

        ComputePFromBetaSEPipeIfNeeded(),
        ComputeMlog10pIfNeededPipe()
    ]
)


@frozen
class SusieStackPlotTask(Task):
    _meta: Meta
    susie_task: Task
    gene_info_task: Task
    gene_info_pipe: DataProcessingPipe
    region_mode: RegionSelect


    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        deps =  [self.susie_task, self.gene_info_task]
        return deps

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:

        susie_asset = fetch(self.susie_task.asset_id)
        assert isinstance(susie_asset, DirectoryAsset)
        susie_dir= susie_asset.path

        gene_info_asset = fetch(self.gene_info_task.asset_id)
        assert isinstance(gene_info_asset, FileAsset)
        gene_info_df = self.gene_info_pipe.process( scan_dataframe_asset(gene_info_asset,self.gene_info_task.meta)).collect().to_polars()


        chrom, start, end = get_region(self.region_mode, susie_output_path=susie_dir)




        fig = plot_locus_tracks_matplotlib(
         gwas_df=  _gwas_pipe.process_eager_polars( pl.read_parquet(susie_dir/FILTERED_GWAS_FILENAME)),
            susie_cs_df= pl.read_parquet(susie_dir / COMBINED_CS_FILENAME),
            ld_np= np.load(susie_dir/FILTERED_LD_FILENAME),
            gene_df=gene_info_df,
            start_bp=start,
            end_bp=end,
            chrom=chrom,
        )
        out_path = scratch_dir
        write_plots_to_dir(
            out_path,
            {
                "plot": fig
            }
        )
        import pdb; pdb.set_trace()
        return DirectoryAsset(out_path)



def plot_locus_tracks_matplotlib(
    gwas_df: pl.DataFrame,
    susie_cs_df: pl.DataFrame,
    ld_np: np.ndarray,
    gene_df: pl.DataFrame,
    start_bp: int ,
    end_bp: int ,
    chrom: int,
    *,
    gwas_pos_col: str = GWASLAB_POS_COL,
    gwas_mlog10p_col: str = GWASLAB_MLOG10P_COL,
    susie_pos_col: str = GWASLAB_POS_COL,
    susie_pip_col: str = PIP_COLUMN,
    susie_cs_col: str = CS_COLUMN,
    gene_start_col: str = GENE_INFO_START_COL,
    gene_end_col: str = GENE_INFO_END_COL,
    gene_name_col: str = GENE_INFO_NAME_COL,
    gene_strand_col: str = GENE_INFO_STRAND_COL,
    # title: str = "Locus tracks",
        max_mlog10p:float=200,
)-> Figure:
    ld2 = ld_np**2

    gwas_df = gwas_df.with_columns(
        pl.min_horizontal(pl.lit(max_mlog10p), pl.col(GWASLAB_MLOG10P_COL)).alias(GWASLAB_MLOG10P_COL),
    )

    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(
        nrows=4, ncols=2,
        height_ratios=[2.2, 2.2, 1.3, 1.3],
        width_ratios=[1.0, 0.04],
        hspace=0.08,
        wspace=0.05,
    )
    ax_manh = fig.add_subplot(gs[1,0])
    ax_ld = fig.add_subplot(gs[0,0], sharex=ax_manh)
    ax_pip = fig.add_subplot(gs[2,0], sharex=ax_manh)
    ax_gene = fig.add_subplot(gs[3,0], sharex=ax_manh)

    ld_cax = fig.add_subplot(gs[0, 1])
    manh_cax = fig.add_subplot(gs[1, 1])

    #1:  Manhattan
    lead = int(np.argmax(gwas_df[gwas_mlog10p_col]))
    ld2_with_lead = ld2[lead,:]
    sc=ax_manh.scatter(gwas_df[gwas_pos_col].to_numpy(), gwas_df[gwas_mlog10p_col], s=10, linewidths=0,
                    c=ld2_with_lead, cmap="plasma")
    ax_manh.scatter(
        gwas_df[gwas_pos_col][lead],gwas_df[gwas_mlog10p_col][lead],s=35,marker="^",
    )
    ax_manh.set_ylabel(r"$-\log_{10}(p)$")


    cbar = fig.colorbar(sc, cax=manh_cax, shrink=0.8)
    cbar.set_label(r"$r^2$ with lead")


    #2: SUSIE results
    if len(susie_cs_df)>0:
        for cs, sub in susie_cs_df.group_by(susie_cs_col, maintain_order=True):
            ax_pip.vlines(
                sub[susie_pos_col].to_numpy(),
                0.0,
                sub[susie_pip_col].to_numpy(),
                linewidth=0.9,
                label=f"CS {cs}",
            )
        # ax_pip.legend(loc="upper right", fontsize=8, frameon=False, ncols=2)
    ax_pip.set_ylabel("PIP (SUSIE)")
    # ax_manh.set_title(title)

    # #3 ld

    ar, edges= get_array_and_edges_for_ld_heatmap(
        ld2=ld2,
        pos=gwas_df[gwas_pos_col].to_numpy(),
    )

    mesh = ax_ld.pcolormesh(
        edges, edges, ld2,
        shading="auto",
        vmin=0 , vmax=1,cmap="plasma"
    )
    ax_ld.set_ylabel("bp")
    ax_ld.set_ylim(float(edges[0]), float(edges[-1]))

    cbar = fig.colorbar(mesh, cax=ld_cax, shrink=0.8)
    cbar.set_label(r"$r^2$")

    # 4: Gene tracks
    plot_gene_tracks(
        ax=ax_gene,
        gene_df=gene_df,
        start_bp=start_bp,
        end_bp=end_bp,
        gene_start_col=gene_start_col,
        gene_end_col=gene_end_col,
        gene_name_col=gene_name_col,
        gene_strand_col=gene_strand_col,
    )


    # axis config
    ax_gene.set_xlim(start_bp-1, end_bp+1)
    ax_gene.set_xlabel(f"Chromosome {chrom} (bp)")

    for ax in [ax_manh, ax_pip, ax_ld]:
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.set_xticks([])

    for ax in [ax_manh, ax_pip, ax_ld, ax_gene]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if ax == ax_gene:
            ax.spines["left"].set_visible(False)

    return fig


@frozen
class BinOptions:
    num_bins: int

def get_array_and_edges_for_ld_heatmap(
        ld2:np.ndarray,
        pos: np.ndarray,
        bin_options: BinOptions|None=None,
) -> tuple[np.ndarray, np.ndarray]:
    da = xr.DataArray(
        ld2,
        coords={
            "x":pos,
            "y":pos,
        },
        dims=["x", "y"],
    )
    if bin_options is not None :
        da= da.groupby_bins(
            "x",bins=bin_options.num_bins,
        ).mean().groupby_bins("y",bins=bin_options.num_bins)
    else:
        da = da.groupby("x").mean().groupby("y").mean()

    new_pos = da.coords["x"].to_numpy()
    mid = (new_pos[:-1] + new_pos[1:]) / 2.0
    first = new_pos[0] - (mid[0] - new_pos[0])
    last = new_pos[-1] + (new_pos[-1] - mid[-1])
    edges= np.concatenate([[first], mid, [last]])
    return da.to_numpy(), edges


def plot_gene_tracks(
        ax,
        gene_df: pl.DataFrame,
        start_bp: int,
        end_bp: int,
        gene_start_col: str,
        gene_end_col: str,
        gene_name_col: str,
        gene_strand_col: str,
        font_size: int = 9,
        min_dist_between_genes: float = 0.02,  # Fraction of x-axis
):
    """
    Gene track plotting function generated with help from Gemini
    """
    # 1. Filter genes within the window
    # We allow genes that partially overlap the window
    region_df = gene_df.filter(
        (pl.col(gene_end_col) >= start_bp) &
        (pl.col(gene_start_col) <= end_bp)
    ).sort(gene_start_col)

    if region_df.height == 0:
        ax.text(0.5, 0.5, "No genes in region", transform=ax.transAxes, ha='center')
        ax.set_yticks([])
        return

    # 2. Heuristic for text width in coordinate space (bp)
    # This assumes a standard matplotlib figure width.
    # If the text is overlapping, increase 'char_width_factor'.
    total_bp = end_bp - start_bp
    char_width_factor = 0.012  # approx 1.2% of plot width per character
    bp_per_char = total_bp * char_width_factor
    min_gap = total_bp * min_dist_between_genes

    # 3. Greedy Packing Algorithm
    # lanes is a list of end-coordinates for each y-level.
    # We store the 'occupied_end' which includes the text width.
    lanes = []

    # Store plot data to draw later
    gene_placements = []

    for row in region_df.iter_rows(named=True):
        g_start = max(start_bp, row[gene_start_col])
        g_end = min(end_bp, row[gene_end_col])
        g_name = row[gene_name_col]

        # Calculate the "Visual End" of this object
        # It's the maximum of the gene structure OR the text label length
        text_end = g_start + (len(g_name) * bp_per_char)
        visual_end = max(g_end, text_end) + min_gap

        # Find the first lane where this gene fits
        y_level = -1
        for i, lane_end in enumerate(lanes):
            if lane_end < g_start:
                y_level = i
                lanes[i] = visual_end  # Update lane with new end
                break

        # If no lane fits, create a new one
        if y_level == -1:
            lanes.append(visual_end)
            y_level = len(lanes) - 1

        gene_placements.append({
            "start": row[gene_start_col],  # Use real start for drawing
            "end": row[gene_end_col],
            "name": g_name,
            "strand": row[gene_strand_col],
            "y": y_level
        })

    # 4. Plotting
    # Auto-scale y-axis based on number of lanes
    n_lanes = len(lanes)
    ax.set_ylim(-0.5, n_lanes + 0.5)

    # Optional: If many lanes, shrink font
    if n_lanes > 10:
        font_size = max(6, font_size - 2)

    for g in gene_placements:
        y = g["y"]
        # Draw the main gene line
        ax.plot([g["start"], g["end"]], [y, y], color='navy', lw=1.5)

        # Draw vertical ticks at ends (like exon boundaries)
        tick_h = 0.2
        ax.plot([g["start"], g["start"]], [y - tick_h, y + tick_h], color='navy', lw=1)
        ax.plot([g["end"], g["end"]], [y - tick_h, y + tick_h], color='navy', lw=1)

        # Draw Label (placed slightly above the line, aligned left)
        # Using a white outline (path_effects) helps readability if overlapping crowded lines
        ax.text(
            g["start"], y + 0.3,
            f"  {g['name']}",  # space for padding
            ha='left', va='center',
            fontsize=font_size,
            style='italic',
            path_effects=[pe.withStroke(linewidth=2, foreground="white")]
        )

        # Draw Strand Arrows
        # We place an arrow at the midpoint or end
        mid = (max(start_bp, g["start"]) + min(end_bp, g["end"])) / 2
        marker = '>' if g["strand"] == '+' else '<'
        ax.scatter([mid], [y], marker=marker, color='navy', s=30, zorder=10)

    # Clean up axes
    ax.set_yticks([])
    # ax.set_xlim(start_bp, end_bp)

    # Invert y-axis so first gene is at top (optional, matches standard genome browsers)
    ax.invert_yaxis()