"""
Stacked locus plots of the displayed trait x protein colocalizations.

One figure per row of the display table, saved as PNG and SVG, with three panels sharing the
GRCh38 x-axis over the locus window: the trait's -log10 p, the protein's -log10 p (both over the
aligned variants the colocalization used, each with its lead variant highlighted), and the genes
in the window. The title carries the gene symbol, cis/trans status and PP.H4.

The output directory also holds plot_index.parquet (PROTEIN_ID, GENE, LOCUS_ID, PNG_FILENAME,
SVG_FILENAME), which is empty when no pair passed the display threshold.
"""

from pathlib import Path

import numpy as np
import polars as pl
from attrs import frozen
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.protein_coloc.coloc_display_table_task import (
    ColocDisplayTableTask,
    read_raw_results,
)
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    two_sided_mlog10p,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
    GENE_INFO_STRAND_COL,
    GWAS_SIGNIFICANCE_MLOG10P,
    plot_gene_tracks,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_NAME_FOR_CODE,
    GWASLAB_POS_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    BETA_PROTEIN_COL,
    BETA_TRAIT_COL,
    CIS_STATUS_COL,
    DISPLAY_LOCUS_COL,
    GENE_COL,
    LOCUS_CHR_COL,
    LOCUS_END_COL,
    LOCUS_ID_COL,
    LOCUS_START_COL,
    PLOT_PNG_FILENAME_COL,
    PLOT_SVG_FILENAME_COL,
    PP_H4_COL,
    PROTEIN_ID_COL,
    SE_PROTEIN_COL,
    SE_TRAIT_COL,
    WINDOWS_DIRNAME,
    pair_file_stem,
)

PLOT_INDEX_FILENAME = "plot_index.parquet"
_FIGURE_WIDTH_IN = 12.0
_PANEL_HEIGHT_IN = 2.6
_PNG_DPI = 150

PLOT_INDEX_SCHEMA = pl.Schema(
    {
        PROTEIN_ID_COL: pl.String,
        GENE_COL: pl.String,
        LOCUS_ID_COL: pl.String,
        PLOT_PNG_FILENAME_COL: pl.String,
        PLOT_SVG_FILENAME_COL: pl.String,
    }
)


@frozen(slots=True)
class _PlotSpec:
    """One displayed pair and where its aligned window and figure files live."""

    protein_id: str
    gene: str
    locus_id: str
    chrom: int
    start: int
    end: int
    cis_status: str
    pp_h4: float

    @property
    def stem(self) -> str:
        return pair_file_stem(self.locus_id, self.protein_id)


def _plot_specs(display: pl.DataFrame, raw: pl.DataFrame) -> list[_PlotSpec]:
    loci = raw.select(
        LOCUS_ID_COL, LOCUS_CHR_COL, LOCUS_START_COL, LOCUS_END_COL
    ).unique(LOCUS_ID_COL)
    joined = display.join(loci, left_on=DISPLAY_LOCUS_COL, right_on=LOCUS_ID_COL)
    return [
        _PlotSpec(
            protein_id=row[PROTEIN_ID_COL],
            gene=row[GENE_COL],
            locus_id=row[DISPLAY_LOCUS_COL],
            chrom=row[LOCUS_CHR_COL],
            start=row[LOCUS_START_COL],
            end=row[LOCUS_END_COL],
            cis_status=row[CIS_STATUS_COL],
            pp_h4=row[PP_H4_COL],
        )
        for row in joined.iter_rows(named=True)
    ]


def _plot_association(
    ax: Axes, window: pl.DataFrame, beta: str, se: str, label: str
) -> None:
    pos = window[GWASLAB_POS_COL].to_numpy()
    with np.errstate(invalid="ignore", divide="ignore"):
        mlog10p = two_sided_mlog10p((window[beta] / window[se]).to_numpy())
    ax.scatter(pos, mlog10p, s=8, color="tab:blue").set_rasterized(True)
    lead = int(np.nanargmax(mlog10p))
    ax.scatter(pos[lead], mlog10p[lead], s=40, marker="^", color="black")
    ax.axhline(y=GWAS_SIGNIFICANCE_MLOG10P, color="grey", linestyle="--", linewidth=1)
    ax.set_ylabel(f"{label}\n-log10 p")


def _plot_genes(ax: Axes, genes: pl.DataFrame, spec: _PlotSpec) -> None:
    chrom_name = GWASLAB_CHROM_NAME_FOR_CODE.get(spec.chrom, None) or str(spec.chrom)
    plot_gene_tracks(
        ax=ax,
        gene_df=genes.filter(pl.col(GENE_INFO_CHROM_COL).cast(pl.String) == chrom_name),
        start_bp=spec.start,
        end_bp=spec.end,
        gene_start_col=GENE_INFO_START_COL,
        gene_end_col=GENE_INFO_END_COL,
        gene_name_col=GENE_INFO_NAME_COL,
        gene_strand_col=GENE_INFO_STRAND_COL,
    )
    ax.set_ylabel("genes")
    ax.set_xlabel(f"GRCh38 chr{chrom_name} position (bp)")


def _render(
    spec: _PlotSpec, window: pl.DataFrame, genes: pl.DataFrame, out_dir: Path
) -> None:
    fig = Figure(figsize=(_FIGURE_WIDTH_IN, _PANEL_HEIGHT_IN * 3))
    ax_trait = fig.add_subplot(3, 1, 1)
    ax_protein = fig.add_subplot(3, 1, 2, sharex=ax_trait)
    ax_genes = fig.add_subplot(3, 1, 3, sharex=ax_trait)
    _plot_association(ax_trait, window, BETA_TRAIT_COL, SE_TRAIT_COL, "trait")
    _plot_association(ax_protein, window, BETA_PROTEIN_COL, SE_PROTEIN_COL, spec.gene)
    _plot_genes(ax_genes, genes, spec)
    ax_trait.set_xlim(spec.start, spec.end)
    for ax in (ax_trait, ax_protein):
        ax.tick_params(axis="x", labelbottom=False)
    for ax in (ax_trait, ax_protein, ax_genes):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(f"{spec.gene} ({spec.cis_status}) PP.H4={spec.pp_h4:.2f}")
    fig.savefig(out_dir / f"{spec.stem}.png", dpi=_PNG_DPI, bbox_inches="tight")
    fig.savefig(out_dir / f"{spec.stem}.svg", bbox_inches="tight")


def _plot_index(specs: list[_PlotSpec]) -> pl.DataFrame:
    return pl.DataFrame(
        [
            {
                PROTEIN_ID_COL: spec.protein_id,
                GENE_COL: spec.gene,
                LOCUS_ID_COL: spec.locus_id,
                PLOT_PNG_FILENAME_COL: f"{spec.stem}.png",
                PLOT_SVG_FILENAME_COL: f"{spec.stem}.svg",
            }
            for spec in specs
        ],
        schema=PLOT_INDEX_SCHEMA,
    )


def _render_all(
    specs: list[_PlotSpec], coloc_dir: Path, genes: pl.DataFrame, out_dir: Path
) -> None:
    for spec in specs:
        window = pl.read_parquet(coloc_dir / WINDOWS_DIRNAME / f"{spec.stem}.parquet")
        _render(spec, window, genes, out_dir)


@frozen(slots=True)
class ColocLocusPlotTask(Task):
    """
    coloc_task: the TraitProteinColocTask whose aligned windows are plotted.
    display_task: selects which pairs to plot.
    gene_info_task: a GRCh38 gene table with the GENE_INFO_* columns.
    """

    meta: ResultDirectoryMeta
    coloc_task: Task
    display_task: ColocDisplayTableTask
    gene_info_task: Task

    @property
    def deps(self) -> list[Task]:
        return [self.coloc_task, self.display_task, self.gene_info_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        coloc_asset = fetch(self.coloc_task.asset_id)
        assert isinstance(coloc_asset, DirectoryAsset)
        display_asset = fetch(self.display_task.asset_id)
        assert isinstance(display_asset, FileAsset)
        specs = _plot_specs(
            pl.read_parquet(display_asset.path),
            read_raw_results(fetch, self.coloc_task),
        )
        genes = (
            scan_dataframe_asset(
                fetch(self.gene_info_task.asset_id), meta=self.gene_info_task.meta
            )
            .collect()
            .to_polars()
        )
        out_dir = scratch_dir / "plots"
        out_dir.mkdir()
        _render_all(specs, coloc_asset.path, genes, out_dir)
        _plot_index(specs).write_parquet(out_dir / PLOT_INDEX_FILENAME)
        return DirectoryAsset(out_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        coloc_task: Task,
        display_task: ColocDisplayTableTask,
        gene_info_task: Task,
    ) -> "ColocLocusPlotTask":
        source_meta = coloc_task.meta
        assert isinstance(source_meta, ResultDirectoryMeta), (
            f"coloc_task must carry ResultDirectoryMeta, got {type(source_meta).__name__}"
        )
        return cls(
            meta=ResultDirectoryMeta(
                id=asset_id, trait=source_meta.trait, project=source_meta.project
            ),
            coloc_task=coloc_task,
            display_task=display_task,
            gene_info_task=gene_info_task,
        )
