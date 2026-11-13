"""
Stacked locus plots of the colocalizing proteins at each locus of a coloc display table.

The loci are the distinct loci of the display table (pairs with PP.H4 above its threshold),
ranked by their best PP.H4 (the posterior probability of a shared causal variant), and at most
max_loci are plotted. One figure per locus, saved as PNG and SVG, with rows sharing the GRCh38
x-axis over the locus window:
1. The trait: -log10 p of every usable variant (finite beta, finite SE > 0) of the full trait
   GWAS in the window, with a dashed genome-wide-significance line, and the recombination rate (cM/Mb,
   from a genetic map) on a second y-axis as rough LD context.
2. One row per display-table protein at the locus, in PP.H4 order: -log10 p over the aligned
   variants coloc used (the pair's saved window), and a dotted vertical
   line at the trait lead position so offsets between the peaks are visible.
3. The genes in the window.

The output directory also holds stack_plot_index.parquet (LOCUS_ID, N_PROTEINS, PNG_FILENAME,
SVG_FILENAME). An empty display table gives an empty index and no figures.
"""

from pathlib import Path

import numpy as np
import polars as pl
from attrs import frozen
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.locus_plot_helpers import (
    load_recomb,
    plot_association,
    plot_gene_track_panel,
    plot_recomb,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_display_table_task import (
    coloc_result_dir,
    ranked_display_loci,
)
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    two_sided_mlog10p,
)
from mecfs_bio.build_system.task.protein_coloc.trait_protein_coloc_task import (
    scan_trait,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    BETA_PROTEIN_COL,
    DISPLAY_CIS_STATUS_COL,
    DISPLAY_GENE_COL,
    DISPLAY_LOCUS_COL,
    DISPLAY_PP_H4_COL,
    DISPLAY_PROTEIN_COL,
    LOCUS_CHR_COL,
    LOCUS_END_COL,
    LOCUS_ID_COL,
    LOCUS_START_COL,
    N_PROTEINS_COL,
    PLOT_PNG_FILENAME_COL,
    PLOT_SVG_FILENAME_COL,
    RAW_RESULTS_FILENAME,
    SE_PROTEIN_COL,
    WINDOWS_DIRNAME,
    pair_file_stem,
)

STACK_PLOT_INDEX_FILENAME = "stack_plot_index.parquet"
_FIGURE_WIDTH_IN = 12.0
_PANEL_HEIGHT_IN = 2.2
_PNG_DPI = 150

STACK_PLOT_INDEX_SCHEMA = pl.Schema(
    {
        LOCUS_ID_COL: pl.String,
        N_PROTEINS_COL: pl.Int32,
        PLOT_PNG_FILENAME_COL: pl.String,
        PLOT_SVG_FILENAME_COL: pl.String,
    }
)


@frozen(slots=True)
class _ProteinRow:
    gene: str
    protein_id: str
    cis_status: str
    pp_h4: float


@frozen(slots=True)
class _LocusSpec:
    locus_id: str
    chrom: int
    start: int
    end: int
    proteins: tuple[_ProteinRow, ...]

    @property
    def stem(self) -> str:
        return self.locus_id.replace(":", "_")


@frozen(slots=True)
class _TrackData:
    """Positions and -log10 p of one Manhattan row."""

    pos: np.ndarray
    mlog10p: np.ndarray


def _mlog10p(beta: pl.Series, se: pl.Series) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"):
        return two_sided_mlog10p((beta / se).to_numpy())


def _read_table(fetch: Fetch, task: Task) -> pl.DataFrame:
    return (
        scan_dataframe_asset(fetch(task.asset_id), meta=task.meta).collect().to_polars()
    )


def _locus_specs(
    display: pl.DataFrame, raw: pl.DataFrame, locus_ids: list[str]
) -> list[_LocusSpec]:
    bounds = raw.select(
        LOCUS_ID_COL, LOCUS_CHR_COL, LOCUS_START_COL, LOCUS_END_COL
    ).unique(LOCUS_ID_COL)
    specs = []
    for locus_id in locus_ids:
        bound = bounds.filter(pl.col(LOCUS_ID_COL) == locus_id)
        assert bound.height == 1, f"locus {locus_id} is not in the coloc raw results"
        rows = display.filter(pl.col(DISPLAY_LOCUS_COL) == locus_id).sort(
            DISPLAY_PP_H4_COL, descending=True
        )
        specs.append(
            _LocusSpec(
                locus_id=locus_id,
                chrom=int(bound[LOCUS_CHR_COL].item()),
                start=int(bound[LOCUS_START_COL].item()),
                end=int(bound[LOCUS_END_COL].item()),
                proteins=tuple(
                    _ProteinRow(
                        gene=row[DISPLAY_GENE_COL],
                        protein_id=row[DISPLAY_PROTEIN_COL],
                        cis_status=row[DISPLAY_CIS_STATUS_COL],
                        pp_h4=row[DISPLAY_PP_H4_COL],
                    )
                    for row in rows.iter_rows(named=True)
                ),
            )
        )
    return specs


def _trait_track(trait: pl.LazyFrame, spec: _LocusSpec) -> _TrackData:
    window = (
        trait.filter(
            pl.col(GWASLAB_CHROM_COL) == spec.chrom,
            pl.col(GWASLAB_POS_COL).is_between(spec.start, spec.end),
            pl.col(GWASLAB_BETA_COL).is_finite(),
            pl.col(GWASLAB_SE_COL).is_finite() & (pl.col(GWASLAB_SE_COL) > 0),
        )
        .select(GWASLAB_POS_COL, GWASLAB_BETA_COL, GWASLAB_SE_COL)
        .collect()
        .sort(GWASLAB_POS_COL)
    )
    assert window.height > 0, f"no trait variants in {spec.locus_id}"
    return _TrackData(
        pos=window[GWASLAB_POS_COL].to_numpy(),
        mlog10p=_mlog10p(window[GWASLAB_BETA_COL], window[GWASLAB_SE_COL]),
    )


def _protein_track(windows_dir: Path, locus_id: str, protein_id: str) -> _TrackData:
    window = pl.read_parquet(
        windows_dir / f"{pair_file_stem(locus_id, protein_id)}.parquet"
    ).sort(GWASLAB_POS_COL)
    return _TrackData(
        pos=window[GWASLAB_POS_COL].to_numpy(),
        mlog10p=_mlog10p(window[BETA_PROTEIN_COL], window[SE_PROTEIN_COL]),
    )


def _strip_spines(ax: Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _render_locus(
    spec: _LocusSpec,
    trait_track: _TrackData,
    protein_tracks: list[_TrackData],
    recomb: pl.DataFrame,
    genes: pl.DataFrame,
    out_dir: Path,
) -> None:
    n_rows = 2 + len(protein_tracks)
    fig = Figure(figsize=(_FIGURE_WIDTH_IN, _PANEL_HEIGHT_IN * n_rows))
    ax_trait = fig.add_subplot(n_rows, 1, 1)
    lead = plot_association(
        ax_trait,
        trait_track.pos,
        trait_track.mlog10p,
        label="trait",
        significance_line=True,
    )
    trait_lead_pos = float(trait_track.pos[lead])
    plot_recomb(ax_trait.twinx(), recomb)
    ax_trait.set_xlim(spec.start, spec.end)
    _strip_spines(ax_trait)
    for i, (protein, track) in enumerate(
        zip(spec.proteins, protein_tracks, strict=True)
    ):
        ax = fig.add_subplot(n_rows, 1, 2 + i, sharex=ax_trait)
        plot_association(
            ax, track.pos, track.mlog10p, label=protein.gene, significance_line=False
        )
        ax.axvline(trait_lead_pos, linestyle=":", color="grey", linewidth=1)
        ax.set_title(
            f"{protein.gene} ({protein.cis_status}) PP.H4={protein.pp_h4:.2f}",
            fontsize=9,
            loc="left",
        )
        _strip_spines(ax)
    ax_genes = fig.add_subplot(n_rows, 1, n_rows, sharex=ax_trait)
    plot_gene_track_panel(ax_genes, genes, spec.chrom, spec.start, spec.end)
    _strip_spines(ax_genes)
    for ax in fig.axes:
        if ax is not ax_genes:
            ax.tick_params(axis="x", labelbottom=False)
    fig.suptitle(spec.locus_id)
    fig.savefig(out_dir / f"{spec.stem}.png", dpi=_PNG_DPI, bbox_inches="tight")
    fig.savefig(out_dir / f"{spec.stem}.svg", bbox_inches="tight")


def _plot_index(specs: list[_LocusSpec]) -> pl.DataFrame:
    return pl.DataFrame(
        [
            {
                LOCUS_ID_COL: spec.locus_id,
                N_PROTEINS_COL: len(spec.proteins),
                PLOT_PNG_FILENAME_COL: f"{spec.stem}.png",
                PLOT_SVG_FILENAME_COL: f"{spec.stem}.svg",
            }
            for spec in specs
        ],
        schema=STACK_PLOT_INDEX_SCHEMA,
    )


@frozen(slots=True)
class ColocLocusStackPlotTask(Task):
    """
    coloc_task: the TraitProteinColocTask whose saved windows give the protein rows.
    trait_task: that coloc task's trait GWAS (GRCh38), plotted in full over each locus window.
    display_task: the ColocDisplayTableTask of the same coloc run; it chooses the loci and
        proteins.
    gene_info_task: a GRCh38 gene table with the GENE_INFO_* columns.
    genetic_map_task: a GRCh38 genetic map with the GMAP_* columns.
    max_loci: at most this many loci are plotted, best PP.H4 first.
    """

    meta: ResultDirectoryMeta
    coloc_task: Task
    trait_task: Task
    display_task: Task
    gene_info_task: Task
    genetic_map_task: Task
    max_loci: int

    def __attrs_post_init__(self) -> None:
        assert self.max_loci >= 1, f"max_loci must be at least 1: {self.max_loci}"

    @property
    def deps(self) -> list[Task]:
        return [
            self.coloc_task,
            self.trait_task,
            self.display_task,
            self.gene_info_task,
            self.genetic_map_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        coloc_dir = coloc_result_dir(fetch, self.coloc_task)
        display = _read_table(fetch, self.display_task)
        specs = _locus_specs(
            display,
            pl.read_parquet(coloc_dir / RAW_RESULTS_FILENAME),
            ranked_display_loci(display, self.max_loci),
        )
        out_dir = scratch_dir / "stack_plots"
        out_dir.mkdir()
        _render_all(self, specs, fetch, coloc_dir / WINDOWS_DIRNAME, out_dir)
        _plot_index(specs).write_parquet(out_dir / STACK_PLOT_INDEX_FILENAME)
        return DirectoryAsset(out_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        coloc_task: Task,
        trait_task: Task,
        display_task: Task,
        gene_info_task: Task,
        genetic_map_task: Task,
        max_loci: int = 10,
    ) -> "ColocLocusStackPlotTask":
        source_meta = coloc_task.meta
        assert isinstance(source_meta, ResultDirectoryMeta), (
            f"coloc_task must carry ResultDirectoryMeta, got {type(source_meta).__name__}"
        )
        return cls(
            meta=ResultDirectoryMeta(
                id=asset_id, trait=source_meta.trait, project=source_meta.project
            ),
            coloc_task=coloc_task,
            trait_task=trait_task,
            display_task=display_task,
            gene_info_task=gene_info_task,
            genetic_map_task=genetic_map_task,
            max_loci=max_loci,
        )


def _render_all(
    task: ColocLocusStackPlotTask,
    specs: list[_LocusSpec],
    fetch: Fetch,
    windows_dir: Path,
    out_dir: Path,
) -> None:
    if not specs:
        return
    trait = scan_trait(task.trait_task, fetch)
    genes = _read_table(fetch, task.gene_info_task)
    for spec in specs:
        _render_locus(
            spec,
            trait_track=_trait_track(trait, spec),
            protein_tracks=[
                _protein_track(windows_dir, spec.locus_id, protein.protein_id)
                for protein in spec.proteins
            ],
            recomb=load_recomb(
                fetch, task.genetic_map_task, spec.chrom, spec.start, spec.end
            ),
            genes=genes,
            out_dir=out_dir,
        )
