"""
Stacked locus plots of chosen trait x protein colocalizations.

One figure per requested (locus, protein) pair, saved as PNG and SVG, with three panels sharing the
GRCh38 x-axis over the locus window: the trait's -log10 p, the protein's -log10 p (both over the
aligned variants the colocalization used, each with its lead variant highlighted), and the genes
in the window. The title carries the gene symbol, cis/trans status and PP.H4.

The output directory also holds plot_index.parquet (PROTEIN_ID, GENE, LOCUS_ID, PNG_FILENAME,
SVG_FILENAME). Pairs are chosen by hand rather than by a PP.H4 threshold, so only loci worth
a close look get figures; any tested pair can be plotted, including H3-dominant ones.
"""

from pathlib import Path

import numpy as np
import polars as pl
from attrs import frozen
from matplotlib.figure import Figure

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.locus_plot_helpers import (
    plot_association,
    plot_gene_track_panel,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_display_table_task import (
    read_raw_results,
)
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    two_sided_mlog10p,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_POS_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    BETA_PROTEIN_COL,
    BETA_TRAIT_COL,
    CIS_STATUS_COL,
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
    STATUS_COL,
    STATUS_TESTED,
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
class ColocPair:
    """One trait locus x protein pair of a coloc run, as named in its raw results."""

    locus_id: str
    protein_id: str


@frozen(slots=True)
class _PlotSpec:
    """One requested pair and where its aligned window and figure files live."""

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


def _plot_specs(raw: pl.DataFrame, pairs: tuple[ColocPair, ...]) -> list[_PlotSpec]:
    tested = raw.filter(pl.col(STATUS_COL) == STATUS_TESTED)
    specs = []
    for pair in pairs:
        rows = tested.filter(
            pl.col(LOCUS_ID_COL) == pair.locus_id,
            pl.col(PROTEIN_ID_COL) == pair.protein_id,
        )
        assert rows.height == 1, f"{pair} is not a tested pair of the coloc results"
        row = rows.row(0, named=True)
        specs.append(
            _PlotSpec(
                protein_id=row[PROTEIN_ID_COL],
                gene=row[GENE_COL],
                locus_id=row[LOCUS_ID_COL],
                chrom=row[LOCUS_CHR_COL],
                start=row[LOCUS_START_COL],
                end=row[LOCUS_END_COL],
                cis_status=row[CIS_STATUS_COL],
                pp_h4=row[PP_H4_COL],
            )
        )
    return specs


def _window_mlog10p(window: pl.DataFrame, beta: str, se: str) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"):
        return two_sided_mlog10p((window[beta] / window[se]).to_numpy())


def _render(
    spec: _PlotSpec, window: pl.DataFrame, genes: pl.DataFrame, out_dir: Path
) -> None:
    fig = Figure(figsize=(_FIGURE_WIDTH_IN, _PANEL_HEIGHT_IN * 3))
    ax_trait = fig.add_subplot(3, 1, 1)
    ax_protein = fig.add_subplot(3, 1, 2, sharex=ax_trait)
    ax_genes = fig.add_subplot(3, 1, 3, sharex=ax_trait)
    pos = window[GWASLAB_POS_COL].to_numpy()
    plot_association(
        ax_trait,
        pos,
        _window_mlog10p(window, BETA_TRAIT_COL, SE_TRAIT_COL),
        label="trait",
        significance_line=True,
    )
    plot_association(
        ax_protein,
        pos,
        _window_mlog10p(window, BETA_PROTEIN_COL, SE_PROTEIN_COL),
        label=spec.gene,
        significance_line=True,
    )
    plot_gene_track_panel(ax_genes, genes, spec.chrom, spec.start, spec.end)
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
    gene_info_task: a GRCh38 gene table with the GENE_INFO_* columns.
    pairs: the tested pairs to plot, one figure each.
    """

    meta: ResultDirectoryMeta
    coloc_task: Task
    gene_info_task: Task
    pairs: tuple[ColocPair, ...]

    def __attrs_post_init__(self) -> None:
        assert len(self.pairs) > 0, "no pairs to plot"
        assert len(set(self.pairs)) == len(self.pairs), f"duplicate pairs: {self.pairs}"

    @property
    def deps(self) -> list[Task]:
        return [self.coloc_task, self.gene_info_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        coloc_asset = fetch(self.coloc_task.asset_id)
        assert isinstance(coloc_asset, DirectoryAsset)
        specs = _plot_specs(read_raw_results(fetch, self.coloc_task), self.pairs)
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
        gene_info_task: Task,
        pairs: tuple[ColocPair, ...],
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
            gene_info_task=gene_info_task,
            pairs=pairs,
        )
