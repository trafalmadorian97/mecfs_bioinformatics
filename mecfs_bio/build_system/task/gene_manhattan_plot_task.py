"""
Task to produce an interactive gene-level Manhattan plot as an HTML file.

Supports two source types:

- :class:`MagmaGeneSource`: read a MAGMA gene-level analysis output directory
  (the ``.genes.out`` file produced by :class:`MagmaGeneAnalysisTask`) and join
  a gene thesaurus to translate Ensembl IDs into human-readable gene names.
- :class:`GenePValueTableSource`: read an arbitrary table of
  ``(gene_ensembl_id, p_value)`` rows and look up chromosomal locations and
  human-readable gene names from a gene-locations reference (such as
  ``MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW``). Intended for
  rare-variant test output or any other gene-level result table.

The plot uses Plotly's WebGL ``Scattergl`` renderer for performance with
20k-30k gene points and exposes hover text containing the gene name, Ensembl
ID, chromosome, genomic midpoint position (labelled ``Position (hg19)`` or
``Position (hg38)`` according to the source's declared ``genome_build``), and
``-log10(p)``.
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.plot_file_meta import GWASPlotFileMeta
from mecfs_bio.build_system.meta.processed_gwas_data_directory_meta import (
    ProcessedGwasDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameWhiteSpaceSepTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (
    scan_dataframe,
    scan_dataframe_asset,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.magma.magma_gene_analysis_task import (
    GENE_ANALYSIS_OUTPUT_STEM_NAME,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.genomic_coordinate_constants import (
    GenomeBuild,
    extended_mhc_interval,
)
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval
from mecfs_bio.util.plotting.save_fig import PlotlyWriteMode

logger = structlog.get_logger()

# Internal canonical columns used by the plotting routine.
_CHROM = "chrom"
_POS = "pos"
_ENSEMBL_ID = "ensembl_id"
_GENE_NAME = "gene_name"
_P = "p_value"

# Ensembl Biomart gene-thesaurus columns (see
# mecfs_bio/assets/reference_data/ensembl_biomart/gene_thesaurus.py).
_THESAURUS_ENSEMBL_COL = "Gene stable ID"
_THESAURUS_NAME_COL = "Gene name"

# MAGMA ``.genes.out`` columns.
_MAGMA_GENE_COL = "GENE"
_MAGMA_CHR_COL = "CHR"
_MAGMA_START_COL = "START"
_MAGMA_STOP_COL = "STOP"
_MAGMA_P_COL = "P"

# Column name used in the magma gene-locations reference file for the Ensembl
# ID (see magma_ensembl_gene_location_reference_data_build_37.py).
_MAGMA_LOC_ENSEMBL_COL = "ensembl_name"

# Conventional chromosome sort order. Anything else is appended at the end in
# string-sorted order so the plot still renders if exotic contigs sneak in.
_CANONICAL_CHROM_ORDER: list[str] = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]

GeneIdKind = Literal["ensembl_id", "gene_name"]


@frozen
class GeneManhattanData:
    """The genes to plot plus the multiple-testing count for the significance line.

    df holds the rows to plot, after any max_p_value filtering.

    num_genes_for_correction is the number of genes with a valid (positive,
    non-null) p-value before max_p_value filtering. It drives the default
    Bonferroni threshold so that the significance line stays invariant to the
    purely visual max_p_value filter.
    """

    df: pd.DataFrame
    num_genes_for_correction: int


class GeneManhattanSource(ABC):
    """A source that yields rows of (chrom, pos, ensembl_id, gene_name, p) for a Manhattan plot."""

    @property
    @abstractmethod
    def deps(self) -> list[Task]:
        pass

    @property
    @abstractmethod
    def trait(self) -> str:
        """The trait label inherited from the primary input task's metadata."""
        pass

    @property
    @abstractmethod
    def project(self) -> str:
        """The project label inherited from the primary input task's metadata."""
        pass

    @property
    @abstractmethod
    def genome_build(self) -> GenomeBuild:
        """Genome build of the chromosomal positions exposed by load_df.

        Drives the hover-text position label (pos_hg19 vs pos_hg38).
        """
        pass

    @property
    @abstractmethod
    def max_p_value(self) -> float | None:
        """Drop genes whose p-value is at or above this before plotting.

        None disables filtering. Filtering is purely a visual simplification:
        it does not affect the Bonferroni significance threshold, which is based
        on the gene count before filtering.
        """
        pass

    @abstractmethod
    def _load_full_df(self, fetch: Fetch) -> pd.DataFrame:
        """Materialize all candidate rows with columns chrom, pos, ensembl_id, gene_name, p_value.

        Returns every gene the source knows about, with no max_p_value filtering
        applied; the base load_df applies that filter and counts genes for the
        multiple-testing correction.
        """
        pass

    def load_df(self, fetch: Fetch) -> GeneManhattanData:
        """Load the full table, then apply the optional max_p_value filter.

        The multiple-testing count is taken before filtering so that the
        significance threshold is unaffected by max_p_value.
        """
        df = self._load_full_df(fetch)
        valid_p = df[_P].notna() & (df[_P] > 0)
        num_genes_for_correction = int(valid_p.sum())
        if self.max_p_value is not None:
            num_before = len(df)
            df = df[df[_P] < self.max_p_value]
            logger.info(
                "Filtered genes by maximum p-value",
                max_p_value=self.max_p_value,
                num_dropped=num_before - len(df),
                num_kept=len(df),
            )
        return GeneManhattanData(
            df=df, num_genes_for_correction=num_genes_for_correction
        )


@frozen
class MagmaGeneSource(GeneManhattanSource):
    """Load a Manhattan-plot table from a :class:`MagmaGeneAnalysisTask`.

    Chromosomal positions come from the MAGMA output itself. Human-readable
    gene names are joined in from ``gene_thesaurus_task`` by Ensembl ID. When
    a gene is missing from the thesaurus, the Ensembl ID is used as the
    display name.

    max_p_value, when not None, drops genes whose p-value is at or above it
    before plotting, keeping the figure free of the many uninformative
    high-p-value points. The default of 0.1 retains only the nominally
    interesting tail. Filtering does not affect the significance threshold.
    """

    magma_task: Task
    gene_thesaurus_task: Task
    genome_build: GenomeBuild
    max_p_value: float | None = 0.01

    @property
    def deps(self) -> list[Task]:
        return [self.magma_task, self.gene_thesaurus_task]

    @property
    def _magma_meta(self) -> ProcessedGwasDataDirectoryMeta:
        meta = self.magma_task.meta
        assert isinstance(meta, ProcessedGwasDataDirectoryMeta)
        return meta

    @property
    def trait(self) -> str:
        return self._magma_meta.trait

    @property
    def project(self) -> str:
        return self._magma_meta.project

    def _load_full_df(self, fetch: Fetch) -> pd.DataFrame:
        magma_asset = fetch(self.magma_task.asset_id)
        assert isinstance(magma_asset, DirectoryAsset)
        magma_path = magma_asset.path / (GENE_ANALYSIS_OUTPUT_STEM_NAME + ".genes.out")
        magma_df = (
            scan_dataframe(
                path=magma_path,
                spec=DataFrameReadSpec(
                    DataFrameWhiteSpaceSepTextFormat(comment_code="#")
                ),
            )
            .collect()
            .to_pandas()
        )
        magma_df = pd.DataFrame(
            {
                _ENSEMBL_ID: magma_df[_MAGMA_GENE_COL].astype(str),
                _CHROM: magma_df[_MAGMA_CHR_COL].astype(str),
                _POS: (
                    magma_df[_MAGMA_START_COL].astype(float)
                    + magma_df[_MAGMA_STOP_COL].astype(float)
                )
                / 2.0,
                _P: magma_df[_MAGMA_P_COL].astype(float),
            }
        )

        thesaurus_asset = fetch(self.gene_thesaurus_task.asset_id)
        thesaurus_df = (
            scan_dataframe_asset(thesaurus_asset, meta=self.gene_thesaurus_task.meta)
            .collect()
            .to_pandas()
        )
        thesaurus_df = (
            thesaurus_df[[_THESAURUS_ENSEMBL_COL, _THESAURUS_NAME_COL]]
            .rename(
                columns={
                    _THESAURUS_ENSEMBL_COL: _ENSEMBL_ID,
                    _THESAURUS_NAME_COL: _GENE_NAME,
                }
            )
            .drop_duplicates(subset=[_ENSEMBL_ID])
        )
        merged = magma_df.merge(thesaurus_df, on=_ENSEMBL_ID, how="left")
        merged[_GENE_NAME] = merged[_GENE_NAME].fillna(merged[_ENSEMBL_ID])
        return merged


@frozen
class GenePValueTableSource(GeneManhattanSource):
    """Load a Manhattan-plot table from an arbitrary (gene, p-value) table.

    Chromosomal positions and the complementary gene identifier (Ensembl ID or
    human-readable gene name) are looked up from ``gene_locations_task``
    (e.g. the MAGMA Ensembl gene-locations reference) by inner join. Genes
    missing from the locations file are dropped because they cannot be placed
    on the x-axis.

    ``gene_id_kind`` declares which identifier the input table uses in
    ``gene_col``. The locations reference must contain a matching column:
    Ensembl IDs (``"ensembl_id"``) join on the reference's Ensembl-ID column,
    gene symbols (``"gene_name"``) join on the reference's gene-name column.

    max_p_value, when not None, drops genes whose p-value is at or above it
    before plotting, keeping the figure free of the many uninformative
    high-p-value points. The default of 0.1 retains only the nominally
    interesting tail. Filtering does not affect the significance threshold.
    """

    table_task: Task
    gene_locations_task: Task
    gene_col: str
    p_col: str
    genome_build: GenomeBuild
    gene_id_kind: GeneIdKind = "ensembl_id"
    max_p_value: float | None = 0.1

    @property
    def deps(self) -> list[Task]:
        return [self.table_task, self.gene_locations_task]

    @property
    def _table_meta(self) -> ResultTableMeta:
        meta = self.table_task.meta
        assert isinstance(meta, ResultTableMeta)
        return meta

    @property
    def trait(self) -> str:
        return self._table_meta.trait

    @property
    def project(self) -> str:
        return self._table_meta.project

    def _load_full_df(self, fetch: Fetch) -> pd.DataFrame:
        join_col = _ENSEMBL_ID if self.gene_id_kind == "ensembl_id" else _GENE_NAME

        table_asset = fetch(self.table_task.asset_id)
        p_df = (
            scan_dataframe_asset(table_asset, meta=self.table_task.meta)
            .collect()
            .to_pandas()
        )
        p_df = pd.DataFrame(
            {
                join_col: p_df[self.gene_col].astype(str),
                _P: p_df[self.p_col].astype(float),
            }
        )

        loc_asset = fetch(self.gene_locations_task.asset_id)
        loc_df = (
            scan_dataframe_asset(loc_asset, meta=self.gene_locations_task.meta)
            .collect()
            .to_pandas()
        )
        loc_df = pd.DataFrame(
            {
                _ENSEMBL_ID: loc_df[_MAGMA_LOC_ENSEMBL_COL].astype(str),
                _CHROM: loc_df[GENE_INFO_CHROM_COL].astype(str),
                _POS: (
                    loc_df[GENE_INFO_START_COL].astype(float)
                    + loc_df[GENE_INFO_END_COL].astype(float)
                )
                / 2.0,
                _GENE_NAME: loc_df[GENE_INFO_NAME_COL].astype(str),
            }
        ).drop_duplicates(subset=[join_col])

        merged = p_df.merge(loc_df, on=join_col, how="inner")
        missing = len(p_df) - len(merged)
        if missing > 0:
            logger.warning(
                "Dropped genes missing from the locations reference",
                num_dropped=missing,
                num_kept=len(merged),
                gene_id_kind=self.gene_id_kind,
            )
        return merged


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    if chrom in _CANONICAL_CHROM_ORDER:
        return (_CANONICAL_CHROM_ORDER.index(chrom), "")
    return (len(_CANONICAL_CHROM_ORDER), chrom)


def build_manhattan_plot(
    df: pd.DataFrame,
    sig_threshold: float | None,
    point_size: int,
    colors: tuple[str, str],
    sig_line_color: str,
    title: str | None,
    genome_build: GenomeBuild,
    num_genes_for_correction: int | None = None,
    y_axis_start: float | None = None,
    hla_interval: GenomicInterval | None = None,
    hla_marker_symbol: str | None = "diamond",
) -> go.Figure:
    """Construct a Plotly figure containing a gene-level Manhattan plot.

    Genes with non-positive or null p-values are dropped (-log10 is undefined).
    If sig_threshold is None, a Bonferroni-corrected threshold
    0.05 / num_genes_for_correction is used and a dashed horizontal line is
    drawn at the corresponding -log10(p). num_genes_for_correction should be the
    number of genes tested, counted before any p-value filtering of df, so that
    the threshold is invariant to such filtering; it falls back to the number of
    plotted rows when not supplied.

    genome_build selects the hover label for the gene's midpoint position
    (Position (hg19) for build 37, Position (hg38) for build 38). Positions in
    df are assumed to already be in the declared build.

    y_axis_start, when not None, anchors the lower bound of the -log10(p) axis
    (drawn vertically). It is intended to be -log10(max_p_value) so that a
    p-value-filtered plot uses its full vertical extent instead of leaving empty
    space below the lowest surviving point. A small padding is subtracted so
    points sitting right at the cutoff are not sliced by the x-axis; the upper
    bound stays data-driven.

    hla_interval, when not None, marks genes falling inside it (matched on
    chromosome and midpoint position) with hla_marker_symbol instead of the
    default circle, so that extended-HLA/MHC-region genes stand out. Those genes
    keep their chromosome's color; only the symbol changes. hla_marker_symbol
    must be a valid Plotly symbol whenever hla_interval is given.
    """
    assert hla_interval is None or hla_marker_symbol is not None, (
        "hla_marker_symbol is required when hla_interval is given"
    )
    df = df.dropna(subset=[_P]).copy()
    df = df[df[_P] > 0]
    df[_CHROM] = df[_CHROM].astype(str)
    if len(df) == 0:
        raise ValueError(
            "No plottable rows: all gene p-values were null or non-positive."
        )

    chroms = sorted(df[_CHROM].unique().tolist(), key=_chrom_sort_key)

    chrom_max_series = df.groupby(_CHROM)[_POS].max().astype(float)
    chrom_max_pos: dict[str, float] = {
        str(chrom): float(value) for chrom, value in chrom_max_series.items()
    }
    chrom_offsets: dict[str, float] = {}
    chrom_centers: dict[str, float] = {}
    running_offset = 0.0
    for chrom in chroms:
        chrom_offsets[chrom] = running_offset
        chrom_centers[chrom] = running_offset + chrom_max_pos[chrom] / 2.0
        running_offset += chrom_max_pos[chrom]

    df = df.assign(
        _x=df[_POS] + df[_CHROM].map(chrom_offsets),
        _mlog10p=-np.log10(df[_P]),
    )

    if sig_threshold is None:
        n_correction = (
            num_genes_for_correction
            if num_genes_for_correction is not None
            else len(df)
        )
        sig_threshold = 0.05 / n_correction
    sig_y = float(-np.log10(sig_threshold))

    pos_label = f"position (hg{genome_build})"
    fig = go.Figure()
    for idx, chrom in enumerate(chroms):
        chrom_df = df[df[_CHROM] == chrom]
        color = colors[idx % 2]
        if hla_interval is not None and chrom == str(hla_interval.chrom):
            assert hla_marker_symbol is not None  # guaranteed by the top assert
            in_hla = chrom_df[_POS].between(
                hla_interval.start, hla_interval.end, inclusive="left"
            )
            symbol: str | list[str] = np.where(
                in_hla.to_numpy(), hla_marker_symbol, "circle"
            ).tolist()
        else:
            symbol = "circle"
        fig.add_trace(
            go.Scattergl(
                x=chrom_df["_x"],
                y=chrom_df["_mlog10p"],
                mode="markers",
                marker=dict(size=point_size, color=color, symbol=symbol),
                name=f"chr{chrom}",
                customdata=list(
                    zip(
                        chrom_df[_GENE_NAME].astype(str),
                        chrom_df[_ENSEMBL_ID].astype(str),
                        chrom_df[_CHROM].astype(str),
                        chrom_df[_POS].astype(float),
                        strict=True,
                    )
                ),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    f"{pos_label}:" + " chr%{customdata[2]} %{customdata[3]:,.0f}<br>"
                    # "Ensembl: %{customdata[1]}<br>"
                    # "Chromosome: %{customdata[2]}<br>"
                    # f"{pos_label}: " + "%{customdata[3]:,.0f}<br>"
                    "-log<sub>10</sub>(p): %{y:.3f}<br>"
                    "<extra></extra>"
                ),
                showlegend=False,
            )
        )

    fig.add_hline(
        y=sig_y,
        line=dict(color=sig_line_color, dash="dash"),
        annotation_text=f"p = {sig_threshold:.2e}",
        annotation_position="top left",
    )

    yaxis: dict[str, object] = dict(title="-log<sub>10</sub>(p)", zeroline=False)
    if y_axis_start is not None:
        y_top = max(float(df["_mlog10p"].max()), sig_y)
        pad = 0.05 * max(y_top - y_axis_start, 1.0)
        # Drop the lower bound a little below y_axis_start so points sitting right
        # at the cutoff are not sliced in half by the x-axis.
        yaxis["range"] = [y_axis_start - pad, y_top + pad]

    fig.update_layout(
        title=title,
        xaxis=dict(
            tickmode="array",
            tickvals=[chrom_centers[c] for c in chroms],
            ticktext=chroms,
            title="Chromosome",
            showgrid=False,
            zeroline=False,
        ),
        yaxis=yaxis,
        plot_bgcolor="white",
        hovermode="closest",
        showlegend=False,
    )
    return fig


@frozen
class GeneManhattanPlotTask(Task):
    """Create an interactive HTML gene-level Manhattan plot.

    Backed by Plotly's WebGL renderer (Scattergl) so that hover stays
    responsive at gene-scale point counts (~20k-30k).

    When the source declares a max_p_value, the -log10(p) axis starts at
    -log10(max_p_value) rather than 0, so the filtered plot uses its full
    vertical extent.

    hla_marker_symbol (a diamond by default) draws genes in the extended
    HLA/MHC region with that Plotly symbol instead of a circle so they stand
    out; set it to None to disable the marking entirely. It relies on
    extended_mhc_interval, which currently only supports genome build 19.
    """

    meta: Meta
    source: GeneManhattanSource
    sig_threshold: float | None = None
    point_size: int = 5
    colors: tuple[str, str] = ("#1f77b4", "#ff7f0e")
    sig_line_color: str = "red"
    title: str | None = None
    plotly_js_mode: bool | PlotlyWriteMode = "cdn"
    hla_marker_symbol: str | None = "diamond"

    @property
    def deps(self) -> list[Task]:
        return self.source.deps

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        data = self.source.load_df(fetch=fetch)
        max_p_value = self.source.max_p_value
        y_axis_start = (
            float(-np.log10(max_p_value)) if max_p_value is not None else None
        )
        hla_interval = (
            extended_mhc_interval(self.source.genome_build)
            if self.hla_marker_symbol is not None
            else None
        )
        fig = build_manhattan_plot(
            df=data.df,
            sig_threshold=self.sig_threshold,
            point_size=self.point_size,
            colors=self.colors,
            sig_line_color=self.sig_line_color,
            title=self.title,
            genome_build=self.source.genome_build,
            num_genes_for_correction=data.num_genes_for_correction,
            y_axis_start=y_axis_start,
            hla_interval=hla_interval,
            hla_marker_symbol=self.hla_marker_symbol,
        )
        out_path = scratch_dir / "gene_manhattan.html"
        fig.write_html(out_path, include_plotlyjs=self.plotly_js_mode)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        source: GeneManhattanSource,
        sig_threshold: float | None = None,
        title: str | None = None,
    ) -> "GeneManhattanPlotTask":
        meta = GWASPlotFileMeta(
            trait=source.trait,
            project=source.project,
            extension=".html",
            id=AssetId(asset_id),
        )
        return cls(
            meta=meta,
            source=source,
            sig_threshold=sig_threshold,
            title=title,
        )
