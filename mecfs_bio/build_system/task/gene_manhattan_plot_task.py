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
ID, and ``-log10(p)``.
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

    @abstractmethod
    def load_df(self, fetch: Fetch) -> pd.DataFrame:
        """Materialize a pandas DataFrame with columns ``chrom``, ``pos``, ``ensembl_id``, ``gene_name``, ``p_value``."""
        pass


@frozen
class MagmaGeneSource(GeneManhattanSource):
    """Load a Manhattan-plot table from a :class:`MagmaGeneAnalysisTask`.

    Chromosomal positions come from the MAGMA output itself. Human-readable
    gene names are joined in from ``gene_thesaurus_task`` by Ensembl ID. When
    a gene is missing from the thesaurus, the Ensembl ID is used as the
    display name.
    """

    magma_task: Task
    gene_thesaurus_task: Task

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

    def load_df(self, fetch: Fetch) -> pd.DataFrame:
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
    """

    table_task: Task
    gene_locations_task: Task
    gene_col: str
    p_col: str
    gene_id_kind: GeneIdKind = "ensembl_id"

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

    def load_df(self, fetch: Fetch) -> pd.DataFrame:
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
) -> go.Figure:
    """Construct a Plotly figure containing a gene-level Manhattan plot.

    Genes with non-positive or null p-values are dropped (``-log10`` is
    undefined). If ``sig_threshold`` is ``None``, a Bonferroni-corrected
    threshold ``0.05 / N_genes`` is used and a dashed horizontal line is drawn
    at the corresponding ``-log10(p)``.
    """
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
        sig_threshold = 0.05 / len(df)
    sig_y = float(-np.log10(sig_threshold))

    fig = go.Figure()
    for idx, chrom in enumerate(chroms):
        chrom_df = df[df[_CHROM] == chrom]
        color = colors[idx % 2]
        fig.add_trace(
            go.Scattergl(
                x=chrom_df["_x"],
                y=chrom_df["_mlog10p"],
                mode="markers",
                marker=dict(size=point_size, color=color),
                name=f"chr{chrom}",
                customdata=np.stack(
                    [
                        chrom_df[_GENE_NAME].astype(str).to_numpy(),
                        chrom_df[_ENSEMBL_ID].astype(str).to_numpy(),
                    ],
                    axis=-1,
                ),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "Ensembl: %{customdata[1]}<br>"
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
        yaxis=dict(
            title="-log<sub>10</sub>(p)",
            zeroline=False,
        ),
        plot_bgcolor="white",
        hovermode="closest",
    )
    return fig


@frozen
class GeneManhattanPlotTask(Task):
    """Create an interactive HTML gene-level Manhattan plot.

    Backed by Plotly's WebGL renderer (``Scattergl``) so that hover stays
    responsive at gene-scale point counts (~20k-30k).
    """

    meta: Meta
    source: GeneManhattanSource
    sig_threshold: float | None = None
    point_size: int = 5
    colors: tuple[str, str] = ("#1f77b4", "#ff7f0e")
    sig_line_color: str = "red"
    title: str | None = None
    plotly_js_mode: bool | PlotlyWriteMode = "cdn"

    @property
    def deps(self) -> list[Task]:
        return self.source.deps

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        df = self.source.load_df(fetch=fetch)
        fig = build_manhattan_plot(
            df=df,
            sig_threshold=self.sig_threshold,
            point_size=self.point_size,
            colors=self.colors,
            sig_line_color=self.sig_line_color,
            title=self.title,
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
