"""
Database annotations of the colocalizing proteins at the top loci of a coloc display table.

The loci are ranked by their best PP.H4 (the posterior probability of a shared causal variant)
and at most max_loci are kept. Each display row at those loci is split into one row per gene
symbol: a protein complex such as IL12A_IL12B is measured as one protein but has one symbol per
gene. Each symbol is mapped to an Ensembl gene id (by default via the Ensembl REST symbol
lookup), and the ids are annotated once, in a single batch, by gget.info (UniProt function, NCBI
summary, subcellular localisation).

A symbol that does not map keeps its row with null annotation columns, and a warning is logged.
An empty display table gives an empty table with the full schema.
"""

from pathlib import Path, PurePath
from typing import Protocol

import numpy as np
import pandas as pd
import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import write_parquet_table
from mecfs_bio.build_system.task.fetch_gget_info_task import fetch_gget_info
from mecfs_bio.build_system.task.protein_coloc.coloc_display_table_task import (
    ranked_display_loci,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gene_thesaurus_constants import (
    THESAURUS_CHROM_COL,
    THESAURUS_ENSEMBL_ID_COL,
    THESAURUS_GENE_NAME_COL,
    THESAURUS_PRIMARY_CHROMOSOMES,
)
from mecfs_bio.constants.protein_coloc_constants import (
    DISPLAY_GENE_COL,
    DISPLAY_LOCUS_COL,
    DISPLAY_PP_H4_COL,
    DISPLAY_PROTEIN_COL,
)

logger = structlog.getLogger()

GENE_SYMBOL_COL = "gene_symbol"
ENSEMBL_ID_COL = "ensembl_id"
UNIPROT_ID_COL = "uniprot_id"
PRIMARY_GENE_NAME_COL = "primary_gene_name"
NCBI_DESCRIPTION_COL = "ncbi_description"
UNIPROT_DESCRIPTION_COL = "uniprot_description"
SUBCELLULAR_LOCALISATION_COL = "subcellular_localisation"
PROTEIN_NAMES_COL = "protein_names"
ENSEMBL_DESCRIPTION_COL = "ensembl_description"
_COMPLEX_SEPARATOR = "_"
_LIST_JOIN = "; "

INFO_SCHEMA = pl.Schema(
    {
        ENSEMBL_ID_COL: pl.String,
        UNIPROT_ID_COL: pl.String,
        PRIMARY_GENE_NAME_COL: pl.String,
        PROTEIN_NAMES_COL: pl.String,
        ENSEMBL_DESCRIPTION_COL: pl.String,
        NCBI_DESCRIPTION_COL: pl.String,
        UNIPROT_DESCRIPTION_COL: pl.String,
        SUBCELLULAR_LOCALISATION_COL: pl.String,
    }
)

_DISPLAY_KEPT_SCHEMA = pl.Schema(
    {
        DISPLAY_GENE_COL: pl.String,
        DISPLAY_PROTEIN_COL: pl.String,
        DISPLAY_LOCUS_COL: pl.String,
        DISPLAY_PP_H4_COL: pl.Float64,
        GENE_SYMBOL_COL: pl.String,
    }
)


class InfoLookup(Protocol):
    def __call__(self, ensembl_ids: list[str]) -> pl.DataFrame:
        """One row per Ensembl id found, in INFO_SCHEMA."""
        ...


def _as_text(value: object) -> str | None:
    """A gget field as text: lists joined by "; ", and missing or blank values as None."""
    if isinstance(value, (list, tuple, np.ndarray)):
        parts = [str(v) for v in value if v is not None and not pd.isna(v)]
        text = _LIST_JOIN.join(parts)
    elif value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    else:
        text = str(value)
    return text if text.strip() else None


class GgetFetch(Protocol):
    def __call__(self, ensembl_ids: list[str]) -> pd.DataFrame:
        """gget.info output: one row per Ensembl id found, indexed by id."""
        ...


@frozen(slots=True)
class GgetInfoLookup:
    """gget.info (through fetch), with list-valued fields joined by "; " and blanks as None.

    gget's ncbi_description is often missing upstream, so protein_names and ensembl_description
    are kept as well.
    """

    fetch: GgetFetch = fetch_gget_info

    def __call__(self, ensembl_ids: list[str]) -> pl.DataFrame:
        frame = self.fetch(ensembl_ids)
        return pl.DataFrame(
            {
                ENSEMBL_ID_COL: [str(i) for i in frame.index],
                **{
                    name: [
                        _as_text(v)
                        for v in (
                            frame[name]
                            if name in frame.columns
                            else [None] * len(frame)
                        )
                    ]
                    for name in INFO_SCHEMA.names()
                    if name != ENSEMBL_ID_COL
                },
            },
            schema=INFO_SCHEMA,
        )


def _top_loci_by_symbol(display: pl.DataFrame, max_loci: int) -> pl.DataFrame:
    loci = ranked_display_loci(display, max_loci)
    return (
        display.filter(pl.col(DISPLAY_LOCUS_COL).is_in(loci))
        .with_columns(
            pl.col(DISPLAY_LOCUS_COL)
            .replace_strict(loci, list(range(len(loci))), return_dtype=pl.Int64)
            .alias("_locus_rank")
        )
        .sort(["_locus_rank", DISPLAY_PP_H4_COL], descending=[False, True])
        .with_columns(
            pl.col(DISPLAY_GENE_COL)
            .str.split(_COMPLEX_SEPARATOR)
            .alias(GENE_SYMBOL_COL)
        )
        .explode(GENE_SYMBOL_COL)
        .select(_DISPLAY_KEPT_SCHEMA.names())
        .cast(_DISPLAY_KEPT_SCHEMA)
    )


def _primary_assembly_ids(thesaurus: pl.DataFrame) -> pl.DataFrame:
    """Distinct (gene symbol, Ensembl id) pairs of the genes on primary-assembly chromosomes."""
    return (
        thesaurus.filter(
            pl.col(THESAURUS_CHROM_COL)
            .cast(pl.String)
            .is_in(sorted(THESAURUS_PRIMARY_CHROMOSOMES))
        )
        .select(
            pl.col(THESAURUS_GENE_NAME_COL).cast(pl.String),
            pl.col(THESAURUS_ENSEMBL_ID_COL).cast(pl.String),
        )
        .drop_nulls()
        .unique()
    )


def _ensembl_ids(
    symbols: list[str], primary_ids: pl.DataFrame
) -> dict[str, str | None]:
    wanted = list(dict.fromkeys(symbols))
    matches = primary_ids.filter(pl.col(THESAURUS_GENE_NAME_COL).is_in(wanted))
    ambiguous = (
        matches.group_by(THESAURUS_GENE_NAME_COL)
        .len()
        .filter(pl.col("len") > 1)[THESAURUS_GENE_NAME_COL]
        .to_list()
    )
    assert not ambiguous, (
        f"gene symbols with more than one primary-assembly Ensembl id: {ambiguous}"
    )
    found = dict(matches.iter_rows())
    mapping = {symbol: found.get(symbol) for symbol in wanted}
    for symbol, ensembl_id in mapping.items():
        if ensembl_id is None:
            logger.warning(f"gene symbol {symbol} is not in the gene thesaurus")
    return mapping


def _annotate(
    rows: pl.DataFrame, primary_ids: pl.DataFrame, info_lookup: InfoLookup
) -> pl.DataFrame:
    mapping = _ensembl_ids(rows[GENE_SYMBOL_COL].to_list(), primary_ids)
    with_ids = rows.with_columns(
        pl.col(GENE_SYMBOL_COL)
        .replace_strict(mapping, default=None, return_dtype=pl.String)
        .alias(ENSEMBL_ID_COL)
    )
    ids = [i for i in dict.fromkeys(mapping.values()) if i is not None]
    info = info_lookup(ids) if ids else pl.DataFrame(schema=INFO_SCHEMA)
    return with_ids.join(
        info.cast(INFO_SCHEMA).unique(ENSEMBL_ID_COL, keep="first"),
        on=ENSEMBL_ID_COL,
        how="left",
        maintain_order="left",
    )


def _read_table(fetch: Fetch, task: Task) -> pl.DataFrame:
    return (
        scan_dataframe_asset(fetch(task.asset_id), meta=task.meta).collect().to_polars()
    )


@frozen(slots=True)
class ColocProteinInfoTask(Task):
    """
    display_task: a ColocDisplayTableTask (any Task producing its table).
    gene_thesaurus_task: the gene thesaurus, mapping gene symbols to Ensembl gene ids.
    max_loci: only the proteins at the top max_loci loci (best PP.H4 first) are annotated.
    info_lookup: annotates Ensembl gene ids.
    """

    meta: Meta
    display_task: Task
    gene_thesaurus_task: Task
    max_loci: int
    info_lookup: InfoLookup

    def __attrs_post_init__(self) -> None:
        assert self.max_loci >= 1, f"max_loci must be at least 1: {self.max_loci}"

    @property
    def deps(self) -> list[Task]:
        return [self.display_task, self.gene_thesaurus_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        display = _read_table(fetch, self.display_task)
        table = _annotate(
            _top_loci_by_symbol(display, self.max_loci),
            _primary_assembly_ids(_read_table(fetch, self.gene_thesaurus_task)),
            self.info_lookup,
        )
        out_path = scratch_dir / f"{self.meta.asset_id}.parquet"
        write_parquet_table(
            table=table.to_arrow(),
            out_path=out_path,
            compression="zstd",
            compression_level=None,
            byte_stream_split_columns=[DISPLAY_PP_H4_COL],
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        display_task: Task,
        gene_thesaurus_task: Task,
        max_loci: int = 10,
        info_lookup: InfoLookup = GgetInfoLookup(),
    ) -> "ColocProteinInfoTask":
        source_meta = display_task.meta
        assert isinstance(source_meta, ResultTableMeta), (
            f"display_task must carry ResultTableMeta, got {type(source_meta).__name__}"
        )
        return cls(
            meta=ResultTableMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                sub_dir=PurePath("analysis"),
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
            display_task=display_task,
            gene_thesaurus_task=gene_thesaurus_task,
            max_loci=max_loci,
            info_lookup=info_lookup,
        )
