"""
Annotation of a table of genes, identified by Ensembl gene id, from static reference tables (no
network calls). Every join is on an Ensembl gene id or a UniProt accession, never on a gene
symbol: symbols go stale and collide.

- HGNC genes (joined on the Ensembl gene id): the gene's current approved symbol and its name,
  for coding and non-coding genes alike (HGNC_ENSEMBL_GENES).
- HGNC UniProt-to-gene table (joined on the Ensembl gene id): the UniProt accessions of the
  gene's protein products (HGNC_UNIPROT_GENES).
- UniProt (joined on the accession): protein names, function and subcellular location, from a
  bulk export of reviewed human UniProt entries (UNIPROT_LOOKUP), with field labels such as
  "FUNCTION:" and evidence code blocks dropped.

The output is the input table, row for row, with the columns of GENE_INFO_COLS appended. A gene
with several distinct values of a column has them joined by GENE_INFO_SEPARATOR, in order of
UniProt accession. A gene missing from HGNC, or without a protein
product, keeps its row with null columns; genes missing from HGNC are logged as a warning.
"""

from pathlib import Path
from typing import Self

import narwhals
import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import (
    CSVOutFormat,
    OutFormat,
    get_extension_and_read_spec_from_format,
    write_df_according_to_format,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.uniprot_comment_text import (
    null_if_blank,
    uniprot_comment_text,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gene_info_constants import (
    GENE_INFO_COLS,
    GENE_INFO_NAME_COL,
    GENE_INFO_PROTEIN_NAMES_COL,
    GENE_INFO_SEPARATOR,
    GENE_INFO_SUBCELLULAR_LOCALISATION_COL,
    GENE_INFO_SYMBOL_COL,
    GENE_INFO_UNIPROT_DESCRIPTION_COL,
    GENE_INFO_UNIPROT_ID_COL,
)
from mecfs_bio.constants.hgnc_constants import (
    HGNC_ENSEMBL_ID_COL,
    HGNC_NAME_COL,
    HGNC_SYMBOL_COL,
    HGNC_UNIPROT_ID_COL,
)
from mecfs_bio.constants.uniprot_constants import (
    UNIPROT_DAT_FUNCTION_COL,
    UNIPROT_DAT_ID_COL,
    UNIPROT_DAT_PROTEIN_NAMES_COL,
    UNIPROT_DAT_SUBCELLULAR_LOC_COL,
)

logger = structlog.getLogger()

_GENE_KEY_COL = "_ensembl_gene_id"
_PROTEIN_COLS = (
    GENE_INFO_UNIPROT_ID_COL,
    GENE_INFO_PROTEIN_NAMES_COL,
    GENE_INFO_UNIPROT_DESCRIPTION_COL,
    GENE_INFO_SUBCELLULAR_LOCALISATION_COL,
)


@frozen(slots=True)
class GeneInfoReferenceTasks:
    """
    The reference tables AnnotateGenesFromReferenceTask reads.

    hgnc_genes_task: HGNC genes by Ensembl gene id (HGNC_ENSEMBL_GENES).
    uniprot_genes_task: the normalized HGNC UniProt-to-gene table (HGNC_UNIPROT_GENES).
    uniprot_task: the bulk UniProt export (UNIPROT_LOOKUP).
    """

    hgnc_genes_task: Task
    uniprot_genes_task: Task
    uniprot_task: Task


def _read_table(fetch: Fetch, task: Task) -> pl.DataFrame:
    return (
        scan_dataframe_asset(fetch(task.asset_id), meta=task.meta).collect().to_polars()
    )


def _joined(column: str, alias: str) -> pl.Expr:
    """The distinct non-null values of column within a group, joined, or null when none."""
    return null_if_blank(
        pl.col(column)
        .drop_nulls()
        .unique(maintain_order=True)
        .str.join(GENE_INFO_SEPARATOR)
    ).alias(alias)


def _genes_by_ensembl_id(
    hgnc_genes: pl.DataFrame, ensembl_ids: list[str]
) -> pl.DataFrame:
    """Symbol and name of each gene, one row per Ensembl gene id."""
    found = (
        hgnc_genes.filter(pl.col(HGNC_ENSEMBL_ID_COL).is_in(ensembl_ids))
        .sort(HGNC_ENSEMBL_ID_COL, HGNC_SYMBOL_COL)
        .group_by(HGNC_ENSEMBL_ID_COL, maintain_order=True)
        .agg(
            _joined(HGNC_SYMBOL_COL, GENE_INFO_SYMBOL_COL),
            _joined(HGNC_NAME_COL, GENE_INFO_NAME_COL),
        )
        .rename({HGNC_ENSEMBL_ID_COL: _GENE_KEY_COL})
    )
    missing = sorted(set(ensembl_ids) - set(found[_GENE_KEY_COL].to_list()))
    if missing:
        logger.warning(f"{len(missing)} Ensembl gene ids are not in HGNC: {missing}")
    return found


def _cleaned_uniprot(uniprot: pl.DataFrame) -> pl.DataFrame:
    return uniprot.select(
        pl.col(UNIPROT_DAT_ID_COL).cast(pl.String).alias(GENE_INFO_UNIPROT_ID_COL),
        null_if_blank(pl.col(UNIPROT_DAT_PROTEIN_NAMES_COL).cast(pl.String)).alias(
            GENE_INFO_PROTEIN_NAMES_COL
        ),
        uniprot_comment_text(UNIPROT_DAT_FUNCTION_COL, "FUNCTION").alias(
            GENE_INFO_UNIPROT_DESCRIPTION_COL
        ),
        uniprot_comment_text(
            UNIPROT_DAT_SUBCELLULAR_LOC_COL, "SUBCELLULAR LOCATION"
        ).alias(GENE_INFO_SUBCELLULAR_LOCALISATION_COL),
    )


def _proteins_by_ensembl_id(
    uniprot_genes: pl.DataFrame, uniprot: pl.DataFrame, ensembl_ids: list[str]
) -> pl.DataFrame:
    """UniProt accessions and annotations of each gene's products, one row per Ensembl gene id."""
    return (
        uniprot_genes.filter(pl.col(HGNC_ENSEMBL_ID_COL).is_in(ensembl_ids))
        .select(
            pl.col(HGNC_ENSEMBL_ID_COL).alias(_GENE_KEY_COL),
            pl.col(HGNC_UNIPROT_ID_COL).alias(GENE_INFO_UNIPROT_ID_COL),
        )
        .unique()
        .join(_cleaned_uniprot(uniprot), on=GENE_INFO_UNIPROT_ID_COL, how="left")
        .sort(_GENE_KEY_COL, GENE_INFO_UNIPROT_ID_COL)
        .group_by(_GENE_KEY_COL, maintain_order=True)
        .agg(*[_joined(col, col) for col in _PROTEIN_COLS])
    )


def annotate_genes(
    table: pl.DataFrame,
    ensembl_id_col: str,
    hgnc_genes: pl.DataFrame,
    uniprot_genes: pl.DataFrame,
    uniprot: pl.DataFrame,
) -> pl.DataFrame:
    """The table, row for row, with the columns of GENE_INFO_COLS appended."""
    clashing = set(GENE_INFO_COLS) & set(table.columns)
    assert not clashing, f"the table already has gene info columns: {clashing}"
    assert table.schema[ensembl_id_col] == pl.String, (
        f"{ensembl_id_col} must hold Ensembl gene id strings, "
        f"got {table.schema[ensembl_id_col]}"
    )
    ensembl_ids = table[ensembl_id_col].drop_nulls().unique().to_list()
    genes = _genes_by_ensembl_id(hgnc_genes, ensembl_ids)
    proteins = _proteins_by_ensembl_id(uniprot_genes, uniprot, ensembl_ids)
    return (
        table.join(
            genes.join(proteins, on=_GENE_KEY_COL, how="left"),
            left_on=ensembl_id_col,
            right_on=_GENE_KEY_COL,
            how="left",
            maintain_order="left",
        )
        .select(*table.columns, *GENE_INFO_COLS)
        .cast({col: pl.String for col in GENE_INFO_COLS})
    )


@frozen(slots=True)
class AnnotateGenesFromReferenceTask(Task):
    """
    source_df_task: a table with a column of Ensembl gene ids.
    ensembl_id_col: that column.
    reference: the HGNC and UniProt reference tables.
    post_pipe: applied to the annotated table before writing it.
    """

    meta: Meta
    source_df_task: Task
    ensembl_id_col: str
    reference: GeneInfoReferenceTasks
    out_format: OutFormat = CSVOutFormat(",")
    post_pipe: DataProcessingPipe = IdentityPipe()

    @property
    def deps(self) -> list[Task]:
        return [
            self.source_df_task,
            self.reference.hgnc_genes_task,
            self.reference.uniprot_genes_task,
            self.reference.uniprot_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        table = annotate_genes(
            _read_table(fetch, self.source_df_task),
            ensembl_id_col=self.ensembl_id_col,
            hgnc_genes=_read_table(fetch, self.reference.hgnc_genes_task),
            uniprot_genes=_read_table(fetch, self.reference.uniprot_genes_task),
            uniprot=_read_table(fetch, self.reference.uniprot_task),
        )
        # No extension: the build system renames the file according to the task's meta.
        out_path = scratch_dir / str(self.meta.asset_id)
        write_df_according_to_format(
            df=self.post_pipe.process(narwhals.from_native(table).lazy()),
            out_path=out_path,
            out_format=self.out_format,
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        source_df_task: Task,
        ensembl_id_col: str,
        reference: GeneInfoReferenceTasks,
        post_pipe: DataProcessingPipe = IdentityPipe(),
        out_format: OutFormat = CSVOutFormat(","),
    ) -> Self:
        source_meta = source_df_task.meta
        extension, read_spec = get_extension_and_read_spec_from_format(out_format)
        meta: Meta
        if isinstance(source_meta, ResultTableMeta):
            meta = ResultTableMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                extension=extension,
                read_spec=read_spec,
            )
        elif isinstance(source_meta, ReferenceFileMeta):
            meta = ReferenceFileMeta(
                id=AssetId(asset_id),
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                extension=extension,
                read_spec=read_spec,
                sub_folder=source_meta.sub_folder,
            )
        else:
            raise ValueError(
                f"unsupported source meta type: {type(source_meta).__name__}"
            )
        return cls(
            meta=meta,
            source_df_task=source_df_task,
            ensembl_id_col=ensembl_id_col,
            reference=reference,
            post_pipe=post_pipe,
            out_format=out_format,
        )
