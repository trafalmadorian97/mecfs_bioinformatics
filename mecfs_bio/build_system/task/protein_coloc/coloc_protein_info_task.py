"""
Annotations of the colocalizing proteins at the top loci of a coloc display table, from static
in-repo reference tables (no network calls). Every join is on a UniProt accession, never on a
gene symbol: symbols go stale (UKB-PPP lists WARS, now WARS1) and collide.

The loci are ranked by their best PP.H4 (the posterior probability of a shared causal variant)
and at most max_loci are kept. Each display row at those loci becomes one row per UniProt
accession of the protein, as the proteomics source's protein table lists them (via
protein_accessions_reader): a protein complex such as EBI3_IL27 has one accession per subunit.

- Genes (joined on the accession): the HGNC complete set gives the gene's current approved
  symbol, its name and its Ensembl gene id. An accession encoded by duplicated genes (for
  example P59665, by DEFA1 and DEFA1B) gets one row per gene. Readthrough genes (a transcript
  spanning two adjacent genes, such as PTGES3L-AARSD1) are skipped.
- UniProt (joined on the accession): protein names, function and subcellular location, from a
  bulk export of reviewed human UniProt entries. Field labels such as "FUNCTION:" and evidence
  codes such as {ECO:0000269|PubMed:123} are dropped.

An accession missing from HGNC or from UniProt keeps its row with null columns, and a warning is
logged. A display protein missing from the protein table fails the task. An empty display table
gives an empty table with the full schema.
"""

from pathlib import Path, PurePath

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
from mecfs_bio.build_system.task.protein_coloc.coloc_display_table_task import (
    ranked_display_loci,
)
from mecfs_bio.build_system.task.protein_coloc.protein_accessions_reader import (
    PROTEIN_ACCESSIONS_PROTEIN_ID_COL,
    PROTEIN_ACCESSIONS_UNIPROT_COL,
    ProteinAccessionsReader,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.hgnc_constants import (
    HGNC_ENSEMBL_ID_COL,
    HGNC_LIST_SEPARATOR,
    HGNC_LOCUS_TYPE_COL,
    HGNC_NAME_COL,
    HGNC_READTHROUGH_LOCUS_TYPE,
    HGNC_SYMBOL_COL,
    HGNC_UNIPROT_IDS_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    DISPLAY_GENE_COL,
    DISPLAY_LOCUS_COL,
    DISPLAY_PP_H4_COL,
    DISPLAY_PROTEIN_COL,
)
from mecfs_bio.constants.uniprot_constants import (
    UNIPROT_DAT_FUNCTION_COL,
    UNIPROT_DAT_ID_COL,
    UNIPROT_DAT_PROTEIN_NAMES_COL,
    UNIPROT_DAT_SUBCELLULAR_LOC_COL,
)

logger = structlog.getLogger()

UNIPROT_ID_COL = PROTEIN_ACCESSIONS_UNIPROT_COL
GENE_SYMBOL_COL = "gene_symbol"
GENE_NAME_COL = "gene_name"
ENSEMBL_ID_COL = "ensembl_id"
PROTEIN_NAMES_COL = "protein_names"
UNIPROT_DESCRIPTION_COL = "uniprot_description"
SUBCELLULAR_LOCALISATION_COL = "subcellular_localisation"
_ROW_INDEX_COL = "_row"

# An evidence code, e.g. " {ECO:0000269|PubMed:12805376}". UniProt writes it between a
# sentence's full stop and the comment's closing full stop, so dropping it can leave "..".
_ECO_PATTERN = r"\s*\{ECO:[^}]*\}"
# Consecutive comments of one field are joined by "; " after a full stop; inside a comment, "; "
# separates items (for example "Cell membrane; Single-pass type I membrane protein").
_COMMENT_SEPARATOR_PATTERN = r"\.;\s+"

_ROWS_SCHEMA = pl.Schema(
    {
        DISPLAY_GENE_COL: pl.String,
        DISPLAY_PROTEIN_COL: pl.String,
        DISPLAY_LOCUS_COL: pl.String,
        DISPLAY_PP_H4_COL: pl.Float64,
        UNIPROT_ID_COL: pl.String,
    }
)

_GENES_SCHEMA = pl.Schema(
    {
        UNIPROT_ID_COL: pl.String,
        GENE_SYMBOL_COL: pl.String,
        ENSEMBL_ID_COL: pl.String,
        GENE_NAME_COL: pl.String,
    }
)

_UNIPROT_SCHEMA = pl.Schema(
    {
        UNIPROT_ID_COL: pl.String,
        PROTEIN_NAMES_COL: pl.String,
        UNIPROT_DESCRIPTION_COL: pl.String,
        SUBCELLULAR_LOCALISATION_COL: pl.String,
    }
)

# The annotation columns added to each display row.
INFO_SCHEMA = pl.Schema(
    {
        UNIPROT_ID_COL: pl.String,
        GENE_SYMBOL_COL: pl.String,
        ENSEMBL_ID_COL: pl.String,
        GENE_NAME_COL: pl.String,
        PROTEIN_NAMES_COL: pl.String,
        UNIPROT_DESCRIPTION_COL: pl.String,
        SUBCELLULAR_LOCALISATION_COL: pl.String,
    }
)

OUTPUT_SCHEMA = pl.Schema({**_ROWS_SCHEMA, **INFO_SCHEMA})


def _read_table(fetch: Fetch, task: Task) -> pl.DataFrame:
    return (
        scan_dataframe_asset(fetch(task.asset_id), meta=task.meta).collect().to_polars()
    )


def _null_if_blank(expr: pl.Expr) -> pl.Expr:
    stripped = expr.str.strip_chars()
    return pl.when(stripped == "").then(None).otherwise(stripped)


def _without_label(column: str, label: str) -> pl.Expr:
    """A UniProt comment field with its label (e.g. "FUNCTION:") and evidence codes dropped."""
    return _null_if_blank(
        pl.col(column)
        .cast(pl.String)
        .str.replace_all(f"{label}:", "", literal=True)
        .str.replace_all(_ECO_PATTERN, "")
        .str.replace_all(r"\.\.", ".")
        .str.replace_all(_COMMENT_SEPARATOR_PATTERN, ". ")
        .str.replace_all(r"\s+", " ")
    )


def _top_loci_by_accession(
    display: pl.DataFrame, protein_accessions: pl.DataFrame, max_loci: int
) -> pl.DataFrame:
    loci = ranked_display_loci(display, max_loci)
    top = (
        display.filter(pl.col(DISPLAY_LOCUS_COL).is_in(loci))
        .with_columns(
            pl.col(DISPLAY_LOCUS_COL)
            .replace_strict(loci, list(range(len(loci))), return_dtype=pl.Int64)
            .alias("_locus_rank")
        )
        .sort(["_locus_rank", DISPLAY_PP_H4_COL], descending=[False, True])
    )
    missing = set(top[DISPLAY_PROTEIN_COL].to_list()) - set(
        protein_accessions[PROTEIN_ACCESSIONS_PROTEIN_ID_COL].to_list()
    )
    assert not missing, f"display proteins missing from the protein table: {missing}"
    return (
        top.join(
            protein_accessions,
            left_on=DISPLAY_PROTEIN_COL,
            right_on=PROTEIN_ACCESSIONS_PROTEIN_ID_COL,
            how="left",
            maintain_order="left",
        )
        .select(_ROWS_SCHEMA.names())
        .cast(_ROWS_SCHEMA)
    )


def _genes_by_accession(hgnc: pl.DataFrame, accessions: list[str]) -> pl.DataFrame:
    """The non-readthrough HGNC genes encoding each accession, one row per (accession, gene)."""
    found = (
        hgnc.filter(pl.col(HGNC_LOCUS_TYPE_COL) != HGNC_READTHROUGH_LOCUS_TYPE)
        .select(
            pl.col(HGNC_UNIPROT_IDS_COL)
            .cast(pl.String)
            .str.split(HGNC_LIST_SEPARATOR)
            .alias(UNIPROT_ID_COL),
            pl.col(HGNC_SYMBOL_COL).cast(pl.String).alias(GENE_SYMBOL_COL),
            _null_if_blank(pl.col(HGNC_ENSEMBL_ID_COL).cast(pl.String)).alias(
                ENSEMBL_ID_COL
            ),
            _null_if_blank(pl.col(HGNC_NAME_COL).cast(pl.String)).alias(GENE_NAME_COL),
        )
        .explode(UNIPROT_ID_COL)
        .with_columns(pl.col(UNIPROT_ID_COL).str.strip_chars())
        .filter(pl.col(UNIPROT_ID_COL).is_in(accessions))
        .unique()
    )
    for accession in sorted(set(accessions) - set(found[UNIPROT_ID_COL].to_list())):
        logger.warning(f"UniProt accession {accession} is not in the HGNC table")
    return found.cast(_GENES_SCHEMA)


def _uniprot_by_accession(uniprot: pl.DataFrame, accessions: list[str]) -> pl.DataFrame:
    found = uniprot.filter(
        pl.col(UNIPROT_DAT_ID_COL).cast(pl.String).is_in(accessions)
    ).select(
        pl.col(UNIPROT_DAT_ID_COL).cast(pl.String).alias(UNIPROT_ID_COL),
        _null_if_blank(pl.col(UNIPROT_DAT_PROTEIN_NAMES_COL).cast(pl.String)).alias(
            PROTEIN_NAMES_COL
        ),
        _without_label(UNIPROT_DAT_FUNCTION_COL, "FUNCTION").alias(
            UNIPROT_DESCRIPTION_COL
        ),
        _without_label(UNIPROT_DAT_SUBCELLULAR_LOC_COL, "SUBCELLULAR LOCATION").alias(
            SUBCELLULAR_LOCALISATION_COL
        ),
    )
    for accession in sorted(set(accessions) - set(found[UNIPROT_ID_COL].to_list())):
        logger.warning(f"UniProt accession {accession} is not in the UniProt table")
    return found.cast(_UNIPROT_SCHEMA)


def _annotate(
    rows: pl.DataFrame, hgnc: pl.DataFrame, uniprot: pl.DataFrame
) -> pl.DataFrame:
    accessions = rows[UNIPROT_ID_COL].drop_nulls().unique().to_list()
    return (
        rows.with_row_index(_ROW_INDEX_COL)
        .join(_genes_by_accession(hgnc, accessions), on=UNIPROT_ID_COL, how="left")
        .join(_uniprot_by_accession(uniprot, accessions), on=UNIPROT_ID_COL, how="left")
        .sort([_ROW_INDEX_COL, GENE_SYMBOL_COL], nulls_last=True)
        .select(OUTPUT_SCHEMA.names())
        .cast(OUTPUT_SCHEMA)
    )


@frozen(slots=True)
class ColocProteinInfoTask(Task):
    """
    display_task: a ColocDisplayTableTask (any Task producing its table).
    protein_table_task: the proteomics source's protein table (for UKB-PPP, the Sun et al. 2023
        ST3 sheet), read by protein_accessions_reader.
    protein_accessions_reader: lists the UniProt accessions of each protein.
    uniprot_task: the bulk UniProt export (for example UNIPROT_LOOKUP).
    hgnc_task: the HGNC complete set, mapping UniProt accessions to genes.
    max_loci: only the proteins at the top max_loci loci (best PP.H4 first) are annotated.
    """

    meta: Meta
    display_task: Task
    protein_table_task: Task
    protein_accessions_reader: ProteinAccessionsReader
    uniprot_task: Task
    hgnc_task: Task
    max_loci: int

    def __attrs_post_init__(self) -> None:
        assert self.max_loci >= 1, f"max_loci must be at least 1: {self.max_loci}"

    @property
    def deps(self) -> list[Task]:
        return [
            self.display_task,
            self.protein_table_task,
            self.uniprot_task,
            self.hgnc_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        rows = _top_loci_by_accession(
            _read_table(fetch, self.display_task),
            self.protein_accessions_reader(_read_table(fetch, self.protein_table_task)),
            self.max_loci,
        )
        table = _annotate(
            rows,
            hgnc=_read_table(fetch, self.hgnc_task),
            uniprot=_read_table(fetch, self.uniprot_task),
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
        protein_table_task: Task,
        protein_accessions_reader: ProteinAccessionsReader,
        uniprot_task: Task,
        hgnc_task: Task,
        max_loci: int = 10,
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
            protein_table_task=protein_table_task,
            protein_accessions_reader=protein_accessions_reader,
            uniprot_task=uniprot_task,
            hgnc_task=hgnc_task,
            max_loci=max_loci,
        )
