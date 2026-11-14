"""
Annotations of the colocalizing proteins at the top loci of a coloc display table, from static
in-repo reference tables (no network calls). Every join is on a UniProt accession, never on a
gene symbol: symbols go stale (UKB-PPP lists WARS, now WARS1) and collide.

The loci are ranked by their best PP.H4 (the posterior probability of a shared causal variant)
and at most max_loci are kept. Each display row at those loci becomes one row per UniProt
accession of the protein, as the proteomics source's protein table lists them (via
protein_accessions_reader): a protein complex such as EBI3_IL27 has one accession per subunit.

- Genes (joined on the accession): a normalized UniProt-to-gene table from the HGNC complete set
  (HGNC_UNIPROT_GENES) gives the gene's current approved symbol, its name and its Ensembl gene
  id. An accession encoded by duplicated genes (for example P59665, by DEFA1 and DEFA1B) gets
  one row per gene.
- UniProt (joined on the accession): protein names, function and subcellular location, from a
  bulk export of reviewed human UniProt entries. Field labels such as "FUNCTION:" are dropped.
  Evidence code blocks such as {ECO:0000269|PubMed:123, ECO:0000250|UniProtKB:Q8R4Y0} are moved
  out of the text into one evidence column, grouped by readable label with their sources, e.g.
  "experimental (PubMed:123); by similarity (UniProtKB:Q8R4Y0)". An ECO (Evidence and Conclusion
  Ontology) code without a label in _ECO_LABELS is shown as the code itself.

The output has the display columns (protein, assay_id, locus, PP.H4), the ids (UniProt accession,
approved symbol, Ensembl gene id) and the descriptive columns (gene name, protein names, UniProt
function, subcellular location as "locale", and evidence).

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
    HGNC_NAME_COL,
    HGNC_SYMBOL_COL,
    HGNC_UNIPROT_ID_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    DISPLAY_ASSAY_ID_COL,
    DISPLAY_LOCUS_COL,
    DISPLAY_PP_H4_COL,
    DISPLAY_PROTEIN_COL,
    INFO_ENSEMBL_ID_COL,
    INFO_EVIDENCE_COL,
    INFO_GENE_NAME_COL,
    INFO_GENE_SYMBOL_COL,
    INFO_PROTEIN_NAMES_COL,
    INFO_SUBCELLULAR_LOCALISATION_COL,
    INFO_UNIPROT_DESCRIPTION_COL,
    INFO_UNIPROT_ID_COL,
)
from mecfs_bio.constants.uniprot_constants import (
    UNIPROT_DAT_FUNCTION_COL,
    UNIPROT_DAT_ID_COL,
    UNIPROT_DAT_PROTEIN_NAMES_COL,
    UNIPROT_DAT_SUBCELLULAR_LOC_COL,
)

logger = structlog.getLogger()

UNIPROT_ID_COL = INFO_UNIPROT_ID_COL
GENE_SYMBOL_COL = INFO_GENE_SYMBOL_COL
GENE_NAME_COL = INFO_GENE_NAME_COL
ENSEMBL_ID_COL = INFO_ENSEMBL_ID_COL
PROTEIN_NAMES_COL = INFO_PROTEIN_NAMES_COL
UNIPROT_DESCRIPTION_COL = INFO_UNIPROT_DESCRIPTION_COL
SUBCELLULAR_LOCALISATION_COL = INFO_SUBCELLULAR_LOCALISATION_COL
EVIDENCE_COL = INFO_EVIDENCE_COL
_ROW_INDEX_COL = "_row"

# Consecutive comments of one field are joined by "; " after a full stop; inside a comment, "; "
# separates items (for example "Cell membrane; Single-pass type I membrane protein").
_COMMENT_SEPARATOR_PATTERN = r"\.;\s+"

# One evidence block, e.g. {ECO:0000269|PubMed:5, ECO:0000250|UniProtKB:Q8R4Y0}, with the space
# before it. Inside a block, ", " separates items; each item is "ECO:<code>" optionally followed
# by "|<source>".
_EVIDENCE_BLOCK_PATTERN = r"\s*\{ECO:[^}]*\}"
# A block between two full stops ("Binds X. {ECO:...}."), dropped with one of them.
_EVIDENCE_SENTENCE_PATTERN = r"\.\s*\{ECO:[^}]*\}\."
_EVIDENCE_ITEM_SEPARATOR = ", "
_EVIDENCE_SOURCE_SEPARATOR = "|"

# Readable labels of the ECO codes UniProt uses in its comment fields.
_ECO_LABELS: dict[str, str] = {
    "ECO:0000269": "experimental",
    "ECO:0000250": "by similarity",
    "ECO:0000255": "sequence model",
    "ECO:0000256": "sequence model, automatic",
    "ECO:0000303": "non-traceable author statement",
    "ECO:0000305": "curator inference",
    "ECO:0000312": "imported",
    "ECO:0000313": "imported, automatic",
    "ECO:0007744": "combinatorial",
    "ECO:0007829": "combinatorial, automatic",
}

_ROWS_SCHEMA = pl.Schema(
    {
        DISPLAY_PROTEIN_COL: pl.String,
        DISPLAY_ASSAY_ID_COL: pl.String,
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
        EVIDENCE_COL: pl.String,
    }
)

OUTPUT_SCHEMA = pl.Schema(
    {
        **_ROWS_SCHEMA,
        GENE_SYMBOL_COL: pl.String,
        ENSEMBL_ID_COL: pl.String,
        GENE_NAME_COL: pl.String,
        PROTEIN_NAMES_COL: pl.String,
        UNIPROT_DESCRIPTION_COL: pl.String,
        SUBCELLULAR_LOCALISATION_COL: pl.String,
        EVIDENCE_COL: pl.String,
    }
)


def _read_table(fetch: Fetch, task: Task) -> pl.DataFrame:
    return (
        scan_dataframe_asset(fetch(task.asset_id), meta=task.meta).collect().to_polars()
    )


def _null_if_blank(expr: pl.Expr) -> pl.Expr:
    stripped = expr.str.strip_chars()
    return pl.when(stripped == "").then(None).otherwise(stripped)


def _without_label(column: str, label: str) -> pl.Expr:
    """A UniProt comment field with its label (e.g. "FUNCTION:") and evidence blocks dropped."""
    return _null_if_blank(
        pl.col(column)
        .cast(pl.String)
        .str.replace_all(f"{label}:", "", literal=True)
        .str.replace_all(_EVIDENCE_SENTENCE_PATTERN, ".")
        .str.replace_all(_EVIDENCE_BLOCK_PATTERN, "")
        .str.replace_all(_COMMENT_SEPARATOR_PATTERN, ". ")
        .str.replace_all(r"\s+", " ")
    )


def _evidence_blocks(columns: list[str]) -> pl.Expr:
    """All evidence blocks of the given UniProt comment fields, in order, as a list."""
    return pl.concat_list(
        [
            pl.col(column)
            .cast(pl.String)
            .str.extract_all(_EVIDENCE_BLOCK_PATTERN)
            .fill_null(pl.lit([], dtype=pl.List(pl.String)))
            for column in columns
        ]
    )


def _evidence_label(blocks: list[str]) -> str | None:
    """Readable evidence for one entry, e.g. "experimental (PubMed:5, PubMed:6); curator inference".

    Items are grouped by the label of their ECO code, in order of first appearance; each label
    lists its distinct sources. None when there are no evidence blocks.
    """
    sources_by_label: dict[str, list[str]] = {}
    for block in blocks:
        for item in block.strip().strip("{}").split(_EVIDENCE_ITEM_SEPARATOR):
            code, _, source = item.strip().partition(_EVIDENCE_SOURCE_SEPARATOR)
            sources = sources_by_label.setdefault(_ECO_LABELS.get(code, code), [])
            if source and source not in sources:
                sources.append(source)
    if not sources_by_label:
        return None
    return "; ".join(
        f"{label} ({', '.join(sources)})" if sources else label
        for label, sources in sources_by_label.items()
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
    missing = set(top[DISPLAY_ASSAY_ID_COL].to_list()) - set(
        protein_accessions[PROTEIN_ACCESSIONS_PROTEIN_ID_COL].to_list()
    )
    assert not missing, f"display proteins missing from the protein table: {missing}"
    return (
        top.join(
            protein_accessions,
            left_on=DISPLAY_ASSAY_ID_COL,
            right_on=PROTEIN_ACCESSIONS_PROTEIN_ID_COL,
            how="left",
            maintain_order="left",
        )
        .rename({PROTEIN_ACCESSIONS_UNIPROT_COL: UNIPROT_ID_COL})
        .select(_ROWS_SCHEMA.names())
        .cast(_ROWS_SCHEMA)
    )


def _genes_by_accession(
    uniprot_genes: pl.DataFrame, accessions: list[str]
) -> pl.DataFrame:
    """The genes encoding each accession, one row per (accession, gene)."""
    found = uniprot_genes.filter(pl.col(HGNC_UNIPROT_ID_COL).is_in(accessions)).select(
        pl.col(HGNC_UNIPROT_ID_COL).alias(UNIPROT_ID_COL),
        pl.col(HGNC_SYMBOL_COL).alias(GENE_SYMBOL_COL),
        pl.col(HGNC_ENSEMBL_ID_COL).alias(ENSEMBL_ID_COL),
        pl.col(HGNC_NAME_COL).alias(GENE_NAME_COL),
    )
    for accession in sorted(set(accessions) - set(found[UNIPROT_ID_COL].to_list())):
        logger.warning(f"UniProt accession {accession} is not in the HGNC table")
    return found.cast(_GENES_SCHEMA)


def _uniprot_by_accession(uniprot: pl.DataFrame, accessions: list[str]) -> pl.DataFrame:
    found = uniprot.filter(pl.col(UNIPROT_DAT_ID_COL).cast(pl.String).is_in(accessions))
    evidence = [
        _evidence_label(blocks)
        for blocks in found.select(
            _evidence_blocks(
                [UNIPROT_DAT_FUNCTION_COL, UNIPROT_DAT_SUBCELLULAR_LOC_COL]
            )
        )
        .to_series()
        .to_list()
    ]
    found = found.select(
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
    ).with_columns(pl.Series(EVIDENCE_COL, evidence, dtype=pl.String))
    for accession in sorted(set(accessions) - set(found[UNIPROT_ID_COL].to_list())):
        logger.warning(f"UniProt accession {accession} is not in the UniProt table")
    return found.cast(_UNIPROT_SCHEMA)


def _annotate(
    rows: pl.DataFrame, uniprot_genes: pl.DataFrame, uniprot: pl.DataFrame
) -> pl.DataFrame:
    accessions = rows[UNIPROT_ID_COL].drop_nulls().unique().to_list()
    return (
        rows.with_row_index(_ROW_INDEX_COL)
        .join(
            _genes_by_accession(uniprot_genes, accessions),
            on=UNIPROT_ID_COL,
            how="left",
        )
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
    uniprot_genes_task: the normalized HGNC UniProt-to-gene table (HGNC_UNIPROT_GENES).
    max_loci: only the proteins at the top max_loci loci (best PP.H4 first) are annotated.
    """

    meta: Meta
    display_task: Task
    protein_table_task: Task
    protein_accessions_reader: ProteinAccessionsReader
    uniprot_task: Task
    uniprot_genes_task: Task
    max_loci: int

    def __attrs_post_init__(self) -> None:
        assert self.max_loci >= 1, f"max_loci must be at least 1: {self.max_loci}"

    @property
    def deps(self) -> list[Task]:
        return [
            self.display_task,
            self.protein_table_task,
            self.uniprot_task,
            self.uniprot_genes_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        rows = _top_loci_by_accession(
            _read_table(fetch, self.display_task),
            self.protein_accessions_reader(_read_table(fetch, self.protein_table_task)),
            self.max_loci,
        )
        table = _annotate(
            rows,
            uniprot_genes=_read_table(fetch, self.uniprot_genes_task),
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
        uniprot_genes_task: Task,
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
            uniprot_genes_task=uniprot_genes_task,
            max_loci=max_loci,
        )
