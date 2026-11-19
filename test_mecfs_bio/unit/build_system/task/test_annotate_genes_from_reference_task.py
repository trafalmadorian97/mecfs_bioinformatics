from pathlib import Path, PurePath

import polars as pl

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.task.annotate_genes_from_reference_task import (
    AnnotateGenesFromReferenceTask,
    GeneInfoReferenceTasks,
)
from mecfs_bio.build_system.task.dataframe_output import ParquetOutFormat
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import make_wf
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

_ID_COL = "Ensembl ID"

# ENSG_A: one protein. ENSG_B: two proteins, one missing from UniProt. ENSG_LNC: non-coding,
# listed under two symbols. ENSG_GONE: not in HGNC.
_SOURCE_IDS = ["ENSG_B", "ENSG_A", None, "ENSG_LNC", "ENSG_GONE"]

# (Ensembl id, symbol, name)
_HGNC_GENES_ROWS = [
    ("ENSG_A", "AAA", "gene A"),
    ("ENSG_B", "BBB", "gene B"),
    ("ENSG_LNC", "LINC2", "lncRNA 2"),
    ("ENSG_LNC", "LINC1", "lncRNA 1"),
    ("ENSG_OTHER", "OTHER", "unrelated"),
]
# (accession, symbol, name, Ensembl id)
_UNIPROT_GENES_ROWS = [
    ("P00001", "AAA", "gene A", "ENSG_A"),
    ("P00003", "BBB", "gene B", "ENSG_B"),
    ("P00002", "BBB", "gene B", "ENSG_B"),
]
# (accession, protein names, function, subcellular location)
_UNIPROT_ROWS = [
    (
        "P00001",
        "Protein A",
        "FUNCTION: Binds things. {ECO:0000269|PubMed:1}.",
        "SUBCELLULAR LOCATION: Secreted {ECO:0000305}.",
    ),
    ("P00002", "Protein B", "FUNCTION: Does B.", None),
]


def _table_task(asset_id: str) -> FakeTask:
    return FakeTask(
        ResultTableMeta(
            id=AssetId(asset_id),
            trait="synthetic",
            project="test",
            sub_dir=PurePath("analysis"),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
    )


def _write(tmp_path: Path, name: str, df: pl.DataFrame) -> tuple[AssetId, FileAsset]:
    path = tmp_path / f"{name}.parquet"
    df.write_parquet(path)
    return AssetId(name), FileAsset(path)


def _annotated(tmp_path: Path) -> pl.DataFrame:
    assets = dict(
        [
            _write(
                tmp_path,
                "source",
                pl.DataFrame(
                    {_ID_COL: _SOURCE_IDS, "rank": list(range(len(_SOURCE_IDS)))},
                    schema={_ID_COL: pl.String, "rank": pl.Int64},
                ),
            ),
            _write(
                tmp_path,
                "hgnc_genes",
                pl.DataFrame(
                    _HGNC_GENES_ROWS,
                    schema=[HGNC_ENSEMBL_ID_COL, HGNC_SYMBOL_COL, HGNC_NAME_COL],
                    orient="row",
                ),
            ),
            _write(
                tmp_path,
                "uniprot_genes",
                pl.DataFrame(
                    _UNIPROT_GENES_ROWS,
                    schema=[
                        HGNC_UNIPROT_ID_COL,
                        HGNC_SYMBOL_COL,
                        HGNC_NAME_COL,
                        HGNC_ENSEMBL_ID_COL,
                    ],
                    orient="row",
                ),
            ),
            _write(
                tmp_path,
                "uniprot",
                pl.DataFrame(
                    _UNIPROT_ROWS,
                    schema=[
                        UNIPROT_DAT_ID_COL,
                        UNIPROT_DAT_PROTEIN_NAMES_COL,
                        UNIPROT_DAT_FUNCTION_COL,
                        UNIPROT_DAT_SUBCELLULAR_LOC_COL,
                    ],
                    orient="row",
                ),
            ),
        ]
    )
    task = AnnotateGenesFromReferenceTask.create(
        asset_id="annotated",
        source_df_task=_table_task("source"),
        ensembl_id_col=_ID_COL,
        reference=GeneInfoReferenceTasks(
            hgnc_genes_task=_table_task("hgnc_genes"),
            uniprot_genes_task=_table_task("uniprot_genes"),
            uniprot_task=_table_task("uniprot"),
        ),
        out_format=ParquetOutFormat(),
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(
        scratch, fetch=lambda asset_id: assets[asset_id], wf=make_wf()
    )
    assert isinstance(result, FileAsset)
    return pl.read_parquet(result.path)


def test_keeps_source_rows_in_order_and_appends_gene_info(tmp_path: Path) -> None:
    annotated = _annotated(tmp_path)
    assert annotated[_ID_COL].to_list() == _SOURCE_IDS
    assert annotated.columns == [_ID_COL, "rank", *GENE_INFO_COLS]


def test_single_protein_gene_has_cleaned_uniprot_text(tmp_path: Path) -> None:
    row = _annotated(tmp_path).filter(pl.col(_ID_COL) == "ENSG_A").row(0, named=True)
    assert row[GENE_INFO_SYMBOL_COL] == "AAA"
    assert row[GENE_INFO_NAME_COL] == "gene A"
    assert row[GENE_INFO_UNIPROT_ID_COL] == "P00001"
    assert row[GENE_INFO_PROTEIN_NAMES_COL] == "Protein A"
    assert row[GENE_INFO_UNIPROT_DESCRIPTION_COL] == "Binds things."
    assert row[GENE_INFO_SUBCELLULAR_LOCALISATION_COL] == "Secreted."


def test_multi_protein_gene_joins_values_in_accession_order(tmp_path: Path) -> None:
    row = _annotated(tmp_path).filter(pl.col(_ID_COL) == "ENSG_B").row(0, named=True)
    # P00003 is missing from UniProt: its accession is listed, its text is not.
    assert row[GENE_INFO_UNIPROT_ID_COL] == GENE_INFO_SEPARATOR.join(
        ["P00002", "P00003"]
    )
    assert row[GENE_INFO_PROTEIN_NAMES_COL] == "Protein B"
    assert row[GENE_INFO_SUBCELLULAR_LOCALISATION_COL] is None


def test_non_coding_gene_has_symbols_but_no_protein_info(tmp_path: Path) -> None:
    row = _annotated(tmp_path).filter(pl.col(_ID_COL) == "ENSG_LNC").row(0, named=True)
    assert row[GENE_INFO_SYMBOL_COL] == GENE_INFO_SEPARATOR.join(["LINC1", "LINC2"])
    assert row[GENE_INFO_UNIPROT_ID_COL] is None


def test_gene_missing_from_hgnc_and_null_id_get_null_info(tmp_path: Path) -> None:
    annotated = _annotated(tmp_path)
    unannotated = annotated.filter(
        pl.col(_ID_COL).is_null() | (pl.col(_ID_COL) == "ENSG_GONE")
    )
    assert unannotated.height == 2
    assert unannotated.select(GENE_INFO_COLS).null_count().row(0) == (2,) * len(
        GENE_INFO_COLS
    )
