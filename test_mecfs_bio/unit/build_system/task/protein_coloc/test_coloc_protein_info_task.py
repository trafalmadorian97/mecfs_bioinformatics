from pathlib import Path, PurePath

import polars as pl
import pytest

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.ppp_ldsc.gene_coords import (
    ST3_OLINK_ID_COL,
    ST3_UNIPROT_COL,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_protein_info_task import (
    ENSEMBL_ID_COL,
    GENE_NAME_COL,
    GENE_SYMBOL_COL,
    INFO_SCHEMA,
    PROTEIN_NAMES_COL,
    SUBCELLULAR_LOCALISATION_COL,
    UNIPROT_DESCRIPTION_COL,
    UNIPROT_ID_COL,
    ColocProteinInfoTask,
)
from mecfs_bio.build_system.task.protein_coloc.protein_accessions_reader import (
    St3ProteinAccessionsReader,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.hgnc_constants import (
    HGNC_ENSEMBL_ID_COL,
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

# (Olink id, ";"-separated UniProt accessions)
_PROTEIN_TABLE_ROWS = [
    ("OID1", "P00001"),
    ("OID2", "P00002"),
    ("OID3", "P00003"),
    ("OID4", "P23381"),
    ("OID5", "Q99999"),
    ("OID6", "P59665;P59665"),
    ("OID7", "Q9BTE6"),
    ("OID9", "Q14213;Q8NEV9"),
]

# (accession, protein names, function, subcellular location)
_UNIPROT_ROWS = [
    (
        "P00001",
        "Protein A",
        "FUNCTION: Binds things. {ECO:0000269|PubMed:1}.",
        "SUBCELLULAR LOCATION: Secreted {ECO:0000305}.",
    ),
    (
        "P00002",
        "Protein B",
        None,
        "SUBCELLULAR LOCATION: [Cleaved B]: Secreted. {ECO:0000269|PubMed:2}; "
        "SUBCELLULAR LOCATION: Cell membrane; Single-pass type I membrane protein.",
    ),
    ("P00003", "Protein C", None, None),
    ("P23381", "Tryptophan--tRNA ligase", "FUNCTION: Ligase.", None),
]

_PROTEIN_CODING = "gene with protein product"

# (symbol, name, locus type, Ensembl id, "|"-separated UniProt accessions)
_HGNC_ROWS = [
    ("AAA", "gene A", _PROTEIN_CODING, "ENSG_AAA", "P00001"),
    ("BBB", "gene B", _PROTEIN_CODING, "ENSG_BBB", "P00002"),
    ("CCC", "gene C", _PROTEIN_CODING, "ENSG_CCC", "P99998|P00003"),
    (
        "WARS1",
        "tryptophanyl-tRNA synthetase 1",
        _PROTEIN_CODING,
        "ENSG_WARS1",
        "P23381",
    ),
    ("DEFA1B", "defensin alpha 1B", _PROTEIN_CODING, "ENSG_DEFA1B", "P59665"),
    ("DEFA1", "defensin alpha 1", _PROTEIN_CODING, "ENSG_DEFA1", "P59665"),
    ("AARSD1", "alanyl-tRNA editing", _PROTEIN_CODING, "ENSG_AARSD1", "Q9BTE6"),
    (
        "PTGES3L-AARSD1",
        "readthrough",
        HGNC_READTHROUGH_LOCUS_TYPE,
        "ENSG_RT",
        "Q9BTE6",
    ),
    ("EBI3", "EBI3 subunit", _PROTEIN_CODING, "ENSG_EBI3", "Q14213"),
    ("IL27", "IL27 subunit", _PROTEIN_CODING, "ENSG_IL27", "Q8NEV9"),
]

_DISPLAY_SCHEMA = {
    DISPLAY_GENE_COL: pl.String,
    DISPLAY_PROTEIN_COL: pl.String,
    DISPLAY_LOCUS_COL: pl.String,
    DISPLAY_PP_H4_COL: pl.Float64,
}


def _table_task(asset_id: str, read_spec: DataFrameReadSpec) -> FakeTask:
    return FakeTask(
        ResultTableMeta(
            id=AssetId(asset_id),
            trait="synthetic",
            project="test",
            sub_dir=PurePath("analysis"),
            extension=".parquet",
            read_spec=read_spec,
        )
    )


def _write_inputs(
    tmp_path: Path, display_rows: list[tuple[str, str, str, float]]
) -> dict[AssetId, FileAsset]:
    paths = {
        "display": tmp_path / "display.parquet",
        "protein_table": tmp_path / "protein_table.parquet",
        "uniprot": tmp_path / "uniprot.parquet",
        "hgnc": tmp_path / "hgnc.tsv",
    }
    pl.DataFrame(
        [list(r) for r in display_rows], schema=_DISPLAY_SCHEMA, orient="row"
    ).write_parquet(paths["display"])
    pl.DataFrame(
        _PROTEIN_TABLE_ROWS, schema=[ST3_OLINK_ID_COL, ST3_UNIPROT_COL], orient="row"
    ).write_parquet(paths["protein_table"])
    pl.DataFrame(
        _UNIPROT_ROWS,
        schema=[
            UNIPROT_DAT_ID_COL,
            UNIPROT_DAT_PROTEIN_NAMES_COL,
            UNIPROT_DAT_FUNCTION_COL,
            UNIPROT_DAT_SUBCELLULAR_LOC_COL,
        ],
        orient="row",
    ).write_parquet(paths["uniprot"])
    pl.DataFrame(
        _HGNC_ROWS,
        schema=[
            HGNC_SYMBOL_COL,
            HGNC_NAME_COL,
            HGNC_LOCUS_TYPE_COL,
            HGNC_ENSEMBL_ID_COL,
            HGNC_UNIPROT_IDS_COL,
        ],
        orient="row",
    ).write_csv(paths["hgnc"], separator="\t")
    return {AssetId(name): FileAsset(path) for name, path in paths.items()}


def _info(
    tmp_path: Path,
    display_rows: list[tuple[str, str, str, float]],
    max_loci: int = 10,
) -> pl.DataFrame:
    assets = _write_inputs(tmp_path, display_rows)
    parquet = DataFrameReadSpec(DataFrameParquetFormat())
    task = ColocProteinInfoTask.create(
        asset_id="info",
        display_task=_table_task("display", parquet),
        protein_table_task=_table_task("protein_table", parquet),
        protein_accessions_reader=St3ProteinAccessionsReader(),
        uniprot_task=_table_task("uniprot", parquet),
        hgnc_task=_table_task(
            "hgnc", DataFrameReadSpec(DataFrameTextFormat(separator="\t"))
        ),
        max_loci=max_loci,
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(
        scratch, fetch=lambda asset_id: assets[asset_id], wf=make_wf()
    )
    assert isinstance(result, FileAsset)
    return pl.read_parquet(result.path)


def test_one_row_per_gene_of_top_loci(tmp_path: Path) -> None:
    info = _info(
        tmp_path,
        [
            ("AAA", "OID1", "L1", 0.95),
            ("BBB", "OID2", "L1", 0.85),
            ("CCC", "OID3", "L2", 0.9),
        ],
        max_loci=1,
    )
    assert info.select(
        DISPLAY_PROTEIN_COL, UNIPROT_ID_COL, GENE_SYMBOL_COL, ENSEMBL_ID_COL
    ).rows() == [
        ("OID1", "P00001", "AAA", "ENSG_AAA"),
        ("OID2", "P00002", "BBB", "ENSG_BBB"),
    ]
    assert info[GENE_NAME_COL].to_list() == ["gene A", "gene B"]
    assert info[PROTEIN_NAMES_COL].to_list() == ["Protein A", "Protein B"]


def test_uniprot_text_is_cleaned(tmp_path: Path) -> None:
    """Field labels and evidence codes are dropped."""
    info = _info(tmp_path, [("AAA", "OID1", "L1", 0.95)])
    assert info[UNIPROT_DESCRIPTION_COL].to_list() == ["Binds things."]
    assert info[SUBCELLULAR_LOCALISATION_COL].to_list() == ["Secreted."]


def test_multiple_uniprot_comments_join_as_sentences(tmp_path: Path) -> None:
    info = _info(tmp_path, [("BBB", "OID2", "L1", 0.95)])
    assert info[SUBCELLULAR_LOCALISATION_COL].to_list() == [
        "[Cleaved B]: Secreted. Cell membrane; Single-pass type I membrane protein."
    ]


def test_gene_found_among_several_hgnc_accessions(tmp_path: Path) -> None:
    info = _info(tmp_path, [("CCC", "OID3", "L1", 0.9)])
    assert info[GENE_SYMBOL_COL].to_list() == ["CCC"]


def test_symbol_comes_from_hgnc_not_the_display_table(tmp_path: Path) -> None:
    info = _info(tmp_path, [("WARS", "OID4", "L1", 0.9)])
    assert info.select(GENE_SYMBOL_COL, ENSEMBL_ID_COL).rows() == [
        ("WARS1", "ENSG_WARS1")
    ]
    assert info[UNIPROT_DESCRIPTION_COL].to_list() == ["Ligase."]


def test_complex_gives_one_row_per_subunit_gene(tmp_path: Path) -> None:
    info = _info(tmp_path, [("EBI3_IL27", "OID9", "L1", 0.9)])
    assert info.select(UNIPROT_ID_COL, GENE_SYMBOL_COL, ENSEMBL_ID_COL).rows() == [
        ("Q14213", "EBI3", "ENSG_EBI3"),
        ("Q8NEV9", "IL27", "ENSG_IL27"),
    ]
    assert info[DISPLAY_PROTEIN_COL].to_list() == ["OID9", "OID9"]


def test_duplicated_genes_of_one_protein_each_get_a_row(tmp_path: Path) -> None:
    info = _info(tmp_path, [("DEFA1_DEFA1B", "OID6", "L1", 0.9)])
    assert info[GENE_SYMBOL_COL].to_list() == ["DEFA1", "DEFA1B"]


def test_readthrough_gene_is_dropped(tmp_path: Path) -> None:
    info = _info(tmp_path, [("AARSD1", "OID7", "L1", 0.9)])
    assert info[GENE_SYMBOL_COL].to_list() == ["AARSD1"]


def test_accession_missing_from_references_keeps_row(tmp_path: Path) -> None:
    info = _info(tmp_path, [("OBSOLETE", "OID5", "L1", 0.9)])
    assert info[UNIPROT_ID_COL].to_list() == ["Q99999"]
    assert info.select(GENE_SYMBOL_COL, ENSEMBL_ID_COL, PROTEIN_NAMES_COL).rows() == [
        (None, None, None)
    ]


def test_protein_missing_from_protein_table_fails(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        _info(tmp_path, [("ZZZ", "OID404", "L1", 0.95)])


def test_empty_display_gives_empty_table(tmp_path: Path) -> None:
    info = _info(tmp_path, [])
    assert info.height == 0
    assert set(INFO_SCHEMA.names()) <= set(info.columns)
