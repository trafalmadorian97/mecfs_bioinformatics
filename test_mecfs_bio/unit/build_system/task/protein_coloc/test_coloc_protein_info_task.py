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
    ST3_GENE_SYMBOL_COL,
    ST3_OLINK_ID_COL,
    ST3_UNIPROT_COL,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_protein_info_task import (
    ENSEMBL_DESCRIPTION_COL,
    ENSEMBL_ID_COL,
    GENE_SYMBOL_COL,
    INFO_SCHEMA,
    PRIMARY_GENE_NAME_COL,
    PROTEIN_NAMES_COL,
    SUBCELLULAR_LOCALISATION_COL,
    UNIPROT_DESCRIPTION_COL,
    UNIPROT_ID_COL,
    ColocProteinInfoTask,
)
from mecfs_bio.build_system.task.protein_coloc.protein_genes_reader import (
    St3ProteinGenesReader,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gene_thesaurus_constants import (
    THESAURUS_CHROM_COL,
    THESAURUS_DESCRIPTION_COL,
    THESAURUS_ENSEMBL_ID_COL,
    THESAURUS_GENE_NAME_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    DISPLAY_GENE_COL,
    DISPLAY_LOCUS_COL,
    DISPLAY_PP_H4_COL,
    DISPLAY_PROTEIN_COL,
)
from mecfs_bio.constants.uniprot_constants import (
    UNIPROT_DAT_FUNCTION_COL,
    UNIPROT_DAT_GENE_NAMES_COL,
    UNIPROT_DAT_ID_COL,
    UNIPROT_DAT_PROTEIN_NAMES_COL,
    UNIPROT_DAT_SUBCELLULAR_LOC_COL,
)

# (Olink id, ";"-separated gene symbols, ";"-separated UniProt accessions)
_PROTEIN_TABLE_ROWS = [
    ("OID1", "AAA", "P00001"),
    ("OID2", "BBB", "P00002"),
    ("OID3", "CCC", "P00003"),
    ("OID4", "WARS", "P23381"),
    ("OID5", "OBSOLETE", "Q99999"),
    ("OID9", "EBI3;IL27", "Q14213;Q8NEV9"),
]

# (accession, gene names, protein names, function, subcellular location)
_UNIPROT_ROWS = [
    (
        "P00001",
        "AAA AAA_ALIAS",
        "Protein A",
        "FUNCTION: Binds things. {ECO:0000269|PubMed:1}.",
        "SUBCELLULAR LOCATION: Secreted {ECO:0000305}.",
    ),
    (
        "P00002",
        "BBB",
        "Protein B",
        None,
        "SUBCELLULAR LOCATION: [Cleaved B]: Secreted. {ECO:0000269|PubMed:2}; "
        "SUBCELLULAR LOCATION: Cell membrane; Single-pass type I membrane protein.",
    ),
    ("P00003", "CCC", "Protein C", None, None),
    (
        "P23381",
        "WARS1 IFI53 WARS",
        "Tryptophan--tRNA ligase",
        "FUNCTION: Ligase.",
        None,
    ),
    ("Q14213", "EBI3 IL27B", "Interleukin-27 subunit beta", None, None),
    ("Q8NEV9", "IL27 IL27A", "Interleukin-27 subunit alpha", None, None),
]

# (gene symbol, Ensembl id, chromosome, description). AAA also sits on an alternate-haplotype
# scaffold.
_THESAURUS_ROWS = [
    ("AAA", "ENSG_AAA", "1", "gene A [Source:HGNC Symbol;Acc:HGNC:1]"),
    ("AAA", "ENSG_AAA_ALT", "HSCHR1_CTG1", "gene A"),
    ("BBB", "ENSG_BBB", "2", "gene B"),
    ("CCC", "ENSG_CCC", "X", "gene C"),
    ("EBI3", "ENSG_EBI3", "19", "EBI3"),
    ("IL27", "ENSG_IL27", "16", "IL27"),
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
    tmp_path: Path,
    display_rows: list[tuple[str, str, str, float]],
    thesaurus_rows: list[tuple[str, str, str, str]],
) -> dict[AssetId, FileAsset]:
    paths = {
        "display": tmp_path / "display.parquet",
        "protein_table": tmp_path / "protein_table.parquet",
        "uniprot": tmp_path / "uniprot.parquet",
        "thesaurus": tmp_path / "thesaurus.csv",
    }
    pl.DataFrame(
        [list(r) for r in display_rows], schema=_DISPLAY_SCHEMA, orient="row"
    ).write_parquet(paths["display"])
    pl.DataFrame(
        _PROTEIN_TABLE_ROWS,
        schema=[ST3_OLINK_ID_COL, ST3_GENE_SYMBOL_COL, ST3_UNIPROT_COL],
        orient="row",
    ).write_parquet(paths["protein_table"])
    pl.DataFrame(
        _UNIPROT_ROWS,
        schema=[
            UNIPROT_DAT_ID_COL,
            UNIPROT_DAT_GENE_NAMES_COL,
            UNIPROT_DAT_PROTEIN_NAMES_COL,
            UNIPROT_DAT_FUNCTION_COL,
            UNIPROT_DAT_SUBCELLULAR_LOC_COL,
        ],
        orient="row",
    ).write_parquet(paths["uniprot"])
    pl.DataFrame(
        thesaurus_rows,
        schema=[
            THESAURUS_GENE_NAME_COL,
            THESAURUS_ENSEMBL_ID_COL,
            THESAURUS_CHROM_COL,
            THESAURUS_DESCRIPTION_COL,
        ],
        orient="row",
    ).write_csv(paths["thesaurus"])
    return {AssetId(name): FileAsset(path) for name, path in paths.items()}


def _info(
    tmp_path: Path,
    display_rows: list[tuple[str, str, str, float]],
    max_loci: int = 10,
    thesaurus_rows: list[tuple[str, str, str, str]] = _THESAURUS_ROWS,
) -> pl.DataFrame:
    assets = _write_inputs(tmp_path, display_rows, thesaurus_rows)
    parquet = DataFrameReadSpec(DataFrameParquetFormat())
    task = ColocProteinInfoTask.create(
        asset_id="info",
        display_task=_table_task("display", parquet),
        protein_table_task=_table_task("protein_table", parquet),
        protein_genes_reader=St3ProteinGenesReader(),
        uniprot_task=_table_task("uniprot", parquet),
        gene_thesaurus_task=_table_task(
            "thesaurus", DataFrameReadSpec(DataFrameTextFormat(separator=","))
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
    assert info[GENE_SYMBOL_COL].to_list() == ["AAA", "BBB"]
    assert info[UNIPROT_ID_COL].to_list() == ["P00001", "P00002"]
    assert info[ENSEMBL_ID_COL].to_list() == ["ENSG_AAA", "ENSG_BBB"]
    assert info[PRIMARY_GENE_NAME_COL].to_list() == ["AAA", "BBB"]
    assert info[PROTEIN_NAMES_COL].to_list() == ["Protein A", "Protein B"]


def test_annotation_text_is_cleaned(tmp_path: Path) -> None:
    """Field labels, evidence codes and the thesaurus source tag are dropped."""
    info = _info(tmp_path, [("AAA", "OID1", "L1", 0.95)])
    assert info[UNIPROT_DESCRIPTION_COL].to_list() == ["Binds things."]
    assert info[SUBCELLULAR_LOCALISATION_COL].to_list() == ["Secreted."]
    assert info[ENSEMBL_DESCRIPTION_COL].to_list() == ["gene A"]


def test_multiple_uniprot_comments_join_as_sentences(tmp_path: Path) -> None:
    info = _info(tmp_path, [("BBB", "OID2", "L1", 0.95)])
    assert info[SUBCELLULAR_LOCALISATION_COL].to_list() == [
        "[Cleaved B]: Secreted. Cell membrane; Single-pass type I membrane protein."
    ]


def test_complex_gives_one_row_per_gene_with_its_accession(tmp_path: Path) -> None:
    info = _info(tmp_path, [("EBI3_IL27", "OID9", "L1", 0.9)])
    assert info.select(GENE_SYMBOL_COL, UNIPROT_ID_COL, ENSEMBL_ID_COL).rows() == [
        ("EBI3", "Q14213", "ENSG_EBI3"),
        ("IL27", "Q8NEV9", "ENSG_IL27"),
    ]
    assert info[DISPLAY_PROTEIN_COL].to_list() == ["OID9", "OID9"]


def test_outdated_symbol_keeps_uniprot_annotations(tmp_path: Path) -> None:
    info = _info(tmp_path, [("WARS", "OID4", "L1", 0.9)])
    assert info[ENSEMBL_ID_COL].to_list() == [None]
    assert info[PRIMARY_GENE_NAME_COL].to_list() == ["WARS1"]
    assert info[UNIPROT_DESCRIPTION_COL].to_list() == ["Ligase."]


def test_accession_missing_from_uniprot_keeps_row(tmp_path: Path) -> None:
    info = _info(tmp_path, [("OBSOLETE", "OID5", "L1", 0.9)])
    assert info[UNIPROT_ID_COL].to_list() == ["Q99999"]
    assert info[PROTEIN_NAMES_COL].to_list() == [None]


def test_alternate_scaffold_copy_of_a_gene_is_ignored(tmp_path: Path) -> None:
    info = _info(tmp_path, [("AAA", "OID1", "L1", 0.95)])
    assert info[ENSEMBL_ID_COL].to_list() == ["ENSG_AAA"]


def test_symbol_with_two_primary_assembly_ids_fails(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        _info(
            tmp_path,
            [("AAA", "OID1", "L1", 0.95)],
            thesaurus_rows=[
                ("AAA", "ENSG_AAA", "1", "gene A"),
                ("AAA", "ENSG_AAA2", "7", "gene A"),
            ],
        )


def test_protein_missing_from_protein_table_fails(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        _info(tmp_path, [("ZZZ", "OID404", "L1", 0.95)])


def test_empty_display_gives_empty_table(tmp_path: Path) -> None:
    info = _info(tmp_path, [])
    assert info.height == 0
    assert set(INFO_SCHEMA.names()) <= set(info.columns)
