from pathlib import Path, PurePath

import numpy as np
import pandas as pd
import polars as pl
import pytest
from attrs import define, field

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.protein_coloc.coloc_protein_info_task import (
    ENSEMBL_DESCRIPTION_COL,
    ENSEMBL_ID_COL,
    GENE_SYMBOL_COL,
    INFO_SCHEMA,
    PROTEIN_NAMES_COL,
    SUBCELLULAR_LOCALISATION_COL,
    UNIPROT_DESCRIPTION_COL,
    ColocProteinInfoTask,
    GgetInfoLookup,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gene_thesaurus_constants import (
    THESAURUS_CHROM_COL,
    THESAURUS_ENSEMBL_ID_COL,
    THESAURUS_GENE_NAME_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    DISPLAY_GENE_COL,
    DISPLAY_LOCUS_COL,
    DISPLAY_PP_H4_COL,
    DISPLAY_PROTEIN_COL,
)


@define(slots=True)
class _FakeInfoLookup:
    calls: list[list[str]] = field(factory=list)

    def __call__(self, ensembl_ids: list[str]) -> pl.DataFrame:
        self.calls.append(list(ensembl_ids))
        return pl.DataFrame(
            {
                name: (
                    ensembl_ids
                    if name == ENSEMBL_ID_COL
                    else [f"{name} of {i}" for i in ensembl_ids]
                )
                for name in INFO_SCHEMA.names()
            },
            schema=INFO_SCHEMA,
        )


# (gene symbol, Ensembl id, chromosome). AAA also sits on an alternate-haplotype scaffold.
_THESAURUS_ROWS = [
    ("AAA", "ENSG_AAA", "1"),
    ("AAA", "ENSG_AAA_ALT", "HSCHR1_CTG1"),
    ("BBB", "ENSG_BBB", "2"),
    ("CCC", "ENSG_CCC", "X"),
    ("IL12A", "ENSG_IL12A", "3"),
    ("IL12B", "ENSG_IL12B", "5"),
    ("BAR", "ENSG_BAR", "4"),
]


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


def _write_thesaurus(path: Path, rows: list[tuple[str, str, str]]) -> None:
    pl.DataFrame(
        {
            THESAURUS_GENE_NAME_COL: [r[0] for r in rows],
            THESAURUS_ENSEMBL_ID_COL: [r[1] for r in rows],
            THESAURUS_CHROM_COL: [r[2] for r in rows],
        }
    ).write_csv(path)


def _info(
    tmp_path: Path,
    rows: list[tuple[str, str, str, float]],
    max_loci: int,
    info_lookup: _FakeInfoLookup,
    thesaurus_rows: list[tuple[str, str, str]] = _THESAURUS_ROWS,
) -> pl.DataFrame:
    display_path = tmp_path / "display.parquet"
    thesaurus_path = tmp_path / "thesaurus.csv"
    _write_thesaurus(thesaurus_path, thesaurus_rows)
    pl.DataFrame(
        {
            DISPLAY_GENE_COL: [r[0] for r in rows],
            DISPLAY_PROTEIN_COL: [r[1] for r in rows],
            DISPLAY_LOCUS_COL: [r[2] for r in rows],
            DISPLAY_PP_H4_COL: [r[3] for r in rows],
        },
        schema={
            DISPLAY_GENE_COL: pl.String,
            DISPLAY_PROTEIN_COL: pl.String,
            DISPLAY_LOCUS_COL: pl.String,
            DISPLAY_PP_H4_COL: pl.Float64,
        },
    ).write_parquet(display_path)
    task = ColocProteinInfoTask.create(
        asset_id="info",
        display_task=_table_task(
            "display", DataFrameReadSpec(DataFrameParquetFormat())
        ),
        gene_thesaurus_task=_table_task(
            "thesaurus", DataFrameReadSpec(DataFrameTextFormat(separator=","))
        ),
        max_loci=max_loci,
        info_lookup=info_lookup,
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(
        scratch,
        fetch=lambda asset_id: {
            AssetId("display"): FileAsset(display_path),
            AssetId("thesaurus"): FileAsset(thesaurus_path),
        }[asset_id],
        wf=make_wf(),
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
        info_lookup=_FakeInfoLookup(),
    )
    assert info[GENE_SYMBOL_COL].to_list() == ["AAA", "BBB"]
    assert info[ENSEMBL_ID_COL].to_list() == ["ENSG_AAA", "ENSG_BBB"]
    assert info[UNIPROT_DESCRIPTION_COL].to_list() == [
        "uniprot_description of ENSG_AAA",
        "uniprot_description of ENSG_BBB",
    ]


def test_multi_gene_complex_split(tmp_path: Path) -> None:
    info = _info(
        tmp_path,
        [("IL12A_IL12B", "OID9", "L1", 0.9)],
        max_loci=10,
        info_lookup=_FakeInfoLookup(),
    )
    assert info[GENE_SYMBOL_COL].to_list() == ["IL12A", "IL12B"]
    assert info[DISPLAY_PROTEIN_COL].to_list() == ["OID9", "OID9"]


def test_unresolved_symbol_keeps_row_with_null_annotations(tmp_path: Path) -> None:
    info_lookup = _FakeInfoLookup()
    info = _info(
        tmp_path,
        [("FOO", "OID1", "L1", 0.95), ("BAR", "OID2", "L1", 0.9)],
        max_loci=10,
        info_lookup=info_lookup,
    )
    foo = info.filter(pl.col(GENE_SYMBOL_COL) == "FOO")
    assert foo.height == 1
    assert foo[ENSEMBL_ID_COL].item() is None
    assert foo[UNIPROT_DESCRIPTION_COL].item() is None
    assert info_lookup.calls == [["ENSG_BAR"]]


def test_alternate_scaffold_copy_of_a_gene_is_ignored(tmp_path: Path) -> None:
    info = _info(
        tmp_path,
        [("AAA", "OID1", "L1", 0.95)],
        max_loci=10,
        info_lookup=_FakeInfoLookup(),
    )
    assert info[ENSEMBL_ID_COL].to_list() == ["ENSG_AAA"]


def test_symbol_with_two_primary_assembly_ids_fails(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        _info(
            tmp_path,
            [("AAA", "OID1", "L1", 0.95)],
            max_loci=10,
            info_lookup=_FakeInfoLookup(),
            thesaurus_rows=[("AAA", "ENSG_AAA", "1"), ("AAA", "ENSG_AAA2", "7")],
        )


def test_empty_display_gives_empty_table(tmp_path: Path) -> None:
    info_lookup = _FakeInfoLookup()
    info = _info(
        tmp_path,
        [],
        max_loci=10,
        info_lookup=info_lookup,
    )
    assert info.height == 0
    assert set(INFO_SCHEMA.names()) <= set(info.columns)
    assert info_lookup.calls == []


def test_gget_info_lookup_flattens_lists_and_blanks() -> None:
    def fake_fetch(ensembl_ids: list[str]) -> pd.DataFrame:
        return pd.DataFrame(
            {
                UNIPROT_DESCRIPTION_COL: ["", "Binds things."],
                SUBCELLULAR_LOCALISATION_COL: [np.nan, ["Secreted", "Cell membrane"]],
                PROTEIN_NAMES_COL: ["Protein A", "Protein B"],
            },
            index=ensembl_ids,
        )

    info = GgetInfoLookup(fetch=fake_fetch)(["ENSG1", "ENSG2"])
    assert info.schema == INFO_SCHEMA
    assert info[ENSEMBL_ID_COL].to_list() == ["ENSG1", "ENSG2"]
    assert info[UNIPROT_DESCRIPTION_COL].to_list() == [None, "Binds things."]
    assert info[SUBCELLULAR_LOCALISATION_COL].to_list() == [
        None,
        "Secreted; Cell membrane",
    ]
    assert info[PROTEIN_NAMES_COL].to_list() == ["Protein A", "Protein B"]
    assert info[ENSEMBL_DESCRIPTION_COL].to_list() == [None, None]
