"""
Tests implemented by Claude
"""

from pathlib import Path

import pandas as pd
import pytest

from mecfs_bio.constants.msigdb_columns import (
    EXACT_SOURCE,
    GENE_SYMBOLS,
    NCBI_IDS,
    STANDARD_NAME,
    SYSTEMATIC_NAME,
)
from mecfs_bio.constants.vocabulary_classes.gene_set import MSigDBGeneSetSpec
from mecfs_bio.util.gene_set.msigdb_lookup import (
    JACCARD_COL_INDEX,
    JACCARD_COL_N_INTERSECTION,
    JACCARD_COL_N_UNION,
    JACCARD_COL_SET_1,
    JACCARD_COL_SET_2,
    gene_set_jaccard_index,
    verify_gene_set_specs,
)


def _write_parquet(tmp_path: Path, rows: list[dict]) -> Path:
    p = tmp_path / "gene_sets.parquet"
    pd.DataFrame(rows).to_parquet(p)
    return p


_ROWS = [
    {
        STANDARD_NAME: "HALLMARK_HYPOXIA",
        SYSTEMATIC_NAME: "M5891",
        EXACT_SOURCE: None,
        GENE_SYMBOLS: ["A", "B", "C"],
        NCBI_IDS: [1, 2, 3],
    },
    {
        STANDARD_NAME: "GOBP_RESPONSE_TO_VIRUS",
        SYSTEMATIC_NAME: "M16779",
        EXACT_SOURCE: "GO:0009615",
        GENE_SYMBOLS: ["B", "C", "D"],
        NCBI_IDS: [2, 3, 4],
    },
    {
        STANDARD_NAME: "KEGG_GLYCOLYSIS",
        SYSTEMATIC_NAME: "M123",
        EXACT_SOURCE: "hsa00010",
        GENE_SYMBOLS: ["X", "Y"],
        NCBI_IDS: [99, 100],
    },
    # duplicate standard_name to test ambiguity detection
    {
        STANDARD_NAME: "DUPLICATE_SET",
        SYSTEMATIC_NAME: "M001",
        EXACT_SOURCE: "SRC:001",
        GENE_SYMBOLS: [],
        NCBI_IDS: [],
    },
    {
        STANDARD_NAME: "DUPLICATE_SET",
        SYSTEMATIC_NAME: "M002",
        EXACT_SOURCE: "SRC:002",
        GENE_SYMBOLS: [],
        NCBI_IDS: [],
    },
]


def test_valid_specs_pass(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    specs = [
        MSigDBGeneSetSpec(
            standard_name="HALLMARK_HYPOXIA", exact_source=None, systematic_name="M5891"
        ),
        MSigDBGeneSetSpec(
            standard_name="GOBP_RESPONSE_TO_VIRUS",
            exact_source="GO:0009615",
            systematic_name="M16779",
        ),
    ]
    verify_gene_set_specs(path, specs)  # should not raise


def test_no_match_raises(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    specs = [MSigDBGeneSetSpec(standard_name="NONEXISTENT_SET", exact_source=None)]
    with pytest.raises(ValueError, match="no match"):
        verify_gene_set_specs(path, specs)


def test_wrong_exact_source_raises(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    specs = [
        MSigDBGeneSetSpec(
            standard_name="GOBP_RESPONSE_TO_VIRUS", exact_source="GO:WRONG"
        )
    ]
    with pytest.raises(ValueError, match="no match"):
        verify_gene_set_specs(path, specs)


def test_wrong_systematic_name_raises(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    specs = [
        MSigDBGeneSetSpec(
            standard_name="GOBP_RESPONSE_TO_VIRUS",
            exact_source="GO:0009615",
            systematic_name="MWRONG",
        )
    ]
    with pytest.raises(ValueError, match="no match"):
        verify_gene_set_specs(path, specs)


def test_ambiguous_standard_name_resolved_by_exact_source(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    # exact_source disambiguates
    specs = [MSigDBGeneSetSpec(standard_name="DUPLICATE_SET", exact_source="SRC:001")]
    verify_gene_set_specs(path, specs)  # should not raise


def test_ambiguous_match_raises(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    # no disambiguating fields — both rows match
    specs = [MSigDBGeneSetSpec(standard_name="DUPLICATE_SET", exact_source=None)]
    with pytest.raises(ValueError, match="2 matches"):
        verify_gene_set_specs(path, specs)


def test_multiple_failures_reported_together(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    specs = [
        MSigDBGeneSetSpec(standard_name="MISSING_ONE", exact_source=None),
        MSigDBGeneSetSpec(standard_name="MISSING_TWO", exact_source=None),
    ]
    with pytest.raises(ValueError) as exc_info:
        verify_gene_set_specs(path, specs)
    msg = str(exc_info.value)
    assert "MISSING_ONE" in msg
    assert "MISSING_TWO" in msg
    assert "2 spec(s)" in msg


def test_empty_spec_list_passes(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    verify_gene_set_specs(path, [])  # should not raise


# --- gene_set_jaccard_index ---

_SPEC_HYPOXIA = MSigDBGeneSetSpec(
    standard_name="HALLMARK_HYPOXIA", exact_source=None, systematic_name="M5891"
)
_SPEC_VIRUS = MSigDBGeneSetSpec(
    standard_name="GOBP_RESPONSE_TO_VIRUS", exact_source="GO:0009615"
)
_SPEC_GLYCOLYSIS = MSigDBGeneSetSpec(
    standard_name="KEGG_GLYCOLYSIS", exact_source="hsa00010"
)


def test_jaccard_correct_value(tmp_path: Path):
    # HYPOXIA={A,B,C}, VIRUS={B,C,D}  -> intersection=2, union=4, jaccard=0.5
    path = _write_parquet(tmp_path, _ROWS)
    df = gene_set_jaccard_index(path, [_SPEC_HYPOXIA, _SPEC_VIRUS])
    assert len(df) == 1
    row = df.iloc[0]
    assert row[JACCARD_COL_N_INTERSECTION] == 2
    assert row[JACCARD_COL_N_UNION] == 4
    assert row[JACCARD_COL_INDEX] == pytest.approx(0.5)


def test_jaccard_output_columns(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    df = gene_set_jaccard_index(path, [_SPEC_HYPOXIA, _SPEC_VIRUS])
    assert set(df.columns) == {
        JACCARD_COL_SET_1,
        JACCARD_COL_SET_2,
        JACCARD_COL_INDEX,
        JACCARD_COL_N_INTERSECTION,
        JACCARD_COL_N_UNION,
    }


def test_jaccard_sorted_descending(tmp_path: Path):
    # HYPOXIA/VIRUS share B,C (jaccard=0.5); HYPOXIA/GLYCOLYSIS share nothing (0.0)
    path = _write_parquet(tmp_path, _ROWS)
    df = gene_set_jaccard_index(path, [_SPEC_HYPOXIA, _SPEC_VIRUS, _SPEC_GLYCOLYSIS])
    assert df[JACCARD_COL_INDEX].is_monotonic_decreasing


def test_jaccard_n_pairs(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    specs = [_SPEC_HYPOXIA, _SPEC_VIRUS, _SPEC_GLYCOLYSIS]
    df = gene_set_jaccard_index(path, specs)
    assert len(df) == 3  # C(3,2) = 3


def test_jaccard_disjoint_sets_is_zero(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    df = gene_set_jaccard_index(path, [_SPEC_HYPOXIA, _SPEC_GLYCOLYSIS])
    assert df.iloc[0][JACCARD_COL_INDEX] == pytest.approx(0.0)
    assert df.iloc[0][JACCARD_COL_N_INTERSECTION] == 0


def test_jaccard_empty_spec_list_returns_empty_df(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    df = gene_set_jaccard_index(path, [])
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 0


def test_jaccard_single_spec_returns_empty_df(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    df = gene_set_jaccard_index(path, [_SPEC_HYPOXIA])
    assert len(df) == 0


def test_jaccard_bad_spec_raises(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    bad = MSigDBGeneSetSpec(standard_name="NONEXISTENT", exact_source=None)
    with pytest.raises(ValueError):
        gene_set_jaccard_index(path, [_SPEC_HYPOXIA, bad])
