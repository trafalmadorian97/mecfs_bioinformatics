from pathlib import Path

import pandas as pd
import pytest

from mecfs_bio.constants.msigdb_columns import (
    EXACT_SOURCE,
    NCBI_IDS,
    STANDARD_NAME,
    SYSTEMATIC_NAME,
)
from mecfs_bio.constants.vocabulary_classes.gene_set import MSigDBGeneSetSpec
from mecfs_bio.util.gene_set.msigdb_diversity import select_diverse_subset


def _spec(name: str) -> MSigDBGeneSetSpec:
    return MSigDBGeneSetSpec(standard_name=name, exact_source=None)


def _write_parquet(tmp_path: Path, rows: list[dict]) -> Path:
    p = tmp_path / "gene_sets.parquet"
    pd.DataFrame(rows).to_parquet(p)
    return p


# ---------------------------------------------------------------------------
# Fixture: four disjoint sets plus two with partial overlap
#
#   A = {1,2,3}      B = {4,5,6}      C = {7,8,9}      D = {10,11,12}
#   X = {1,4,7,10,13}        (one gene from each of A,B,C,D → uniform low Jaccard)
#   Y = {1,2,3,13,14,15,16,17,18,19,20}  (heavily overlaps A, disjoint from B/C/D)
# ---------------------------------------------------------------------------
_A = _spec("SET_A")
_B = _spec("SET_B")
_C = _spec("SET_C")
_D = _spec("SET_D")
_X = _spec("SET_X")
_Y = _spec("SET_Y")

_ROWS = [
    {
        STANDARD_NAME: "SET_A",
        SYSTEMATIC_NAME: "M1",
        EXACT_SOURCE: None,
        NCBI_IDS: [1, 2, 3],
    },
    {
        STANDARD_NAME: "SET_B",
        SYSTEMATIC_NAME: "M2",
        EXACT_SOURCE: None,
        NCBI_IDS: [4, 5, 6],
    },
    {
        STANDARD_NAME: "SET_C",
        SYSTEMATIC_NAME: "M3",
        EXACT_SOURCE: None,
        NCBI_IDS: [7, 8, 9],
    },
    {
        STANDARD_NAME: "SET_D",
        SYSTEMATIC_NAME: "M4",
        EXACT_SOURCE: None,
        NCBI_IDS: [10, 11, 12],
    },
    # X: one gene from each of A/B/C/D — uniformly low Jaccard with all
    {
        STANDARD_NAME: "SET_X",
        SYSTEMATIC_NAME: "M5",
        EXACT_SOURCE: None,
        NCBI_IDS: [1, 4, 7, 10, 13],
    },
    # Y: heavy overlap with A only — low MEAN but high MAX Jaccard
    {
        STANDARD_NAME: "SET_Y",
        SYSTEMATIC_NAME: "M6",
        EXACT_SOURCE: None,
        NCBI_IDS: [1, 2, 3, 13, 14, 15, 16, 17, 18, 19, 20],
    },
]


def test_returns_target_size(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    result = select_diverse_subset([_A, _B, _C, _D], path, [], target_size=3)
    assert len(result) == 3


def test_initial_gene_sets_always_included(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    result = select_diverse_subset([_A, _B, _C, _D], path, [_A, _B], target_size=3)
    names = [s.standard_name for s in result]
    assert "SET_A" in names
    assert "SET_B" in names


def test_initial_gene_sets_appear_first(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    result = select_diverse_subset([_A, _B, _C, _D], path, [_A, _B], target_size=3)
    assert result[0].standard_name == "SET_A"
    assert result[1].standard_name == "SET_B"


def test_no_duplicates_in_result(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    result = select_diverse_subset([_A, _B, _C, _D, _X], path, [], target_size=4)
    names = [s.standard_name for s in result]
    assert len(names) == len(set(names))


def test_result_is_subset_of_gene_sets(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    all_specs = [_A, _B, _C, _D, _X, _Y]
    result = select_diverse_subset(all_specs, path, [], target_size=3)
    all_names = {s.standard_name for s in all_specs}
    for s in result:
        assert s.standard_name in all_names


def test_target_size_equals_input_returns_all(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    specs = [_A, _B, _C, _D]
    result = select_diverse_subset(specs, path, [], target_size=4)
    assert {s.standard_name for s in result} == {"SET_A", "SET_B", "SET_C", "SET_D"}


def test_target_size_exceeds_available_raises(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    with pytest.raises(ValueError, match="target_size"):
        select_diverse_subset([_A, _B], path, [], target_size=5)


def test_target_size_less_than_initial_raises(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    with pytest.raises(ValueError, match="target_size"):
        select_diverse_subset([_A, _B, _C, _D], path, [_A, _B, _C], target_size=2)


def test_initial_not_in_gene_sets_raises(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    outsider = _spec("SET_OUTSIDER")
    with pytest.raises(ValueError, match="not present in gene_sets"):
        select_diverse_subset([_A, _B, _C], path, [outsider], target_size=2)


def test_strategies_differ_on_designed_example(tmp_path: Path):
    """
    With selected={A,B,C,D} and candidates={X,Y}:

      X = {1,4,7,10,13}: Jaccard 1/7 with each of A,B,C,D → max=1/7, mean=1/7
      Y = {1..3,13..20}:  Jaccard 3/11 with A, 0 with B/C/D  → max=3/11, mean=3/44

    min_max_jaccard prefers X  (1/7 ≈ 0.143 < 3/11 ≈ 0.273)
    min_mean_jaccard prefers Y (3/44 ≈ 0.068 < 1/7 ≈ 0.143)
    """
    path = _write_parquet(tmp_path, _ROWS)
    base = [_A, _B, _C, _D]
    all_specs = base + [_X, _Y]

    result_max = select_diverse_subset(
        all_specs, path, base, target_size=5, strategy="min_max_jaccard"
    )
    result_mean = select_diverse_subset(
        all_specs, path, base, target_size=5, strategy="min_mean_jaccard"
    )

    fifth_max = result_max[4].standard_name
    fifth_mean = result_mean[4].standard_name

    assert fifth_max == "SET_X", f"min_max_jaccard should pick X first, got {fifth_max}"
    assert fifth_mean == "SET_Y", (
        f"min_mean_jaccard should pick Y first, got {fifth_mean}"
    )


def test_empty_initial_seeds_from_first_gene_set(tmp_path: Path):
    path = _write_parquet(tmp_path, _ROWS)
    result = select_diverse_subset([_A, _B, _C, _D], path, [], target_size=2)
    assert result[0].standard_name == "SET_A"
    assert len(result) == 2
