"""
Utilities for verifying MSigDBGeneSetSpec references against a parquet table.
Implemented by Claude
"""

from itertools import combinations
from pathlib import Path

import pandas as pd
from attrs import frozen

from mecfs_bio.constants.msigdb_columns import (
    EXACT_SOURCE,
    NCBI_IDS,
    STANDARD_NAME,
    SYSTEMATIC_NAME,
)
from mecfs_bio.constants.vocabulary_classes.gene_set import MSigDBGeneSetSpec

# Output column names for gene_set_jaccard_index
JACCARD_COL_SET_1 = "standard_name_1"
JACCARD_COL_SET_2 = "standard_name_2"
JACCARD_COL_INDEX = "jaccard_index"
JACCARD_COL_N_INTERSECTION = "n_intersection"
JACCARD_COL_N_UNION = "n_union"


def apply_spec_mask(df: pd.DataFrame, spec: MSigDBGeneSetSpec) -> pd.Series:
    mask = df[STANDARD_NAME] == spec.standard_name
    if spec.exact_source is not None:
        mask &= df[EXACT_SOURCE] == spec.exact_source
    if spec.systematic_name is not None:
        mask &= df[SYSTEMATIC_NAME] == spec.systematic_name
    return mask


@frozen
class _SpecFailure:
    spec: MSigDBGeneSetSpec
    reason: str
    n_matches: int


def verify_gene_set_specs(
    parquet_path: Path,
    specs: list[MSigDBGeneSetSpec],
) -> None:
    """
    Verify that each spec in `specs` matches exactly one row in the MSigDB parquet table.

    Filters are applied for all non-None fields (standard_name, exact_source,
    systematic_name) combined with AND logic.  Raises ValueError listing every
    spec that matched zero rows or more than one row.
    """
    df = pd.read_parquet(
        parquet_path, columns=[STANDARD_NAME, SYSTEMATIC_NAME, EXACT_SOURCE]
    )

    failures: list[_SpecFailure] = []
    for spec in specs:
        n = int(apply_spec_mask(df, spec).sum())
        if n != 1:
            reason = "no match" if n == 0 else f"{n} matches (ambiguous)"
            failures.append(_SpecFailure(spec=spec, reason=reason, n_matches=n))

    if failures:
        lines = [f"verify_gene_set_specs failed for {len(failures)} spec(s):"]
        for f in failures:
            lines.append(
                f"  [{f.reason}] {STANDARD_NAME}={f.spec.standard_name!r}"
                f"  {EXACT_SOURCE}={f.spec.exact_source!r}"
                f"  {SYSTEMATIC_NAME}={f.spec.systematic_name!r}"
            )
        raise ValueError("\n".join(lines))


def gene_set_jaccard_index(
    parquet_path: Path,
    specs: list[MSigDBGeneSetSpec],
) -> pd.DataFrame:
    """
    Compute the Jaccard index of gene overlap for every pair of gene sets.

    Each spec must match exactly one row in the parquet (same rules as
    verify_gene_set_specs).  Returns a DataFrame with columns
    standard_name_1, standard_name_2, jaccard_index, n_intersection, n_union,
    sorted in descending order of jaccard_index.


    Purpose: evaluate whether a proposed gene set analysis contains redundant sets that can be reduced to
    increase statistical power without loss of granularity.
    """
    df = pd.read_parquet(
        parquet_path,
        columns=[STANDARD_NAME, SYSTEMATIC_NAME, EXACT_SOURCE, NCBI_IDS],
    )

    gene_sets: list[frozenset[int]] = []
    for spec in specs:
        rows = df[apply_spec_mask(df, spec)]
        if len(rows) != 1:
            raise ValueError(
                f"Expected exactly 1 match for {STANDARD_NAME}={spec.standard_name!r}, "
                f"got {len(rows)}"
            )
        gene_sets.append(frozenset(rows.iloc[0][NCBI_IDS]))

    records = []
    for (i, gs_i), (j, gs_j) in combinations(enumerate(gene_sets), 2):
        intersection = len(gs_i & gs_j)
        union = len(gs_i | gs_j)
        jaccard = intersection / union if union > 0 else 0.0
        records.append(
            {
                JACCARD_COL_SET_1: specs[i].standard_name,
                JACCARD_COL_SET_2: specs[j].standard_name,
                JACCARD_COL_INDEX: jaccard,
                JACCARD_COL_N_INTERSECTION: intersection,
                JACCARD_COL_N_UNION: union,
            }
        )

    result = pd.DataFrame(records)
    if not result.empty:
        result = result.sort_values(JACCARD_COL_INDEX, ascending=False).reset_index(
            drop=True
        )
    return result
