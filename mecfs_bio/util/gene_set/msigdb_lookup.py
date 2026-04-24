"""
Utilities for verifying MSigDBGeneSetSpec references against a parquet table.
"""

from pathlib import Path

import pandas as pd
from attrs import frozen

from mecfs_bio.constants.vocabulary_classes.gene_set import MSigDBGeneSetSpec


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
        parquet_path, columns=["standard_name", "systematic_name", "exact_source"]
    )

    failures: list[_SpecFailure] = []
    for spec in specs:
        mask = df["standard_name"] == spec.standard_name
        if spec.exact_source is not None:
            mask &= df["exact_source"] == spec.exact_source
        if spec.systematic_name is not None:
            mask &= df["systematic_name"] == spec.systematic_name

        n = int(mask.sum())
        if n != 1:
            reason = "no match" if n == 0 else f"{n} matches (ambiguous)"
            failures.append(_SpecFailure(spec=spec, reason=reason, n_matches=n))

    if failures:
        lines = [f"verify_gene_set_specs failed for {len(failures)} spec(s):"]
        for f in failures:
            lines.append(
                f"  [{f.reason}] standard_name={f.spec.standard_name!r}"
                f"  exact_source={f.spec.exact_source!r}"
                f"  systematic_name={f.spec.systematic_name!r}"
            )
        raise ValueError("\n".join(lines))
