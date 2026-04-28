"""
Greedy diverse subset selection for MSigDB gene sets.

Added by Claude
"""

from pathlib import Path
from typing import Literal

import pandas as pd

from mecfs_bio.constants.msigdb_columns import (
    EXACT_SOURCE,
    NCBI_IDS,
    STANDARD_NAME,
    SYSTEMATIC_NAME,
)
from mecfs_bio.constants.vocabulary_classes.gene_set import MSigDBGeneSetSpec
from mecfs_bio.util.gene_set.msigdb_lookup import apply_spec_mask

DiversityStrategy = Literal["min_max_jaccard", "min_mean_jaccard"]


def select_diverse_subset(
    gene_sets: list[MSigDBGeneSetSpec],
    db_path: Path,
    initial_gene_sets: list[MSigDBGeneSetSpec],
    target_size: int,
    strategy: DiversityStrategy = "min_max_jaccard",
) -> list[MSigDBGeneSetSpec]:
    """
    Greedily select a diverse subset of `target_size` gene sets from `gene_sets`.

    Starts from `initial_gene_sets` and repeatedly adds the candidate that is
    most dissimilar to the current selection, using one of two strategies:

    - ``"min_max_jaccard"`` (default): at each step pick the candidate whose
      maximum Jaccard similarity to any already-selected set is smallest.
      This is the greedy k-center / farthest-point algorithm; it guarantees
      no two selected sets are more similar than necessary and has a known
      2-approximation bound for the k-center objective.

    - ``"min_mean_jaccard"``: at each step pick the candidate with the lowest
      mean Jaccard similarity to all already-selected sets.  This can favour
      a set that overlaps heavily with one existing set if it is dissimilar to
      all others, so it is more permissive about local clustering.

    ``initial_gene_sets`` must be a subset of ``gene_sets``.  When it is empty
    the first element of ``gene_sets`` is used as the seed.

    Raises ``ValueError`` if ``target_size`` exceeds ``len(gene_sets)``, or if
    ``target_size`` is less than ``len(initial_gene_sets)``, or if any element
    of ``initial_gene_sets`` is not present in ``gene_sets``.
    """
    if target_size > len(gene_sets):
        raise ValueError(
            f"target_size ({target_size}) exceeds the number of available "
            f"gene sets ({len(gene_sets)})"
        )

    gene_set_names = {spec.standard_name for spec in gene_sets}
    for spec in initial_gene_sets:
        if spec.standard_name not in gene_set_names:
            raise ValueError(
                f"initial_gene_sets contains {spec.standard_name!r} "
                f"which is not present in gene_sets"
            )

    if target_size < len(initial_gene_sets):
        raise ValueError(
            f"target_size ({target_size}) is less than the number of "
            f"initial_gene_sets ({len(initial_gene_sets)})"
        )

    df = pd.read_parquet(
        db_path,
        columns=[STANDARD_NAME, SYSTEMATIC_NAME, EXACT_SOURCE, NCBI_IDS],
    )

    ncbi_ids: dict[str, frozenset[int]] = {}
    for spec in gene_sets:
        rows = df[apply_spec_mask(df, spec)]
        if len(rows) != 1:
            raise ValueError(
                f"Expected exactly 1 match for {STANDARD_NAME}={spec.standard_name!r}, "
                f"got {len(rows)}"
            )
        ncbi_ids[spec.standard_name] = frozenset(rows.iloc[0][NCBI_IDS])

    def jaccard(a: frozenset, b: frozenset) -> float:
        union = len(a | b)
        return len(a & b) / union if union > 0 else 0.0

    selected: list[MSigDBGeneSetSpec] = list(initial_gene_sets)
    selected_names: set[str] = {s.standard_name for s in selected}
    candidates: list[MSigDBGeneSetSpec] = [
        s for s in gene_sets if s.standard_name not in selected_names
    ]

    if not selected and candidates:
        seed = candidates.pop(0)
        selected.append(seed)
        selected_names.add(seed.standard_name)

    while len(selected) < target_size:
        best_spec: MSigDBGeneSetSpec | None = None
        best_score = float("inf")

        for candidate in candidates:
            cand_ids = ncbi_ids[candidate.standard_name]
            similarities = [
                jaccard(cand_ids, ncbi_ids[s.standard_name]) for s in selected
            ]
            if strategy == "min_max_jaccard":
                score = max(similarities)
            else:
                score = sum(similarities) / len(similarities)

            if score < best_score:
                best_score = score
                best_spec = candidate

        assert best_spec is not None
        selected.append(best_spec)
        selected_names.add(best_spec.standard_name)
        candidates = [
            c for c in candidates if c.standard_name != best_spec.standard_name
        ]

    return selected
