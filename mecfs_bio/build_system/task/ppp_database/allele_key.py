"""
Order-agnostic allele key, shared across the PPP-database tasks for allele-aware
variant matching (so a swapped effect/non-effect orientation still matches).

The key is valid only for SNVs. Mirrored indels such as T/TCA and TCA/T are distinct
variants, yet they sort to the same {A, B} key and would collide. Callers that key on it
must guarantee SNV-only input; assert_all_snv is the guard for that.

Remaining unordered_allele_key uses, all SNV-safe by assertion or by construction:
  - construct_ppp_variant_index_task: SNV-asserted (the HapMap3 index is SNV-only).
  - common_1kg_membership_task, build_slim_protein_parquet_task: the common-1kg PPP mode
    contains indels and is NOT built today; add filter-then-assert before enabling it.
HarmonizeGWASWithReferenceViaAlleles no longer uses this key but is likewise SNV-only and
asserts it. The annotation path and the SUSIE prior/main joins now join on the exact
(CHR, POS, EA, NEA) tuple instead of this key.
"""

import polars as pl


def unordered_allele_key(a: str, b: str) -> pl.Expr:
    """Sort the two allele columns and join them, so {A, B} == {B, A}.

    Valid only for SNVs: mirrored indels (T/TCA vs TCA/T) sort to the same key though they
    are distinct variants. Guard callers with assert_all_snv."""
    return (
        pl.when(pl.col(a) <= pl.col(b))
        .then(pl.col(a) + pl.lit("_") + pl.col(b))
        .otherwise(pl.col(b) + pl.lit("_") + pl.col(a))
    )


def assert_all_snv(df: pl.DataFrame, *allele_cols: str) -> None:
    """Fail fast if any allele in the named columns is not a single base.

    unordered_allele_key and the PPP variant index are only valid for SNVs: mirrored
    indels (T/TCA vs TCA/T) sort to the same key though they are distinct variants."""
    assert allele_cols, "assert_all_snv requires at least one allele column"
    non_snv = df.filter(
        pl.any_horizontal(pl.col(c).str.len_bytes() != 1 for c in allele_cols)
    )
    assert non_snv.height == 0, (
        f"{non_snv.height} non-SNV variant(s) in columns {allele_cols}; the unordered "
        f"allele key is invalid for indels. First rows:\n{non_snv.select(allele_cols).head()}"
    )
