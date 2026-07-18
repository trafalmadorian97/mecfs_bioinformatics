"""
Align an external trait's signed z-score onto the shared UKB-PPP LDSC context SNP set.

The batched rg kernel needs the trait's z-score in the SAME row order and the SAME effect-allele
orientation as the per-protein z-scores. Protein betas are already oriented to the variant index's
effect allele (see align_protein_to_index), so we orient the trait to the index effect allele too:
join the trait to the context variants on rsID, keep z = BETA/SE where the trait's alleles match
the index orientation, negate it where they are swapped, and drop it (NaN) where the alleles do not
match at all or the trait lacks the variant. Strand-ambiguous variants are already removed from the
context, so no palindrome handling is needed here.

The trait sample size is taken per-SNP from the trait's N column when present, otherwise from a
caller-supplied constant. A collapsed trait/context overlap almost always means a harmonization
failure upstream (wrong genome build, malformed rsIDs, or systematic allele mismatch), so we fail
fast when fewer than min_trait_snps context variants match.
"""

from __future__ import annotations

import numpy as np
import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_rsid import (
    match_flipped_reference_expr,
    match_reference_expr,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

_CTX_ROW = "__ctx_row__"
_REF_EA = "__ref_ea__"
_REF_NEA = "__ref_nea__"
_TRAIT_Z = "__trait_z__"
_TRAIT_N = "__trait_n__"


@frozen
class TraitAligned:
    """The trait aligned to the context SNP set (parallel (S,) arrays in context row order).

    z: signed trait z-score oriented to the index effect allele, NaN where the trait lacks the
    variant or its alleles do not match. n: trait sample size per SNP, NaN where z is NaN.
    """

    z: np.ndarray
    n: np.ndarray

    def __attrs_post_init__(self) -> None:
        s = self.z.shape[0]
        for name, arr in (("z", self.z), ("n", self.n)):
            assert arr.ndim == 1, f"{name} must be 1-D, got shape {arr.shape}"
            assert arr.shape[0] == s, f"{name} has length {arr.shape[0]}, expected {s}"
            assert arr.dtype.kind == "f", f"{name} must be float, got dtype {arr.dtype}"


def align_trait_to_context(
    trait_df: pl.DataFrame,
    context_variants: pl.DataFrame,
    *,
    trait_total_sample_size: int | None,
    min_trait_snps: int,
) -> TraitAligned:
    """Align a trait sumstats dataframe onto the context SNP set.

    trait_df: gwaslab-standard columns rsID, EA, NEA, BETA, SE, and optionally N.
    context_variants: the context SNPs in row order, with rsID and the index EA/NEA.
    trait_total_sample_size: constant N to use when trait_df has no N column (else None).
    """
    has_n = GWASLAB_SAMPLE_SIZE_COLUMN in trait_df.columns
    if has_n:
        n_expr = pl.col(GWASLAB_SAMPLE_SIZE_COLUMN).cast(pl.Float64)
    else:
        assert trait_total_sample_size is not None, (
            "trait sumstats lack an N column; supply trait_total_sample_size in the config"
        )
        n_expr = pl.lit(float(trait_total_sample_size))
    trait = trait_df.select(
        pl.col(GWASLAB_RSID_COL),
        pl.col(GWASLAB_EFFECT_ALLELE_COL),
        pl.col(GWASLAB_NON_EFFECT_ALLELE_COL),
        (pl.col(GWASLAB_BETA_COL) / pl.col(GWASLAB_SE_COL)).alias(_TRAIT_Z),
        n_expr.alias(_TRAIT_N),
    ).unique(subset=[GWASLAB_RSID_COL], keep="first")

    joined = (
        context_variants.select(
            pl.col(GWASLAB_RSID_COL),
            pl.col(GWASLAB_EFFECT_ALLELE_COL).alias(_REF_EA),
            pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).alias(_REF_NEA),
        )
        .with_row_index(_CTX_ROW)
        .join(trait, on=GWASLAB_RSID_COL, how="left")
        .sort(_CTX_ROW)
    )

    match = match_reference_expr(
        ea_col=GWASLAB_EFFECT_ALLELE_COL,
        nea_col=GWASLAB_NON_EFFECT_ALLELE_COL,
        ref_ea_col=_REF_EA,
        ref_nea_col=_REF_NEA,
    )
    flipped = match_flipped_reference_expr(
        ea_col=GWASLAB_EFFECT_ALLELE_COL,
        nea_col=GWASLAB_NON_EFFECT_ALLELE_COL,
        ref_ea_col=_REF_EA,
        ref_nea_col=_REF_NEA,
    )
    aligned = joined.select(
        pl.when(match)
        .then(pl.col(_TRAIT_Z))
        .when(flipped)
        .then(-pl.col(_TRAIT_Z))
        .otherwise(None)
        .alias(_TRAIT_Z),
        pl.when(match | flipped).then(pl.col(_TRAIT_N)).otherwise(None).alias(_TRAIT_N),
    )

    z = aligned[_TRAIT_Z].to_numpy().astype(float)
    n = aligned[_TRAIT_N].to_numpy().astype(float)
    present = int(np.isfinite(z).sum())
    assert present > min_trait_snps, (
        f"trait matched only {present} of {context_variants.height} context SNPs "
        f"(<= min_trait_snps={min_trait_snps}); this usually means a harmonization failure "
        "in the trait input (wrong genome build, malformed rsIDs, or systematic allele mismatch)"
    )
    return TraitAligned(z=z, n=n)
