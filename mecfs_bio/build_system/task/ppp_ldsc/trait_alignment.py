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

Alignment also reports counts for callers' own upstream-bug checks: how many context rsIDs the
trait carries (n_found), how many of those have alleles matching the index directly or swapped
(n_matched), and, when both the context variants and the trait carry CHR and POS, how many matched
variants have a usable trait position (n_position_compared) and how many of those equal the
context's (n_position_concordant).

The trait arrives lazily and is narrowed to the context variants before anything is materialized.
Only a small fraction of a trait's variants are HapMap3 ones, and a GWAS on whole-genome sequencing
data can carry hundreds of millions of rows, so restricting to the context first is what keeps that
input's size off the peak memory of the run.

The trait is a narwhals LazyFrame because it has been through a DataProcessingPipe, which is free
to change the backend. The filter is therefore expressed in narwhals and only the (small) filtered
result is collected into polars, where the allele-orientation helpers live. Converting to polars
any earlier would mean collecting the whole trait first, which is exactly what the filter exists to
avoid.
"""

from __future__ import annotations

import narwhals
import numpy as np
import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_rsid import (
    match_flipped_reference_expr,
    match_reference_expr,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

logger = structlog.get_logger()
_CTX_ROW = "__ctx_row__"
_REF_EA = "__ref_ea__"
_REF_NEA = "__ref_nea__"
_TRAIT_Z = "__trait_z__"
_TRAIT_N = "__trait_n__"
_TRAIT_CHR = "__trait_chr__"
_TRAIT_POS = "__trait_pos__"
_REF_CHR = "__ref_chr__"
_REF_POS = "__ref_pos__"
_FOUND = "__found__"


@frozen(slots=True)
class TraitAligned:
    """The trait aligned to the context SNP set (parallel (S,) arrays in context row order).

    z: signed trait z-score oriented to the index effect allele, NaN where the trait lacks the
    variant or its alleles do not match. n: trait sample size per SNP, NaN where z is NaN.

    n_found: context variants whose rsID the trait carries. n_matched: those whose alleles
    match the index directly or swapped. n_position_compared / n_position_concordant: matched
    variants with a usable trait (CHR, POS), and those whose (CHR, POS) equals the context's;
    both None when positions were not available on both sides.
    """

    z: np.ndarray
    n: np.ndarray
    n_found: int
    n_matched: int
    n_position_compared: int | None
    n_position_concordant: int | None

    def __attrs_post_init__(self) -> None:
        s = self.z.shape[0]
        for name, arr in (("z", self.z), ("n", self.n)):
            assert arr.ndim == 1, f"{name} must be 1-D, got shape {arr.shape}"
            assert arr.shape[0] == s, f"{name} has length {arr.shape[0]}, expected {s}"
            assert arr.dtype.kind == "f", f"{name} must be float, got dtype {arr.dtype}"
        assert 0 <= self.n_matched <= self.n_found <= s, (
            f"need 0 <= n_matched ({self.n_matched}) <= n_found ({self.n_found}) <= {s}"
        )
        assert (self.n_position_compared is None) == (
            self.n_position_concordant is None
        ), "position counts must be both set or both None"
        if (
            self.n_position_compared is not None
            and self.n_position_concordant is not None
        ):
            assert (
                0
                <= self.n_position_concordant
                <= self.n_position_compared
                <= self.n_matched
            ), "position counts must satisfy concordant <= compared <= n_matched"


def align_trait_to_context(
    trait_df: narwhals.LazyFrame,
    context_variants: pl.DataFrame,
    *,
    trait_total_sample_size: int | None,
    min_trait_snps: int,
) -> TraitAligned:
    """Align a trait sumstats dataframe onto the context SNP set.

    trait_df: a LAZY narwhals frame with gwaslab-standard columns rsID, EA, NEA, BETA, SE, and
        optionally N. Kept lazy so the filter down to the context variants happens during the scan,
        and backend-agnostic because a DataProcessingPipe may have produced it.
    context_variants: the context SNPs in row order, with rsID and the index EA/NEA, and
        optionally CHR and POS (needed for the position counts).
    trait_total_sample_size: constant N to use when trait_df has no N column (else None).
    """
    trait_cols = trait_df.collect_schema().names()
    has_n = GWASLAB_SAMPLE_SIZE_COLUMN in trait_cols
    position_cols = (GWASLAB_CHROM_COL, GWASLAB_POS_COL)
    compare_positions = all(c in trait_cols for c in position_cols) and all(
        c in context_variants.columns for c in position_cols
    )
    if has_n:
        n_expr = pl.col(GWASLAB_SAMPLE_SIZE_COLUMN).cast(pl.Float64)
    else:
        assert trait_total_sample_size is not None, (
            "trait sumstats lack an N column; supply trait_total_sample_size in the config"
        )
        n_expr = pl.lit(float(trait_total_sample_size))
    projected = [
        GWASLAB_RSID_COL,
        GWASLAB_EFFECT_ALLELE_COL,
        GWASLAB_NON_EFFECT_ALLELE_COL,
        _TRAIT_Z,
        _TRAIT_N,
    ]
    position_exprs: list[pl.Expr] = []
    if compare_positions:
        # strict=False: non-numeric chromosomes (X, Y) become null and are simply not compared.
        position_exprs = [
            pl.col(GWASLAB_CHROM_COL).cast(pl.Int64, strict=False).alias(_TRAIT_CHR),
            pl.col(GWASLAB_POS_COL).cast(pl.Int64, strict=False).alias(_TRAIT_POS),
        ]
        projected += [_TRAIT_CHR, _TRAIT_POS]
    trait = (
        trait_df
        # Applied while lazy, so a polars-backed scan evaluates it as a parquet SELECTION and only
        # context variants are ever materialized. The rsIDs go in as a plain list: it is the one
        # form every narwhals backend accepts, and it avoids polars' ambiguous-elementwise
        # deprecation for a bare Series.
        .filter(
            narwhals.col(GWASLAB_RSID_COL).is_in(
                context_variants[GWASLAB_RSID_COL].to_list()
            )
        )
        # Small by now -- at most one row per context variant, plus any duplicates -- so
        # collecting into polars here costs little and gives the rest of the function the polars
        # expressions the orientation helpers are written in.
        .collect(backend="polars")
        .to_polars()
        .select(
            pl.col(GWASLAB_RSID_COL),
            pl.col(GWASLAB_EFFECT_ALLELE_COL),
            pl.col(GWASLAB_NON_EFFECT_ALLELE_COL),
            (pl.col(GWASLAB_BETA_COL) / pl.col(GWASLAB_SE_COL)).alias(_TRAIT_Z),
            n_expr.alias(_TRAIT_N),
            *position_exprs,
            pl.lit(True).alias(_FOUND),
        )
        # Sorting on every projected column before dropping duplicate rsIDs makes the survivor a
        # function of the data alone. Taking whichever row the scan happened to emit first would
        # leave the result -- and so the asset's hash -- dependent on how the file was read.
        .sort(projected)
        .unique(subset=[GWASLAB_RSID_COL], keep="first", maintain_order=True)
    )

    ref_position_exprs: list[pl.Expr] = []
    if compare_positions:
        ref_position_exprs = [
            pl.col(GWASLAB_CHROM_COL).cast(pl.Int64).alias(_REF_CHR),
            pl.col(GWASLAB_POS_COL).cast(pl.Int64).alias(_REF_POS),
        ]
    joined = (
        context_variants.select(
            pl.col(GWASLAB_RSID_COL),
            pl.col(GWASLAB_EFFECT_ALLELE_COL).alias(_REF_EA),
            pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).alias(_REF_NEA),
            *ref_position_exprs,
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

    matched = joined.select((match | flipped).fill_null(False)).to_series()
    n_found = int(joined[_FOUND].fill_null(False).sum())
    n_matched = int(matched.sum())
    n_position_compared: int | None = None
    n_position_concordant: int | None = None
    if compare_positions:
        comparable = joined.filter(
            matched
            & pl.col(_TRAIT_CHR).is_not_null()
            & pl.col(_TRAIT_POS).is_not_null()
        )
        n_position_compared = comparable.height
        n_position_concordant = comparable.filter(
            (pl.col(_TRAIT_CHR) == pl.col(_REF_CHR))
            & (pl.col(_TRAIT_POS) == pl.col(_REF_POS))
        ).height

    z = aligned[_TRAIT_Z].to_numpy().astype(float)
    n = aligned[_TRAIT_N].to_numpy().astype(float)
    present = int(np.isfinite(z).sum())
    assert present > min_trait_snps, (
        f"trait matched only {present} of {context_variants.height} context SNPs "
        f"(<= min_trait_snps={min_trait_snps}); this usually means a harmonization failure "
        "in the trait input (wrong genome build, malformed rsIDs, or systematic allele mismatch)"
    )
    logger.debug(f"Trait matched {present} of {context_variants.height} context SNPs")
    return TraitAligned(
        z=z,
        n=n,
        n_found=n_found,
        n_matched=n_matched,
        n_position_compared=n_position_compared,
        n_position_concordant=n_position_concordant,
    )
