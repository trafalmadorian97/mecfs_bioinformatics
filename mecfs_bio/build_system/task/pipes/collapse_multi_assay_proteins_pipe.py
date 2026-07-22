"""
Collapse the per-assay UKB-PPP rg table to one row per protein.

A handful of proteins are measured on more than one Olink Explore panel (e.g. TNF, IL6, CXCL8,
IDO1, LMOD1, SCRIB on all four v1 panels), so they appear as several assay rows -- same UniProt,
different OID -- each with its own genetic-correlation estimate. Keeping all of them distorts a
ranked / thresholded discussion list: a multi-panel protein gets several chances to cross a
significance cut-off, and if it enters via its most significant panel the reported rg is
winner's-curse inflated.

This pipe combines a protein's assay rows into a single row. Because the assay estimates share the
same trait GWAS (and largely the same protein signal), they are strongly positively correlated, so
a naive fixed-effect meta-analysis would badly understate the combined SE. We therefore combine
under the conservative worst-case assumption of perfect correlation between assays, for which
standard deviations add linearly:

    rg_comb = sum(rg_i / se_i^2) / sum(1 / se_i^2)      (inverse-variance weighted mean)
    se_comb = sum(1 / se_i)     / sum(1 / se_i^2)       (same-weighted mean of the SEs)

At perfect correlation Var(sum w_i rg_i) = (sum w_i se_i)^2, and Var is monotincreasing in the
pairwise correlation, so se_comb is never smaller than the truth for any real correlation. The
point estimate is a weighted mean (not the max), so one lucky panel cannot inflate it. z and p are
recomputed from the combined estimate.

Only duplicated OIDs need combining; every other OID forms a singleton group and passes through
unchanged (the formulas above are a no-op on a single row). The pipe carries the small
OID -> shared-key map for the duplicated OIDs only.
"""

from __future__ import annotations

import narwhals
import numpy as np
import polars as pl
import scipy.stats
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.ppp_ldsc_constants import (
    PPP_RG_GCOV_COL,
    PPP_RG_GCOV_INTERCEPT_COL,
    PPP_RG_GENE_COL,
    PPP_RG_H2_PROTEIN_COL,
    PPP_RG_H2_TRAIT_COL,
    PPP_RG_N_ASSAYS_COL,
    PPP_RG_N_SNPS_COL,
    PPP_RG_OID_COL,
    PPP_RG_RG_COL,
    PPP_RG_RG_P_COL,
    PPP_RG_RG_SE_COL,
    PPP_RG_RG_SPREAD_COL,
    PPP_RG_VARIANT_SET_COL,
)

# Temporary group key and partial-sum columns (dropped before returning).
_GROUP = "__protein_group__"
_SUM_RG_OVER_VAR = "__sum_rg_over_var__"
_SUM_INV_VAR = "__sum_inv_var__"
_SUM_INV_SE = "__sum_inv_se__"


@frozen
class CollapseMultiAssayProteinsPipe(DataProcessingPipe):
    """Combine multi-panel assay rows of the same protein into one conservatively-meta-analysed row.

    duplicate_oid_to_group: for each OID that shares its UniProt with another assay, the shared
        group key (its UniProt). OIDs absent from this map form singleton groups. Passing only the
        duplicated OIDs keeps the map (and the task identity) small; supply it from the manifest,
        e.g. the UniProt of every OID whose UniProt occurs more than once.
    """

    duplicate_oid_to_group: tuple[tuple[str, str], ...]

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        df = x.collect().to_polars()

        mapping = dict(self.duplicate_oid_to_group)
        group_expr = (
            pl.col(PPP_RG_OID_COL).replace(mapping)
            if mapping
            else pl.col(PPP_RG_OID_COL)
        )
        collapsed = (
            df.with_columns(group_expr.alias(_GROUP))
            .group_by(_GROUP, maintain_order=True)
            .agg(
                # Representative identity: the row with the smallest OID in the group.
                pl.col(PPP_RG_OID_COL).min().alias(PPP_RG_OID_COL),
                pl.col(PPP_RG_GENE_COL)
                .sort_by(PPP_RG_OID_COL)
                .first()
                .alias(PPP_RG_GENE_COL),
                pl.col(PPP_RG_VARIANT_SET_COL).first().alias(PPP_RG_VARIANT_SET_COL),
                pl.len().alias(PPP_RG_N_ASSAYS_COL),
                (pl.col(PPP_RG_RG_COL).max() - pl.col(PPP_RG_RG_COL).min()).alias(
                    PPP_RG_RG_SPREAD_COL
                ),
                # Conservative (perfect-correlation) inverse-variance combination.
                (pl.col(PPP_RG_RG_COL) / pl.col(PPP_RG_RG_SE_COL) ** 2)
                .sum()
                .alias(_SUM_RG_OVER_VAR),
                (1.0 / pl.col(PPP_RG_RG_SE_COL) ** 2).sum().alias(_SUM_INV_VAR),
                (1.0 / pl.col(PPP_RG_RG_SE_COL)).sum().alias(_SUM_INV_SE),
                # Secondary diagnostics: averaged across assays (h2_trait is constant per run).
                pl.col(PPP_RG_GCOV_COL).mean().alias(PPP_RG_GCOV_COL),
                pl.col(PPP_RG_GCOV_INTERCEPT_COL)
                .mean()
                .alias(PPP_RG_GCOV_INTERCEPT_COL),
                pl.col(PPP_RG_H2_TRAIT_COL).first().alias(PPP_RG_H2_TRAIT_COL),
                pl.col(PPP_RG_H2_PROTEIN_COL).mean().alias(PPP_RG_H2_PROTEIN_COL),
                pl.col(PPP_RG_N_SNPS_COL).min().alias(PPP_RG_N_SNPS_COL),
            )
            .with_columns(
                (pl.col(_SUM_RG_OVER_VAR) / pl.col(_SUM_INV_VAR)).alias(PPP_RG_RG_COL),
                (pl.col(_SUM_INV_SE) / pl.col(_SUM_INV_VAR)).alias(PPP_RG_RG_SE_COL),
            )
            .drop(_GROUP, _SUM_RG_OVER_VAR, _SUM_INV_VAR, _SUM_INV_SE)
        )

        z = (collapsed[PPP_RG_RG_COL] / collapsed[PPP_RG_RG_SE_COL]).to_numpy()
        collapsed = collapsed.with_columns(
            pl.Series(PPP_RG_RG_P_COL, 2.0 * scipy.stats.norm.sf(np.abs(z)))
        )

        ordered = collapsed.select(
            PPP_RG_OID_COL,
            PPP_RG_GENE_COL,
            PPP_RG_VARIANT_SET_COL,
            PPP_RG_RG_COL,
            PPP_RG_RG_SE_COL,
            PPP_RG_RG_P_COL,
            PPP_RG_RG_SPREAD_COL,
            PPP_RG_N_ASSAYS_COL,
            PPP_RG_GCOV_COL,
            PPP_RG_GCOV_INTERCEPT_COL,
            PPP_RG_H2_TRAIT_COL,
            PPP_RG_H2_PROTEIN_COL,
            PPP_RG_N_SNPS_COL,
        ).sort(PPP_RG_RG_P_COL)
        return narwhals.from_native(ordered).lazy()
