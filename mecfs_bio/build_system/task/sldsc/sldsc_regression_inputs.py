"""
Turn gwaslab sumstats plus prepared S-LDSC references into kernel inputs, reproducing the
preprocessing of ldsc's cell_type_specific as wrapped by gwaslab's estimate_h2_cts_by_ldsc:
rows restricted to gwaslab's bundled HapMap3 list (rsID present, and EA/NEA equal to the HapMap3
A1/A2 in either orientation; the list depends on the sumstats' genome build), Z from BETA/SE,
rows with missing SNP/Z/alleles dropped, an inner join onto the baseline SNP order (the prepared baseline
is already restricted to SNPs with a weight LD score, matching ldsc's inner join with w_ld),
N cast to int64 by truncation, and the chi-square filter Z^2 < max(0.001 * max(N), 80).

Fail-fast checks with no gwaslab counterpart: fewer than min_hapmap3_snps sumstats SNPs in the
baseline's HapMap3 regression-SNP universe raises; a regression SNP with missing N raises
(gwaslab fills with the median N); a SNP appearing more than once raises (gwaslab's merge would
duplicate the row), unless the duplicate-rsID policy is drop, in which case every row of such a
SNP is removed.

M_5_50 (passed through to the kernel with each LD-score column): for each annotation, the sum of
its values over reference-panel SNPs with minor allele frequency above 5% (for a binary
annotation, the number of common SNPs it contains). In the LDSC model
E[chi2_j] = N_j * sum_c tau_c * l(j, c) + N_j * a + 1, with l(j, c) the LD score of SNP j with
annotation c and tau_c the per-SNP heritability contributed by annotation c, M_5_50 converts
per-SNP effects into the heritability attributed to an annotation (h2_c = tau_c * M_5_50_c).
"""

import narwhals
import numpy as np
import polars as pl
import structlog
from attrs import frozen
from gwaslab.bd.bd_get_hapmap3 import _get_hapmap_full_polars

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.read_spec.read_sumstats import read_sumstats
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.sldsc.batched_sldsc_cts import (
    SLDSCModel,
    SLDSCRegressionInputs,
)
from mecfs_bio.build_system.task.sldsc.prepared_sldsc_reference import (
    SNP_COL,
    W_LD_COL,
    LdScoreTable,
    PreparedAnnotationSet,
    PreparedBaseline,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
    GWASLAB_Z_COL,
)
from mecfs_bio.constants.sldsc_constants import SLDSCDuplicateRsidPolicy

logger = structlog.get_logger()

_CHISQ_FILTER_N_COEF = 0.001
_CHISQ_FILTER_FLOOR = 80.0
_ROW_COL = "_row"
_ORDER_COL = "_order"
_HAPMAP3_RSID_COL = "rsid"
_HAPMAP3_A1_COL = "A1"
_HAPMAP3_A2_COL = "A2"


@frozen(slots=True)
class AlignedSLDSCRegression:
    inputs: SLDSCRegressionInputs
    models: tuple[SLDSCModel, ...]


def load_sldsc_sumstats(asset: Asset, pre_pipe: DataProcessingPipe) -> pl.DataFrame:
    """SNP, Z, N (Float64, possibly null) for every usable HapMap3 row of the sumstats."""
    sumstats = read_sumstats(asset)
    data = (
        pre_pipe.process(narwhals.from_native(sumstats.data).lazy())
        .collect()
        .to_polars()
    )
    required = (
        GWASLAB_RSID_COL,
        GWASLAB_SAMPLE_SIZE_COLUMN,
        GWASLAB_EFFECT_ALLELE_COL,
        GWASLAB_NON_EFFECT_ALLELE_COL,
    )
    missing = [c for c in required if c not in data.columns]
    assert not missing, f"sumstats need columns {missing} for S-LDSC"
    z = (
        pl.col(GWASLAB_Z_COL)
        if GWASLAB_Z_COL in data.columns
        else pl.col(GWASLAB_BETA_COL) / pl.col(GWASLAB_SE_COL)
    )
    selected = data.select(
        pl.col(GWASLAB_RSID_COL).cast(pl.String).alias(SNP_COL),
        z.cast(pl.Float64).alias(GWASLAB_Z_COL),
        pl.col(GWASLAB_SAMPLE_SIZE_COLUMN).cast(pl.Float64).fill_nan(None),
        pl.col(GWASLAB_EFFECT_ALLELE_COL).cast(pl.String),
        pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).cast(pl.String),
    )
    return (
        _restrict_to_hapmap3(sumstats=selected, build=sumstats.build)
        .drop_nulls(
            subset=[
                SNP_COL,
                GWASLAB_Z_COL,
                GWASLAB_EFFECT_ALLELE_COL,
                GWASLAB_NON_EFFECT_ALLELE_COL,
            ]
        )
        .filter(pl.col(GWASLAB_Z_COL).is_not_nan())
        .select(SNP_COL, GWASLAB_Z_COL, GWASLAB_SAMPLE_SIZE_COLUMN)
    )


def _restrict_to_hapmap3(sumstats: pl.DataFrame, build: str) -> pl.DataFrame:
    """
    gwaslab's _get_hapmap3 as called by estimate_h2_cts_by_ldsc: keep rows whose rsID is in the
    HapMap3 list for this build and whose EA/NEA equal the HapMap3 A1/A2 in either orientation.
    Duplicate sumstats rsIDs are kept (and rejected later), as gwaslab's merge keeps them.
    """
    hapmap3 = _get_hapmap_full_polars(build, include_alleles=True).select(
        pl.col(_HAPMAP3_RSID_COL).alias(SNP_COL),
        pl.col(_HAPMAP3_A1_COL),
        pl.col(_HAPMAP3_A2_COL),
    )
    ea = pl.col(GWASLAB_EFFECT_ALLELE_COL)
    nea = pl.col(GWASLAB_NON_EFFECT_ALLELE_COL)
    a1 = pl.col(_HAPMAP3_A1_COL)
    a2 = pl.col(_HAPMAP3_A2_COL)
    return (
        sumstats.join(hapmap3, on=SNP_COL, how="inner", maintain_order="left")
        .filter(((ea == a1) & (nea == a2)) | ((ea == a2) & (nea == a1)))
        .drop(_HAPMAP3_A1_COL, _HAPMAP3_A2_COL)
    )


def build_regression_inputs(
    sumstats: pl.DataFrame,
    baseline: PreparedBaseline,
    annotation_set: PreparedAnnotationSet,
    min_hapmap3_snps: int,
    duplicate_rsids: SLDSCDuplicateRsidPolicy,
) -> AlignedSLDSCRegression:
    if duplicate_rsids == "drop":
        sumstats = _drop_duplicate_rsids(sumstats)
    joined = _join_baseline(sumstats=sumstats, baseline=baseline)
    _require_min_snps(joined=joined, min_hapmap3_snps=min_hapmap3_snps)
    regression = _filter_chisq(joined=_cast_n(joined=joined))
    baseline_rows = regression[_ROW_COL].to_numpy()
    annotation_rows = _annotation_rows(regression=regression, ld=annotation_set.ld)
    z = regression[GWASLAB_Z_COL].to_numpy()
    return AlignedSLDSCRegression(
        inputs=SLDSCRegressionInputs(
            chi2=z * z,
            n=regression[GWASLAB_SAMPLE_SIZE_COLUMN].to_numpy().astype(np.float64),
            w_ld=baseline.w_ld[W_LD_COL].to_numpy()[baseline_rows].astype(np.float64),
            baseline_ld=baseline.ld.matrix[baseline_rows].astype(np.float64),
            baseline_m=baseline.ld.m_5_50,
            annotation_ld=annotation_set.ld.matrix[annotation_rows],
            annotation_m=annotation_set.ld.m_5_50,
        ),
        models=tuple(
            SLDSCModel(
                label=m.label,
                target_index=annotation_set.ld.column_index(m.target_column),
                control_indices=tuple(
                    annotation_set.ld.column_index(c) for c in m.control_columns
                ),
            )
            for m in annotation_set.models
        ),
    )


def _drop_duplicate_rsids(sumstats: pl.DataFrame) -> pl.DataFrame:
    """Remove every row whose rsID occurs more than once."""
    duplicated = pl.col(SNP_COL).is_duplicated()
    count = sumstats.select(duplicated.sum()).item()
    if count:
        logger.warning(
            "dropping rows with duplicate S-LDSC regression rsIDs", count=count
        )
    return sumstats.filter(~duplicated)


def _join_baseline(sumstats: pl.DataFrame, baseline: PreparedBaseline) -> pl.DataFrame:
    """Sumstats rows in the baseline SNP universe, in baseline row order, with that row index."""
    joined = (
        baseline.ld.snps.select(SNP_COL)
        .with_row_index(_ROW_COL)
        .join(sumstats, on=SNP_COL, how="inner")
        .sort(_ROW_COL)
    )
    assert joined[SNP_COL].is_unique().all(), (
        "sumstats contain duplicate rsIDs among S-LDSC regression SNPs"
    )
    return joined


def _require_min_snps(joined: pl.DataFrame, min_hapmap3_snps: int) -> None:
    assert joined.height >= min_hapmap3_snps, (
        f"only {joined.height} sumstats SNPs are in the S-LDSC HapMap3 regression-SNP set; "
        f"at least {min_hapmap3_snps} are required"
    )


def _cast_n(joined: pl.DataFrame) -> pl.DataFrame:
    n_missing = joined[GWASLAB_SAMPLE_SIZE_COLUMN].null_count()
    assert n_missing == 0, (
        f"{n_missing} S-LDSC regression SNPs have no sample size; this indicates an upstream "
        "data-processing problem"
    )
    # gwaslab casts N with .apply("int64"): truncation toward zero.
    return joined.with_columns(pl.col(GWASLAB_SAMPLE_SIZE_COLUMN).cast(pl.Int64))


def _filter_chisq(joined: pl.DataFrame) -> pl.DataFrame:
    n_max = joined[GWASLAB_SAMPLE_SIZE_COLUMN].max()
    assert isinstance(n_max, int)
    chisq_max = max(_CHISQ_FILTER_N_COEF * n_max, _CHISQ_FILTER_FLOOR)
    return joined.filter(pl.col(GWASLAB_Z_COL) ** 2 < chisq_max)


def _annotation_rows(regression: pl.DataFrame, ld: LdScoreTable) -> np.ndarray:
    """Row of each regression SNP in the annotation set's LD matrix, in regression order."""
    rows = (
        regression.select(SNP_COL)
        .with_row_index(_ORDER_COL)
        .join(ld.snps.select(SNP_COL).with_row_index(_ROW_COL), on=SNP_COL, how="left")
        .sort(_ORDER_COL)[_ROW_COL]
    )
    assert rows.null_count() == 0, (
        f"{rows.null_count()} regression SNPs are missing from the S-LDSC annotation set"
    )
    return rows.to_numpy()
