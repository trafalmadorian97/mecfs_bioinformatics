import narwhals
import numpy as np
import polars as pl
import scipy.stats

from mecfs_bio.build_system.task.pipes.collapse_multi_assay_proteins_pipe import (
    CollapseMultiAssayProteinsPipe,
)
from mecfs_bio.constants.ppp_ldsc_constants import (
    PPP_RG_GENE_COL,
    PPP_RG_H2_PROTEIN_COL,
    PPP_RG_N_ASSAYS_COL,
    PPP_RG_N_SNPS_COL,
    PPP_RG_OID_COL,
    PPP_RG_RG_COL,
    PPP_RG_RG_P_COL,
    PPP_RG_RG_SE_COL,
    PPP_RG_RG_SPREAD_COL,
)

# Columns the pipe sees (the display chain has already dropped rg_z / n_bar_* before it).
_INPUT = {
    PPP_RG_OID_COL: ["OID_B", "OID_A", "OID_C"],
    PPP_RG_GENE_COL: ["TNF", "TNF", "SINGLE"],
    "variant_set": ["all_variants", "all_variants", "all_variants"],
    PPP_RG_RG_COL: [0.30, 0.40, 0.20],
    PPP_RG_RG_SE_COL: [0.20, 0.10, 0.05],
    PPP_RG_RG_P_COL: [0.1, 0.05, 0.001],  # stale per-assay p; must be recomputed
    "gcov": [0.04, 0.05, 0.03],
    "gcov_intercept": [0.02, 0.01, 0.005],
    "h2_trait": [0.16, 0.16, 0.16],
    PPP_RG_H2_PROTEIN_COL: [0.10, 0.09, 0.05],
    PPP_RG_N_SNPS_COL: [1001, 1000, 1002],
}


def _run() -> pl.DataFrame:
    nw_df = narwhals.from_native(pl.DataFrame(_INPUT)).lazy()
    pipe = CollapseMultiAssayProteinsPipe(
        duplicate_oid_to_group=(("OID_A", "P_TNF"), ("OID_B", "P_TNF"))
    )
    return pipe.process(nw_df).collect().to_polars()


def test_duplicated_protein_is_conservatively_combined():
    out = _run()
    tnf = out.filter(pl.col(PPP_RG_OID_COL) == "OID_A").row(0, named=True)

    # Inverse-variance combine with perfect-correlation (linearly added) SE.
    inv_var = np.array([1 / 0.10**2, 1 / 0.20**2])
    rg = np.array([0.40, 0.30])
    se = np.array([0.10, 0.20])
    expect_rg = float((rg * inv_var).sum() / inv_var.sum())
    expect_se = float((1 / se).sum() / inv_var.sum())

    assert tnf[PPP_RG_OID_COL] == "OID_A"  # representative = smallest OID
    assert tnf[PPP_RG_GENE_COL] == "TNF"
    assert tnf[PPP_RG_N_ASSAYS_COL] == 2
    assert np.isclose(tnf[PPP_RG_RG_COL], expect_rg)
    assert np.isclose(tnf[PPP_RG_RG_SE_COL], expect_se)
    # SE is never smaller than the most precise single assay (no manufactured precision).
    assert tnf[PPP_RG_RG_SE_COL] >= se.min()
    assert np.isclose(tnf[PPP_RG_RG_SPREAD_COL], 0.10)
    expect_p = float(2.0 * scipy.stats.norm.sf(abs(expect_rg / expect_se)))
    assert np.isclose(tnf[PPP_RG_RG_P_COL], expect_p)
    # Diagnostics: averaged h2_protein, min n_snps.
    assert np.isclose(tnf[PPP_RG_H2_PROTEIN_COL], 0.095)
    assert tnf[PPP_RG_N_SNPS_COL] == 1000


def test_singleton_protein_passes_through_unchanged():
    out = _run()
    single = out.filter(pl.col(PPP_RG_OID_COL) == "OID_C").row(0, named=True)

    assert single[PPP_RG_N_ASSAYS_COL] == 1
    assert np.isclose(single[PPP_RG_RG_COL], 0.20)
    assert np.isclose(single[PPP_RG_RG_SE_COL], 0.05)
    assert np.isclose(single[PPP_RG_RG_SPREAD_COL], 0.0)
    expect_p = float(2.0 * scipy.stats.norm.sf(abs(0.20 / 0.05)))
    assert np.isclose(single[PPP_RG_RG_P_COL], expect_p)
