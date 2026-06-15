"""
Shared configuration, source descriptions, and constants for the GenomicSEM
task family.

This module is deliberately free of rpy2/R: it holds only the pure-data
dataclasses, type aliases, and string constants that both the R-backed tasks
and the full-Python tasks need. Keeping it R-free lets the Python tasks (and
their tests) depend on it without dragging in rpy2.
"""

from __future__ import annotations

from attrs import frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe

# --- Column names expected by GenomicSEM::munge by default. ---
MUNGE_SNP_COL = "SNP"
MUNGE_A1_COL = "A1"
MUNGE_A2_COL = "A2"
MUNGE_EFFECT_COL = "effect"
MUNGE_SE_COL = "SE"
MUNGE_P_COL = "P"
MUNGE_N_COL = "N"
MUNGE_MAF_COL = "MAF"
MUNGE_INFO_COL = "INFO"  # optional input column (imputation quality)
MUNGE_Z_COL = "Z"  # Z-score column in the munged output

# --- Columns in the LD-score reference files (`<chr>.l2.ldscore.gz`) read by
#     GenomicSEM::ldsc. SNP is shared with the munged columns (MUNGE_SNP_COL). ---
LDSC_CHR_COL = "CHR"
LDSC_BP_COL = "BP"
LDSC_L2_COL = "L2"  # per-SNP LD score

# --- Output structure produced by the basic GenomicSEMTask. ---
MUNGED_SUBDIR = "munged"
LDSC_LOG_PREFIX = "genomic_sem"  # GenomicSEM::ldsc appends "_ldsc.log" to this
LDSC_LOG_FILENAME = f"{LDSC_LOG_PREFIX}_ldsc.log"
LDSC_S_FILENAME = "ldsc_S.csv"
LDSC_V_FILENAME = "ldsc_V.csv"
LDSC_S_STAND_FILENAME = "ldsc_S_Stand.csv"
LAVAAN_MODEL_FILENAME = "lavaan_model.txt"
MODEL_FIT_FILENAME = "model_fit.csv"
MODEL_RESULTS_FILENAME = "model_results.csv"

# --- Output structure produced by the GWAS-extension tasks. ---
GWAS_RESULTS_SUBDIR = "gwas_results"
COMMON_FACTOR_GWAS_FILENAME = "common_factor.parquet"
SUBTRACTION_F_FILENAME = "common_factor.parquet"  # F: factor shared by both traits
SUBTRACTION_R_FILENAME = "remainder_factor.parquet"  # R: residual unique to T1

# --- Columns of the per-factor subtraction result tables (make_result_df).
#     SNP/CHR/BP/MAF/A1/A2 reuse the munge/LDSC column constants; the rest are
#     the subtraction-specific estimate columns. ---
SUBTRACTION_LHS_COL = "lhs"  # left-hand side of the model term (factor name)
SUBTRACTION_OP_COL = "op"  # model operator ("~")
SUBTRACTION_RHS_COL = "rhs"  # right-hand side of the model term ("SNP")
SUBTRACTION_EST_COL = "est"  # per-SNP factor effect estimate
SUBTRACTION_SE_COL = "se_c"  # sandwich standard error
SUBTRACTION_Z_COL = "Z_Estimate"
SUBTRACTION_P_COL = "Pval_Estimate"
SUBTRACTION_N_EFF_COL = "N_eff"
SUBTRACTION_FAIL_COL = "fail"


@frozen
class GenomicSEMSumstatsSource:
    """
    A single trait sumstats source for GWAS-by-subtraction.

    The source Task is expected to produce tabular GWAS summary statistics
    with columns named according to gwaslab conventions
    (see mecfs_bio.constants.gwaslab_constants).

    Every trait is standardised on a single linear (OLS-style) scale and no
    liability conversion is applied, so the source carries only a sample size
    (no case/control prevalences). See ``_genomic_sem_sumstats`` and
    ``genomic_sem_gwas_by_subtraction_full_python_task`` for why the binary /
    logistic special-casing was removed.

    ``sample_size`` is the total GWAS sample size, used to standardise the
    per-SNP effects. ``None`` means "take N from the source's own N column".
    """

    task: Task
    alias: str
    sample_size: float | None = None
    pipe: DataProcessingPipe = IdentityPipe()

    @property
    def asset_id(self) -> AssetId:
        return self.task.asset_id
