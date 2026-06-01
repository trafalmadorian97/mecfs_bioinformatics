"""
Shared configuration, source descriptions, and constants for the GenomicSEM
task family.

This module is deliberately free of rpy2/R: it holds only the pure-data
dataclasses, type aliases, and string constants that both the R-backed tasks
and the full-Python tasks need. Keeping it R-free lets the Python tasks (and
their tests) depend on it without dragging in rpy2.
"""

from __future__ import annotations

from typing import Final, Literal

from attrs import frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    PhenotypeInfo,
)
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

GSEstmationType = Literal["DWLS", "ML"]

# How the source GWAS was estimated. Controls the per-trait flag (one of
# se.logit / OLS / linprob) passed to GenomicSEM::sumstats.
OLS: Final = "ols"
LOGISTIC: Final = "logistic"
LINEAR_PROB: Final = "linear_prob"

# `ty` requires the type arguments to Literal to be inline literals rather
# than Final-typed names, so the strings are repeated here.
GWASMethod = Literal["ols", "logistic", "linear_prob"]


@frozen
class GenomicSEMSumstatsSource:
    """
    A single trait sumstats source for GenomicSEM.

    The source Task is expected to produce tabular GWAS summary statistics
    with columns named according to gwaslab conventions
    (see mecfs_bio.constants.gwaslab_constants).
    """

    task: Task
    alias: str
    sample_info: PhenotypeInfo
    pipe: DataProcessingPipe = IdentityPipe()

    @property
    def asset_id(self) -> AssetId:
        return self.task.asset_id


@frozen
class GenomicSEMConfig:
    """
    Configuration for the basic GenomicSEM workflow (munge + ldsc + usermodel).
    """

    estimation: GSEstmationType = (
        "DWLS"  # GenomicSEM::usermodel estimation: "DWLS" or "ML"
    )
    std_lv: bool = False  # If TRUE, latent variables are standardised
    cfi_calc: bool = True  # Whether to compute CFI fit index
    fix_resid: bool = True  # Whether to fix residual variances
    info_filter: float = 0.9  # munge: minimum INFO score
    maf_filter: float = 0.01  # munge: minimum MAF
    # If non-empty, treated as a basename prefix on the LD reference files
    # (e.g. "LDscore." matches files named "LDscore.<chr>.l2.ldscore.gz").
    # GenomicSEM::ldsc hardcodes "<ld>/<chr>.l2.ldscore.gz", so when a prefix
    # is set we build a symlinked view of the dir with the prefix stripped and
    # pass that to ldsc.
    ld_file_filename_pattern: str = ""


@frozen
class GenomicSEMGWASSumstatsSource:
    """
    A source trait for GenomicSEM GWAS-extension workflows.

    Wraps the existing GenomicSEMSumstatsSource and adds the per-trait
    estimation method that sumstats() needs.
    """

    source: GenomicSEMSumstatsSource
    gwas_method: GWASMethod

    @property
    def task(self) -> Task:
        return self.source.task

    @property
    def alias(self) -> str:
        return self.source.alias

    @property
    def asset_id(self) -> AssetId:
        return self.source.asset_id


@frozen
class GenomicSEMSumstatsConfig:
    """Configuration for GenomicSEM::sumstats."""

    info_filter: float = 0.6  # GenomicSEM default — lower than munge's 0.9
    maf_filter: float = 0.01
    keep_indel: bool = False
    parallel: bool = False
    cores: int | None = None
    ambig: bool = False


@frozen
class GenomicSEMGWASRunConfig:
    """Configuration shared by commonfactorGWAS and userGWAS."""

    estimation: str = "DWLS"
    parallel: bool = True
    cores: int | None = None
    gc_correction: str = "standard"  # "standard", "conserv", or "none"
    toler: float | bool = False
    snp_se: float | bool = False
    smooth_check: bool = False
