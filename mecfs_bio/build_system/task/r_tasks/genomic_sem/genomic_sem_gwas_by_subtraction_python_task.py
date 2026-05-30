"""
GWAS-by-subtraction via munge+ldsc+sumstats in R, with the per-SNP
Cholesky decomposition and delta-method SE computed in numpy.

Convention: source[0] = composite trait (T1), source[1] = reference trait (T2).
Factor F is the genetic factor common to both traits (defined by the reference
trait T2); factor R is the remainder unique to the composite trait T1,
orthogonal to F. See `_gwas_by_subtraction_kernel` for the model.

The output is two parquets, one per factor (F~SNP and R~SNP), matching the
layout that GenomicSEM's userGWAS with sub=c("F~SNP","R~SNP") would produce.
"""

import gc
import tempfile
from pathlib import Path, PurePath
from typing import Sequence

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from attrs import frozen
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    MULTI_TRAIT,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._gwas_by_subtraction_kernel import (
    fit_gwas_by_subtraction,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_task import (
    GenomicSEMConfig,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_user_gwas_task import (
    GWAS_RESULTS_SUBDIR,
    GenomicSEMGWASRunConfig,
    GenomicSEMGWASSumstatsSource,
    GenomicSEMSumstatsConfig,
    _prepare_gwas_inputs,
    _resolve_file_path,
    _resolve_ld_path,
    _validate_sources,
    logger,
)
from mecfs_bio.build_system.wf.base_wf import WF

SUBTRACTION_F_FILENAME = "common_factor.parquet"  # F: factor shared by both traits
SUBTRACTION_R_FILENAME = "remainder_factor.parquet"  # R: residual unique to T1


@frozen
class GenomicSEMGWASBySubtractionPythonTask(Task):
    """
    GWAS-by-subtraction where the per-SNP Cholesky decomposition and SE are
    computed in numpy instead of lavaan.

    source[0] = composite trait (T1, e.g. educational attainment)
    source[1] = reference trait (T2, e.g. cognitive performance)

    Output: two parquets under gwas_results/ — one for the common factor
    (F~SNP) and one for the remainder factor (R~SNP).
    """

    meta: Meta
    sources: Sequence[GenomicSEMGWASSumstatsSource]
    ld_ref_task: Task
    hapmap_snps_task: Task
    sumstats_ref_task: Task
    munge_config: GenomicSEMConfig = GenomicSEMConfig()
    sumstats_config: GenomicSEMSumstatsConfig = GenomicSEMSumstatsConfig()
    run_config: GenomicSEMGWASRunConfig = GenomicSEMGWASRunConfig()

    def __attrs_post_init__(self):
        _validate_sources(self.sources)
        assert len(self.sources) == 2, (
            f"GWAS-by-subtraction requires exactly 2 traits; got {len(self.sources)}"
        )

    @property
    def deps(self) -> list[Task]:
        result: list[Task] = [
            self.ld_ref_task,
            self.hapmap_snps_task,
            self.sumstats_ref_task,
        ]
        for source in self.sources:
            result.append(source.task)
        return result

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        gsem = importr("GenomicSEM")
        ld_path = _resolve_ld_path(self.ld_ref_task, fetch, self.munge_config)
        hapmap_path = _resolve_file_path(self.hapmap_snps_task, fetch)
        sumstats_ref_path = _resolve_file_path(self.sumstats_ref_task, fetch)

        with tempfile.TemporaryDirectory() as tmp_dir_str:
            tmp_dir = Path(tmp_dir_str)
            covstruc_r, snps_r = _prepare_gwas_inputs(
                gsem=gsem,
                sources=self.sources,
                ld_path=ld_path,
                hapmap_path=hapmap_path,
                sumstats_ref_path=sumstats_ref_path,
                munge_config=self.munge_config,
                sumstats_config=self.sumstats_config,
                fetch=fetch,
                scratch_dir=scratch_dir,
                tmp_dir=tmp_dir,
            )
            logger.debug("Running Python GWAS-by-subtraction kernel")
            frames = _run_python_subtraction(
                covstruc_r=covstruc_r,
                snps_r=snps_r,
            )
            out_dir = scratch_dir / GWAS_RESULTS_SUBDIR
            out_dir.mkdir(parents=True, exist_ok=True)
            frames.f_df.to_parquet(out_dir / SUBTRACTION_F_FILENAME, index=False)
            frames.r_df.to_parquet(out_dir / SUBTRACTION_R_FILENAME, index=False)
            logger.debug(
                f"Wrote subtraction results: "
                f"F ({len(frames.f_df)} rows, {int(frames.f_df['fail'].sum())} fails), "
                f"R ({len(frames.r_df)} rows, {int(frames.r_df['fail'].sum())} fails)"
            )

        gc.collect()
        ro.r("gc()")
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        sources: Sequence[GenomicSEMGWASSumstatsSource],
        ld_ref_task: Task,
        hapmap_snps_task: Task,
        sumstats_ref_task: Task,
        munge_config: GenomicSEMConfig = GenomicSEMConfig(),
        sumstats_config: GenomicSEMSumstatsConfig = GenomicSEMSumstatsConfig(),
        run_config: GenomicSEMGWASRunConfig = GenomicSEMGWASRunConfig(),
    ) -> "GenomicSEMGWASBySubtractionPythonTask":
        meta = ResultDirectoryMeta(
            id=AssetId(asset_id),
            trait=MULTI_TRAIT,
            project="genomic_sem",
            sub_dir=PurePath("analysis"),
        )
        return cls(
            meta=meta,
            sources=sources,
            ld_ref_task=ld_ref_task,
            hapmap_snps_task=hapmap_snps_task,
            sumstats_ref_task=sumstats_ref_task,
            munge_config=munge_config,
            sumstats_config=sumstats_config,
            run_config=run_config,
        )


def _r_to_pandas(r_df) -> pd.DataFrame:
    conv = ro.default_converter + pandas2ri.converter
    with localconverter(conv):
        return ro.conversion.get_conversion().rpy2py(r_df)


def _r_matrix_to_numpy(r_matrix) -> np.ndarray:
    conv = ro.default_converter + pandas2ri.converter
    with localconverter(conv):
        arr = np.asarray(ro.conversion.get_conversion().rpy2py(r_matrix))
    if arr.ndim == 1:
        dim = tuple(int(x) for x in r_matrix.dim)
        arr = arr.reshape(dim, order="F")
    return arr


@frozen
class _CovStruct:
    """LDSC covariance structure extracted from the R covstruc list."""

    S_LD: np.ndarray  # (k, k) genetic covariance
    V_LD: np.ndarray  # sampling covariance of vech(S_LD)
    I_LD: np.ndarray  # (k, k) LDSC intercepts


@frozen
class _SubtractionFrames:
    """The two per-factor result tables written by the task."""

    f_df: pd.DataFrame  # common factor F ~ SNP
    r_df: pd.DataFrame  # remainder factor R ~ SNP


def _extract_covstruc_arrays(covstruc_r) -> _CovStruct:
    return _CovStruct(
        S_LD=_r_matrix_to_numpy(covstruc_r.rx2("S")),
        V_LD=_r_matrix_to_numpy(covstruc_r.rx2("V")),
        I_LD=_r_matrix_to_numpy(covstruc_r.rx2("I")),
    )


def _make_result_df(
    snps_df: pd.DataFrame,
    est: np.ndarray,
    se_c: np.ndarray,
    z: np.ndarray,
    p: np.ndarray,
    n_eff: np.ndarray,
    fail: np.ndarray,
    lhs: str,
) -> pd.DataFrame:
    n = len(snps_df)
    return pd.DataFrame(
        {
            "SNP": snps_df["SNP"].to_numpy(),
            "CHR": snps_df.get("CHR", pd.Series([pd.NA] * n)).to_numpy(),
            "BP": snps_df.get("BP", pd.Series([pd.NA] * n)).to_numpy(),
            "MAF": snps_df["MAF"].to_numpy(),
            "A1": snps_df.get("A1", pd.Series([""] * n)).to_numpy(),
            "A2": snps_df.get("A2", pd.Series([""] * n)).to_numpy(),
            "lhs": [lhs] * n,
            "op": ["~"] * n,
            "rhs": ["SNP"] * n,
            "est": est,
            "se_c": se_c,
            "Z_Estimate": z,
            "Pval_Estimate": p,
            "N_eff": n_eff,
            "fail": fail,
        }
    )


def _run_python_subtraction(*, covstruc_r, snps_r) -> _SubtractionFrames:
    cov = _extract_covstruc_arrays(covstruc_r)
    snps_df = _r_to_pandas(snps_r)

    beta_cols = [c for c in snps_df.columns if str(c).startswith("beta.")]
    se_cols = [c for c in snps_df.columns if str(c).startswith("se.")]
    assert len(beta_cols) == len(se_cols) == 2

    beta_SNP = snps_df[beta_cols].to_numpy(dtype=float)
    SE_SNP = snps_df[se_cols].to_numpy(dtype=float)
    maf = snps_df["MAF"].to_numpy(dtype=float)
    varSNP = 2.0 * maf * (1.0 - maf)

    result = fit_gwas_by_subtraction(
        S_LD=cov.S_LD,
        V_LD=cov.V_LD,
        I_LD=cov.I_LD,
        beta_SNP=beta_SNP,
        SE_SNP=SE_SNP,
        varSNP=varSNP,
    )

    f_df = _make_result_df(
        snps_df,
        est=result.beta_F,
        se_c=result.se_c_F,
        z=result.z_F,
        p=result.p_F,
        n_eff=result.n_eff_F,
        fail=result.fail,
        lhs="F",
    )
    r_df = _make_result_df(
        snps_df,
        est=result.beta_R,
        se_c=result.se_c_R,
        z=result.z_R,
        p=result.p_R,
        n_eff=result.n_eff_R,
        fail=result.fail,
        lhs="R",
    )
    return _SubtractionFrames(f_df=f_df, r_df=r_df)
