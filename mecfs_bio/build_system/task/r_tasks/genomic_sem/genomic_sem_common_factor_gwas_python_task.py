"""
GenomicSEM common factor GWAS via munge+ldsc+sumstats in R, but with the
per-SNP commonfactorGWAS loop reimplemented in vectorized numpy/scipy.

This is the Python sibling of `GenomicSEMCommonFactorGWASTask`. The setup
phase (munge + multivariable LDSC + sumstats reference alignment) is shared
unchanged, but the per-SNP DWLS fit + sandwich SE happens in
`_common_factor_kernel.fit_common_factor_gwas` instead of in lavaan. This
removes lavaan's per-call overhead, which on the user's WSL2 laptop turned
~20 ms per SNP into ~3 us per SNP -- ~7000x speedup at the per-SNP level.

v1 limitations (asserted at runtime):
- k = 2 traits only
- GC = "standard"
- DWLS only (no ML)
- No Q heterogeneity statistic (Q, Q_df, Q_pval columns set to NaN)
- No nearPD smoothing of V_Full or S_Fullrun (non-PD SNPs flagged warning=1)
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
from mecfs_bio.build_system.task.r_tasks.genomic_sem._common_factor_kernel import (
    fit_common_factor_gwas,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_task import (
    GenomicSEMConfig,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_user_gwas_task import (
    COMMON_FACTOR_GWAS_FILENAME,
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


@frozen
class GenomicSEMCommonFactorGWASPythonTask(Task):
    """
    Common-factor GWAS where the per-SNP loop runs in numpy instead of lavaan.

    Field signature mirrors `GenomicSEMCommonFactorGWASTask`. Only k=2 is
    supported in v1 -- the kernel will raise if more traits are passed.
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
            f"Python kernel v1 only supports k=2 traits; got {len(self.sources)}"
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
            logger.debug("Running Python commonfactorGWAS kernel")
            results_df = _run_python_common_factor_gwas(
                covstruc_r=covstruc_r,
                snps_r=snps_r,
                varSNPSE2=_resolve_snp_se(self.run_config),
            )
            out_dir = scratch_dir / GWAS_RESULTS_SUBDIR
            out_dir.mkdir(parents=True, exist_ok=True)
            out_path = out_dir / COMMON_FACTOR_GWAS_FILENAME
            results_df.to_parquet(out_path, index=False)
            logger.debug(
                f"Wrote common factor GWAS sumstats to {out_path} "
                f"({len(results_df)} rows; {int(results_df['fail'].sum())} fails)"
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
    ) -> "GenomicSEMCommonFactorGWASPythonTask":
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


def _resolve_snp_se(run_config: GenomicSEMGWASRunConfig) -> float:
    """
    GenomicSEM's SNPSE arg defaults to FALSE, which the R code interprets as
    "use (5e-4)^2". Match that behavior here.
    """
    if run_config.snp_se is False or run_config.snp_se is None:
        return (5e-4) ** 2
    return float(run_config.snp_se)


def _r_to_pandas(r_df) -> pd.DataFrame:
    conv = ro.default_converter + pandas2ri.converter
    with localconverter(conv):
        return ro.conversion.get_conversion().rpy2py(r_df)


def _r_matrix_to_numpy(r_matrix) -> np.ndarray:
    """Convert an R numeric matrix to a (rows, cols) numpy array."""
    conv = ro.default_converter + pandas2ri.converter
    with localconverter(conv):
        arr = np.asarray(ro.conversion.get_conversion().rpy2py(r_matrix))
    if arr.ndim == 1:
        # R sometimes returns 1-d when conversion drops the matrix attr;
        # reshape using the original dim attribute.
        dim = tuple(int(x) for x in r_matrix.dim)
        arr = arr.reshape(dim, order="F")
    return arr


def _extract_covstruc_arrays(
    covstruc_r,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (S_LD, V_LD, I_LD) from GenomicSEM's ldsc output."""
    S_LD = _r_matrix_to_numpy(covstruc_r.rx2("S"))
    V_LD = _r_matrix_to_numpy(covstruc_r.rx2("V"))
    I_LD = _r_matrix_to_numpy(covstruc_r.rx2("I"))
    return S_LD, V_LD, I_LD


def _run_python_common_factor_gwas(
    *, covstruc_r, snps_r, varSNPSE2: float
) -> pd.DataFrame:
    """
    Convert the R covstruc + SNPs dataframe to numpy, run the kernel, and
    package the results in the column layout downstream code expects.
    """
    S_LD, V_LD, I_LD = _extract_covstruc_arrays(covstruc_r)
    snps_df = _r_to_pandas(snps_r)

    # Sumstats() emits one beta.<trait> and se.<trait> per trait, in the order
    # the LDSC was run. Discover them by name to remain robust to extra
    # columns.
    beta_cols = [c for c in snps_df.columns if str(c).startswith("beta.")]
    se_cols = [c for c in snps_df.columns if str(c).startswith("se.")]
    assert len(beta_cols) == len(se_cols) == 2, (
        "Python kernel v1 only supports k=2 traits; sumstats produced "
        f"beta cols {beta_cols} and se cols {se_cols}"
    )

    beta_SNP = snps_df[beta_cols].to_numpy(dtype=float)
    SE_SNP = snps_df[se_cols].to_numpy(dtype=float)
    maf = snps_df["MAF"].to_numpy(dtype=float)
    varSNP = 2.0 * maf * (1.0 - maf)

    result = fit_common_factor_gwas(
        S_LD=S_LD,
        V_LD=V_LD,
        I_LD=I_LD,
        beta_SNP=beta_SNP,
        SE_SNP=SE_SNP,
        varSNP=varSNP,
        varSNPSE2=varSNPSE2,
    )

    # Match the column layout of R commonfactorGWAS so downstream code can
    # consume either output. Columns we don't compute (Q, lavaan-default se,
    # N_estimate) are NaN; lhs/op/rhs are the constant model labels.
    n = len(snps_df)
    out = pd.DataFrame(
        {
            "SNP": snps_df["SNP"].to_numpy(),
            "CHR": snps_df.get("CHR", pd.Series([pd.NA] * n)).to_numpy(),
            "BP": snps_df.get("BP", pd.Series([pd.NA] * n)).to_numpy(),
            "MAF": maf,
            "A1": snps_df.get("A1", pd.Series([""] * n)).to_numpy(),
            "A2": snps_df.get("A2", pd.Series([""] * n)).to_numpy(),
            "lhs": ["F1"] * n,
            "op": ["~"] * n,
            "rhs": ["SNP"] * n,
            "est": result.est,
            "se": np.full(n, np.nan, dtype=float),  # lavaan-default SE; not computed
            "se_c": result.se_c,
            "Q": np.full(n, np.nan, dtype=float),
            "Q_df": np.full(n, np.nan, dtype=float),
            "Q_pval": np.full(n, np.nan, dtype=float),
            "Z_Estimate": result.z,
            "Pval_Estimate": result.p,
            "fail": result.fail,
            "warning": result.warning,
            "n_iter": result.n_iter,
        }
    )
    return out
