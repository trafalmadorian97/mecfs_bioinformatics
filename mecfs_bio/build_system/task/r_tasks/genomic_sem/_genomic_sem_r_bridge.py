"""
rpy2 bridge: the R-calling helpers shared across the R-backed GenomicSEM
tasks (munge, ldsc, sumstats, userGWAS, and the R<->Python conversions).

This is the only shared module that imports rpy2. The full-Python tasks must
not depend on it; they use the polars/numpy ports in `_genomic_sem_munge`,
`_genomic_sem_ldsc`, and `_genomic_sem_sumstats` instead.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import rpy2.robjects as ro
import structlog
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_config import (
    GWAS_RESULTS_SUBDIR,
    LDSC_LOG_PREFIX,
    LDSC_S_FILENAME,
    LDSC_S_STAND_FILENAME,
    LDSC_V_FILENAME,
    MODEL_FIT_FILENAME,
    MODEL_RESULTS_FILENAME,
    MUNGED_SUBDIR,
    GenomicSEMConfig,
    GenomicSEMGWASRunConfig,
    GenomicSEMGWASSumstatsSource,
    GenomicSEMSumstatsConfig,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_inputs import (
    _get_prevs,
    _get_sample_size,
    _gwas_method_flags,
    _ld_dir_with_genomic_sem_naming,
    _sanitize_component_name,
    _write_munge_input,
)

logger = structlog.get_logger()


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


def _run_munge(
    gsem,
    input_files: list[str],
    hapmap_path: str,
    trait_names: list[str],
    sample_sizes: list[float],
    output_dir: Path,
    info_filter: float,
    maf_filter: float,
) -> None:
    """
    Run GenomicSEM::munge. munge writes <trait_name>.sumstats.gz files to the
    current working directory together with a <trait_name>.log file. We chdir
    into output_dir for the call to keep all generated artifacts inside the
    target output directory.
    """
    base = importr("base")
    cwd_before = str(base.getwd()[0])
    try:
        base.setwd(str(output_dir))
        gsem.munge(
            files=ro.StrVector(input_files),
            hm3=hapmap_path,
            trait_names=ro.StrVector(trait_names),
            N=ro.FloatVector(sample_sizes),
            info_filter=info_filter,
            maf_filter=maf_filter,
        )
    finally:
        base.setwd(cwd_before)


def _run_ldsc(
    gsem,
    munged_paths: list[str],
    trait_names: list[str],
    sample_prevs: list[float],
    population_prevs: list[float],
    ld_path: str,
    ld_file_basename_prefix: str,
    ldsc_log_prefix: str,
):
    """
    Run GenomicSEM::ldsc to obtain a genetic covariance matrix. The same path
    is used for both ld and wld since GenomicSEM expects pre-built LD scores
    in standard form. ldsc appends "_ldsc.log" to the supplied prefix.

    The LDSC function returns an object with the following attributes:

        LDSCoutput$S is the covariance matrix (on the liability scale for
            case/control designs).
        LDSCoutput$V which is the sampling covariance matrix in the format
            expected by lavaan.
        LDSCoutput$I is the matrix of LDSC intercepts and cross-trait (i.e.,
            bivariate) intercepts.
        LDSCoutput$N contains the sample sizes (N) for the heritabilities and
            sqrt(N1N2) for the co-heritabilities.
        LDSCoutput$m is the number of SNPs used to construct the LD score.
    """
    with _ld_dir_with_genomic_sem_naming(
        ld_path, ld_file_basename_prefix
    ) as effective_ld:
        return gsem.ldsc(
            traits=ro.StrVector(munged_paths),
            sample_prev=ro.FloatVector(sample_prevs),
            population_prev=ro.FloatVector(population_prevs),
            ld=effective_ld,
            wld=effective_ld,
            trait_names=ro.StrVector(trait_names),
            ldsc_log=ldsc_log_prefix,
            stand=True,
        )


def _save_ldsc_outputs(covstruc, scratch_dir: Path) -> None:
    """
    Convert the LDSC output (an R named list) into CSV files so the outputs are
    inspectable from outside R.
    """
    conv = ro.default_converter + pandas2ri.converter

    s_matrix = covstruc.rx2("S")
    v_matrix = covstruc.rx2("V")
    with localconverter(conv):
        s_df = pd.DataFrame(ro.conversion.get_conversion().rpy2py(s_matrix))
        v_df = pd.DataFrame(ro.conversion.get_conversion().rpy2py(v_matrix))

    s_df.to_csv(scratch_dir / LDSC_S_FILENAME, index=False)
    v_df.to_csv(scratch_dir / LDSC_V_FILENAME, index=False)

    if "S_Stand" in tuple(covstruc.names):
        s_stand = covstruc.rx2("S_Stand")
        with localconverter(conv):
            s_stand_df = pd.DataFrame(ro.conversion.get_conversion().rpy2py(s_stand))
        s_stand_df.to_csv(scratch_dir / LDSC_S_STAND_FILENAME, index=False)


def _save_model_outputs(model_result, scratch_dir: Path) -> None:
    """
    GenomicSEM::usermodel returns a named R list with $modelfit (fit indices)
    and $results (parameter estimates).
    """
    conv = ro.default_converter + pandas2ri.converter
    fit_r = model_result.rx2("modelfit")
    results_r = model_result.rx2("results")
    with localconverter(conv):
        fit_df: pd.DataFrame = ro.conversion.get_conversion().rpy2py(fit_r)
        results_df: pd.DataFrame = ro.conversion.get_conversion().rpy2py(results_r)
    logger.debug(f"GenomicSEM model fit:\n{fit_df}")
    logger.debug(f"GenomicSEM model results:\n{results_df}")
    fit_df.to_csv(scratch_dir / MODEL_FIT_FILENAME, index=False)
    results_df.to_csv(scratch_dir / MODEL_RESULTS_FILENAME, index=False)


def _prepare_gwas_inputs(
    *,
    gsem,
    sources: Sequence[GenomicSEMGWASSumstatsSource],
    ld_path: str,
    hapmap_path: str,
    sumstats_ref_path: str,
    munge_config: GenomicSEMConfig,
    sumstats_config: GenomicSEMSumstatsConfig,
    fetch: Fetch,
    scratch_dir: Path,
    tmp_dir: Path,
):
    """
    Run munge → ldsc → sumstats and return (covstruc, snps_r_df).

    Side effects: writes munged sumstats and LDSC outputs into scratch_dir.
    """
    munged_dir = scratch_dir / MUNGED_SUBDIR
    munged_dir.mkdir(parents=True, exist_ok=True)

    input_files: list[str] = []
    trait_names: list[str] = []
    sample_sizes: list[float] = []
    sample_prevs: list[float] = []
    population_prevs: list[float] = []
    for gwas_src in sources:
        input_path = _write_munge_input(
            source=gwas_src.source, fetch=fetch, tmp_dir=tmp_dir
        )
        input_files.append(str(input_path))
        trait_names.append(gwas_src.alias)
        sample_sizes.append(_get_sample_size(gwas_src.source))
        samp_prev, pop_prev = _get_prevs(gwas_src.source.sample_info)
        sample_prevs.append(samp_prev)
        population_prevs.append(pop_prev)

    logger.debug(f"Running GenomicSEM::munge on {len(input_files)} files")
    _run_munge(
        gsem=gsem,
        input_files=input_files,
        hapmap_path=hapmap_path,
        trait_names=trait_names,
        sample_sizes=sample_sizes,
        output_dir=munged_dir,
        info_filter=munge_config.info_filter,
        maf_filter=munge_config.maf_filter,
    )
    munged_paths = [str(munged_dir / f"{name}.sumstats.gz") for name in trait_names]
    for path in munged_paths:
        assert Path(path).is_file(), f"Munged sumstats not found at {path}"

    logger.debug("Running GenomicSEM::ldsc")
    covstruc = _run_ldsc(
        gsem=gsem,
        munged_paths=munged_paths,
        trait_names=trait_names,
        sample_prevs=sample_prevs,
        population_prevs=population_prevs,
        ld_path=ld_path,
        ld_file_basename_prefix=munge_config.ld_file_filename_pattern,
        ldsc_log_prefix=str(scratch_dir / LDSC_LOG_PREFIX),
    )
    _save_ldsc_outputs(covstruc=covstruc, scratch_dir=scratch_dir)

    logger.debug("Running GenomicSEM::sumstats")
    se_logit_flags, ols_flags, linprob_flags = _gwas_method_flags(sources)
    snps = _run_sumstats(
        gsem=gsem,
        input_files=input_files,
        ref_path=sumstats_ref_path,
        trait_names=trait_names,
        se_logit_flags=se_logit_flags,
        ols_flags=ols_flags,
        linprob_flags=linprob_flags,
        sample_sizes=sample_sizes,
        config=sumstats_config,
    )
    return covstruc, snps


def _run_sumstats(
    *,
    gsem,
    input_files: list[str],
    ref_path: str,
    trait_names: list[str],
    se_logit_flags: list[bool],
    ols_flags: list[bool],
    linprob_flags: list[bool],
    sample_sizes: list[float],
    config: GenomicSEMSumstatsConfig,
):
    """
    From the wiki: https://github.com/GenomicSEM/GenomicSEM/wiki/4.-Common-Factor-GWAS

    The sumstats function aligns each trait's per-SNP effects to a common
    reference allele and rescales coefficients/SEs to unit-variance phenotypes,
    combining them into a single listwise-deleted file. Variables must be
    listed in the same order as for the ldsc function.
    """
    kwargs = dict(
        files=ro.StrVector(input_files),
        ref=ref_path,
        trait_names=ro.StrVector(trait_names),
        se_logit=ro.BoolVector(se_logit_flags),
        OLS=ro.BoolVector(ols_flags),
        linprob=ro.BoolVector(linprob_flags),
        N=ro.FloatVector(sample_sizes),
        info_filter=config.info_filter,
        maf_filter=config.maf_filter,
        keep_indel=config.keep_indel,
        parallel=config.parallel,
        ambig=config.ambig,
    )
    if config.cores is not None:
        kwargs["cores"] = config.cores
    return gsem.sumstats(**kwargs)


def _run_user_gwas(
    *,
    gsem,
    covstruc,
    snps,
    model: str,
    sub_components: Sequence[str],
    fix_measurement: bool,
    q_snp: bool,
    std_lv: bool,
    config: GenomicSEMGWASRunConfig,
):
    """
    See section 5 of https://github.com/GenomicSEM/GenomicSEM/wiki/5.-Multivariate-GWAS
    """
    kwargs = dict(
        covstruc=covstruc,
        SNPs=snps,
        model=model,
        estimation=config.estimation,
        sub=ro.StrVector(list(sub_components)),
        parallel=config.parallel,
        GC=config.gc_correction,
        toler=config.toler,
        SNPSE=config.snp_se,
        smooth_check=config.smooth_check,
        fix_measurement=fix_measurement,
        Q_SNP=q_snp,
        std_lv=std_lv,
    )
    if config.cores is not None:
        kwargs["cores"] = config.cores
    return gsem.userGWAS(**kwargs)


def _save_user_gwas_outputs(
    *, result, sub_components: Sequence[str], scratch_dir: Path
) -> list[Path]:
    """
    userGWAS with `sub` returns an R list of dataframes, one per sub component
    in the same order. Write each as parquet under gwas_results/.
    """
    out_dir = scratch_dir / GWAS_RESULTS_SUBDIR
    out_dir.mkdir(parents=True, exist_ok=True)
    out_paths: list[Path] = []
    for i, component in enumerate(sub_components, start=1):
        r_df = result.rx2(i)
        df = _r_to_pandas(r_df)
        filename = f"{_sanitize_component_name(component)}.parquet"
        out_path = out_dir / filename
        df.to_parquet(out_path, index=False)
        logger.debug(
            f"Wrote userGWAS component '{component}' to {out_path} ({len(df)} rows)"
        )
        out_paths.append(out_path)
    return out_paths
