"""
GenomicSEM GWAS-extension workflows: common factor GWAS and multivariate
(user) GWAS. Both add SNP effects to the GenomicSEM model and produce per-SNP
summary statistics for the latent factors.

Workflow:
  1. munge() — produces .sumstats.gz for LDSC (shared with GenomicSEMTask).
  2. ldsc() — produces covstruc.
  3. sumstats() — aligns per-trait stats to a reference panel.
  4. commonfactorGWAS(covstruc, SNPs) or userGWAS(covstruc, SNPs, model).

References:
  - Common Factor GWAS:  https://github.com/GenomicSEM/GenomicSEM/wiki/4.-Common-Factor-GWAS
  - Multivariate GWAS:   https://github.com/GenomicSEM/GenomicSEM/wiki/5.-Multivariate-GWAS

GenomicSEM is an R library, so it is accessed through Python via rpy2.
"""

import gc
import re
import tempfile
from pathlib import Path, PurePath
from typing import Final, Literal, Sequence

import pandas as pd
import rpy2.robjects as ro
import structlog
from attrs import frozen
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    MULTI_TRAIT,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_task import (
    LAVAAN_MODEL_FILENAME,
    LDSC_LOG_PREFIX,
    MUNGED_SUBDIR,
    GenomicSEMConfig,
    GenomicSEMSumstatsSource,
    _get_prevs,
    _get_sample_size,
    _run_ldsc,
    _run_munge,
    _save_ldsc_outputs,
    _write_munge_input,
)
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()


GWAS_RESULTS_SUBDIR = "gwas_results"
COMMON_FACTOR_GWAS_FILENAME = "common_factor.parquet"


# How the source GWAS was estimated. Controls the per-trait flag (one of
# se.logit / OLS / linprob) passed to GenomicSEM::sumstats.
OLS: Final = "ols"
LOGISTIC: Final = "logistic"
LINEAR_PROB: Final = "linear_prob"

# `ty` requires the type arguments to Literal to be inline literals rather
# than Final-typed names, so the strings are repeated here.
GWASMethod = Literal["ols", "logistic", "linear_prob"]


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


@frozen
class GenomicSEMCommonFactorGWASTask(Task):
    """
    Run GenomicSEM common factor GWAS: munge → ldsc → sumstats →
    common factor GWAS. Output is a parquet of per-SNP common factor effects.
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
            covstruc, snps = _prepare_gwas_inputs(
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
            logger.debug("Running GenomicSEM::commonfactorGWAS")
            result = _run_common_factor_gwas(
                gsem=gsem, covstruc=covstruc, snps=snps, config=self.run_config
            )
            _save_common_factor_gwas_output(result, scratch_dir)

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
    ) -> "GenomicSEMCommonFactorGWASTask":
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


@frozen
class GenomicSEMUserGWASTask(Task):
    """
    Run GenomicSEM multivariate (user) GWAS: munge → ldsc → sumstats →
    userGWAS with a user-supplied lavaan model containing SNP effects.

    sub_components selects which model components are extracted per SNP and
    written to disk, e.g. ["F1~SNP", "F2~SNP"]. One parquet per component is
    written to gwas_results/.
    """

    meta: Meta
    sources: Sequence[GenomicSEMGWASSumstatsSource]
    ld_ref_task: Task
    hapmap_snps_task: Task
    sumstats_ref_task: Task
    factor_model: str
    sub_components: Sequence[str]
    munge_config: GenomicSEMConfig = GenomicSEMConfig()
    sumstats_config: GenomicSEMSumstatsConfig = GenomicSEMSumstatsConfig()
    run_config: GenomicSEMGWASRunConfig = GenomicSEMGWASRunConfig()
    fix_measurement: bool = True
    q_snp: bool = False
    std_lv: bool = False

    def __attrs_post_init__(self):
        _validate_sources(self.sources)
        assert len(self.sub_components) >= 1, (
            "userGWAS requires at least one sub_components entry"
        )
        sanitized = [_sanitize_component_name(c) for c in self.sub_components]
        assert len(set(sanitized)) == len(sanitized), (
            f"sub_components map to colliding filenames: {sanitized}"
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
            covstruc, snps = _prepare_gwas_inputs(
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
            (scratch_dir / LAVAAN_MODEL_FILENAME).write_text(self.factor_model)
            logger.debug("Running GenomicSEM::userGWAS")
            result = _run_user_gwas(
                gsem=gsem,
                covstruc=covstruc,
                snps=snps,
                model=self.factor_model,
                sub_components=self.sub_components,
                fix_measurement=self.fix_measurement,
                q_snp=self.q_snp,
                std_lv=self.std_lv,
                config=self.run_config,
            )
            _save_user_gwas_outputs(
                result=result,
                sub_components=self.sub_components,
                scratch_dir=scratch_dir,
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
        factor_model: str,
        sub_components: Sequence[str],
        munge_config: GenomicSEMConfig = GenomicSEMConfig(),
        sumstats_config: GenomicSEMSumstatsConfig = GenomicSEMSumstatsConfig(),
        run_config: GenomicSEMGWASRunConfig = GenomicSEMGWASRunConfig(),
        fix_measurement: bool = True,
        q_snp: bool = False,
        std_lv: bool = False,
    ) -> "GenomicSEMUserGWASTask":
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
            factor_model=factor_model,
            sub_components=sub_components,
            munge_config=munge_config,
            sumstats_config=sumstats_config,
            run_config=run_config,
            fix_measurement=fix_measurement,
            q_snp=q_snp,
            std_lv=std_lv,
        )


def _validate_sources(sources: Sequence[GenomicSEMGWASSumstatsSource]) -> None:
    assert len(sources) >= 2, "GenomicSEM GWAS requires at least two traits"
    aliases = [s.alias for s in sources]
    assert len(set(aliases)) == len(aliases), (
        f"Trait aliases must be unique. Got: {aliases}"
    )


def _gwas_method_flags(
    sources: Sequence[GenomicSEMGWASSumstatsSource],
) -> tuple[list[bool], list[bool], list[bool]]:
    """
    GenomicSEM::sumstats expects a (se.logit, OLS, linprob) BoolVector triple.
    Exactly one is TRUE per trait.
    """
    se_logit: list[bool] = []
    ols: list[bool] = []
    linprob: list[bool] = []
    for s in sources:
        se_logit.append(s.gwas_method == LOGISTIC)
        ols.append(s.gwas_method == OLS)
        linprob.append(s.gwas_method == LINEAR_PROB)
    return se_logit, ols, linprob


def _sanitize_component_name(name: str) -> str:
    """Map a lavaan sub-component (e.g. 'F1~SNP') to a filename stem ('F1_SNP')."""
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


def _resolve_ld_path(
    ld_ref_task: Task, fetch: Fetch, munge_config: GenomicSEMConfig
) -> str:
    asset = fetch(ld_ref_task.asset_id)
    assert isinstance(asset, DirectoryAsset)
    return str(asset.path) + munge_config.ld_file_filename_pattern


def _resolve_file_path(task: Task, fetch: Fetch) -> str:
    asset = fetch(task.asset_id)
    assert isinstance(asset, FileAsset)
    return str(asset.path)


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


def _run_common_factor_gwas(*, gsem, covstruc, snps, config: GenomicSEMGWASRunConfig):
    kwargs = dict(
        covstruc=covstruc,
        SNPs=snps,
        estimation=config.estimation,
        parallel=config.parallel,
        GC=config.gc_correction,
        toler=config.toler,
        SNPSE=config.snp_se,
        smooth_check=config.smooth_check,
    )
    if config.cores is not None:
        kwargs["cores"] = config.cores
    return gsem.commonfactorGWAS(**kwargs)


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


def _r_to_pandas(r_df) -> pd.DataFrame:
    conv = ro.default_converter + pandas2ri.converter
    with localconverter(conv):
        return ro.conversion.get_conversion().rpy2py(r_df)


def _save_common_factor_gwas_output(result, scratch_dir: Path) -> Path:
    out_dir = scratch_dir / GWAS_RESULTS_SUBDIR
    out_dir.mkdir(parents=True, exist_ok=True)
    df = _r_to_pandas(result)
    out_path = out_dir / COMMON_FACTOR_GWAS_FILENAME
    df.to_parquet(out_path, index=False)
    logger.debug(f"Wrote common factor GWAS sumstats to {out_path} ({len(df)} rows)")
    return out_path


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
