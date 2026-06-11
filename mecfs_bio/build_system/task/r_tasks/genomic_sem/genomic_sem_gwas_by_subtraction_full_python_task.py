"""
GWAS-by-subtraction with the *entire* pipeline in Python -- no rpy2/R.

This is the full-Python sibling of `GenomicSEMGWASBySubtractionPythonTask`,
which still calls R for munge/ldsc/sumstats and only runs the Cholesky kernel
in numpy. Here every stage is Python:

    munge_sumstats   (polars)  -> per-trait .sumstats tables
    run_ldsc         (numpy)   -> S / V / I genetic covariance structure
    run_sumstats     (polars)  -> reference-aligned per-SNP betas/SEs
    fit_gwas_by_subtraction (numpy) -> per-SNP factor effects + delta-method SE

Each port has been validated against GenomicSEM to machine precision on real
data (see `_genomic_sem_ldsc`, `_genomic_sem_munge`, `_genomic_sem_sumstats`).
The win over the R-backed task is runtime: R munge+sumstats took ~12 min
combined on the user's laptop; the polars ports run in seconds.

Convention: the two traits are named explicitly (composite_trait_source = T1,
reference_trait_source = T2) to avoid ordering mistakes. Factor F is the genetic
factor common to both traits (defined by the reference trait T2); factor R is
the remainder unique to the composite trait T1, orthogonal to F. See
`_gwas_by_subtraction_kernel` for the model.

Output: two parquets under gwas_results/ -- one for the common factor (F~SNP)
and one for the remainder factor (R~SNP) -- in the same column layout as the
R-backed task, so downstream consumers are interchangeable.
"""

import gc
from pathlib import Path, PurePath

import numpy as np
import pandas as pd
import polars as pl
import structlog
from attrs import frozen

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
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_config import (
    GWAS_RESULTS_SUBDIR,
    LDSC_S_FILENAME,
    LDSC_S_STAND_FILENAME,
    LDSC_V_FILENAME,
    LINEAR_PROB,
    LOGISTIC,
    MUNGED_SUBDIR,
    OLS,
    SUBTRACTION_F_FILENAME,
    SUBTRACTION_R_FILENAME,
    GenomicSEMConfig,
    GenomicSEMGWASRunConfig,
    GenomicSEMGWASSumstatsSource,
    GenomicSEMSumstatsConfig,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_inputs import (
    build_munge_input_df,
    get_prevs,
    get_sample_size,
    ld_dir_with_genomic_sem_naming,
    read_dataframe_from_task,
    resolve_ld_path,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_ldsc import (
    LDSCResult,
    run_ldsc,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_munge import (
    munge_sumstats,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_sumstats import (
    SumstatsTrait,
    run_sumstats,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._gwas_by_subtraction_kernel import (
    fit_gwas_by_subtraction,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._subtraction_result import (
    SubtractionFrames,
    make_result_df,
)
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()

# The munged sumstats are written uncompressed (polars has no csv-gzip writer);
# run_ldsc reads them straight back via pandas, which is format-agnostic.
_MUNGED_SUFFIX = ".sumstats"


@frozen
class GenomicSEMGWASBySubtractionFullPythonTask(Task):
    """
    GWAS-by-subtraction with munge, LDSC, sumstats, and the Cholesky kernel all
    computed in Python (no rpy2/R at execution time).

    The two traits are named explicitly to remove any ordering ambiguity:
    - composite_trait_source: T1 (loads on both factors, e.g. educational
      attainment); factor R is its residual.
    - reference_trait_source: T2 (a pure indicator of the common factor F,
      e.g. cognitive performance).

    Output: two parquets under gwas_results/ -- one for the common factor
    (F~SNP) and one for the remainder factor (R~SNP).
    """

    meta: Meta
    composite_trait_source: GenomicSEMGWASSumstatsSource
    reference_trait_source: GenomicSEMGWASSumstatsSource
    ld_ref_task: Task
    hapmap_snps_task: Task
    sumstats_ref_task: Task
    munge_config: GenomicSEMConfig = GenomicSEMConfig()
    sumstats_config: GenomicSEMSumstatsConfig = GenomicSEMSumstatsConfig()
    run_config: GenomicSEMGWASRunConfig = GenomicSEMGWASRunConfig()

    def __attrs_post_init__(self):
        assert self.composite_trait_source.alias != self.reference_trait_source.alias, (
            "composite and reference trait aliases must differ"
        )

    @property
    def _ordered_sources(
        self,
    ) -> tuple[GenomicSEMGWASSumstatsSource, GenomicSEMGWASSumstatsSource]:
        # Order defines the trait axis for ldsc/sumstats and the kernel:
        # index 0 = T1 (composite), index 1 = T2 (reference).
        return tuple([self.composite_trait_source, self.reference_trait_source])

    @property
    def deps(self) -> list[Task]:
        return [
            self.ld_ref_task,
            self.hapmap_snps_task,
            self.sumstats_ref_task,
            self.composite_trait_source.task,
            self.reference_trait_source.task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        ld_path = resolve_ld_path(self.ld_ref_task, fetch)
        ref_hm3 = read_dataframe_from_task(self.hapmap_snps_task, fetch)
        ref_1kg = read_dataframe_from_task(self.sumstats_ref_task, fetch)

        with ld_dir_with_genomic_sem_naming(
            ld_path, self.munge_config.ld_file_filename_pattern
        ) as effective_ld:
            inputs = _prepare_python_inputs(
                sources=self._ordered_sources,
                ld_dir=effective_ld,
                ref_hm3=ref_hm3,
                ref_1kg=ref_1kg,
                munge_config=self.munge_config,
                sumstats_config=self.sumstats_config,
                fetch=fetch,
                scratch_dir=scratch_dir,
            )
            logger.debug("Running Python GWAS-by-subtraction kernel")
            frames = _build_subtraction_frames(inputs)
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
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        composite_trait_source: GenomicSEMGWASSumstatsSource,
        reference_trait_source: GenomicSEMGWASSumstatsSource,
        ld_ref_task: Task,
        hapmap_snps_task: Task,
        sumstats_ref_task: Task,
        munge_config: GenomicSEMConfig = GenomicSEMConfig(),
        sumstats_config: GenomicSEMSumstatsConfig = GenomicSEMSumstatsConfig(),
        run_config: GenomicSEMGWASRunConfig = GenomicSEMGWASRunConfig(),
    ) -> "GenomicSEMGWASBySubtractionFullPythonTask":
        meta = ResultDirectoryMeta(
            id=AssetId(asset_id),
            trait=MULTI_TRAIT,
            project="genomic_sem",
            sub_dir=PurePath("analysis"),
        )
        return cls(
            meta=meta,
            composite_trait_source=composite_trait_source,
            reference_trait_source=reference_trait_source,
            ld_ref_task=ld_ref_task,
            hapmap_snps_task=hapmap_snps_task,
            sumstats_ref_task=sumstats_ref_task,
            munge_config=munge_config,
            sumstats_config=sumstats_config,
            run_config=run_config,
        )


@frozen
class _PythonGWASInputs:
    """Everything the kernel needs, produced by the Python munge/ldsc/sumstats."""

    S_LD: np.ndarray  # (2, 2) genetic covariance
    V_LD: np.ndarray  # (3, 3) sampling covariance of vech(S_LD)
    I_LD: np.ndarray  # (2, 2) LDSC intercepts
    snps_df: pl.DataFrame  # reference-aligned per-SNP betas/SEs (run_sumstats)
    trait_names: list[str]  # T1, T2 order (defines the beta./se. column order)


def sumstats_trait(
    gwas_src: GenomicSEMGWASSumstatsSource, df: pl.DataFrame, n: float
) -> SumstatsTrait:
    """Map a source's GWAS method onto the SumstatsTrait standardization flag."""
    method = gwas_src.gwas_method
    return SumstatsTrait(
        df=df,
        name=gwas_src.alias,
        n=n,
        ols=method == OLS,
        se_logit=method == LOGISTIC,
        linprob=method == LINEAR_PROB,
    )


def _save_python_ldsc_outputs(result: LDSCResult, scratch_dir: Path) -> None:
    """Write S/V/S_Stand to CSV for inspection, matching the R task's layout."""
    pd.DataFrame(result.S).to_csv(scratch_dir / LDSC_S_FILENAME, index=False)
    pd.DataFrame(result.V).to_csv(scratch_dir / LDSC_V_FILENAME, index=False)
    if result.S_Stand is not None:
        pd.DataFrame(result.S_Stand).to_csv(
            scratch_dir / LDSC_S_STAND_FILENAME, index=False
        )


def _prepare_python_inputs(
    *,
    sources: tuple[GenomicSEMGWASSumstatsSource, GenomicSEMGWASSumstatsSource],
    ld_dir: Path,
    ref_hm3: pl.DataFrame,
    ref_1kg: pl.DataFrame,
    munge_config: GenomicSEMConfig,
    sumstats_config: GenomicSEMSumstatsConfig,
    fetch: Fetch,
    scratch_dir: Path,
) -> _PythonGWASInputs:
    """
    Run munge -> ldsc -> sumstats entirely in Python and return the kernel
    inputs. ``sources`` is ordered [composite (T1), reference (T2)]; ``ref_hm3``
    and ``ref_1kg`` are the HapMap3 and 1000G reference panels. Side effects:
    writes the munged sumstats and LDSC CSVs into scratch_dir.
    """
    munged_dir = scratch_dir / MUNGED_SUBDIR
    munged_dir.mkdir(parents=True, exist_ok=True)

    trait_names: list[str] = []
    sample_sizes: list[float] = []
    sample_prevs: list[float] = []
    population_prevs: list[float] = []
    munged_paths: list[Path] = []
    sumstats_traits: list[SumstatsTrait] = []

    for gwas_src in sources:
        df = build_munge_input_df(gwas_src.source, fetch)
        name = gwas_src.alias
        n = get_sample_size(gwas_src.source)
        samp_prev, pop_prev = get_prevs(gwas_src.source.sample_info)

        logger.debug(f"Munging '{name}' ({df.height} variants)")
        munged = munge_sumstats(
            df,
            ref_hm3,
            n=n,
            info_filter=munge_config.info_filter,
            maf_filter=munge_config.maf_filter,
        )
        munged_path = munged_dir / f"{name}{_MUNGED_SUFFIX}"
        munged.write_csv(munged_path, separator="\t")

        trait_names.append(name)
        sample_sizes.append(n)
        sample_prevs.append(samp_prev)
        population_prevs.append(pop_prev)
        munged_paths.append(munged_path)
        sumstats_traits.append(sumstats_trait(gwas_src, df, n))

    logger.debug("Running Python LDSC")
    ldsc_result = run_ldsc(
        munged_paths=munged_paths,
        ld_dir=ld_dir,
        sample_prev=sample_prevs,
        population_prev=population_prevs,
    )
    _save_python_ldsc_outputs(ldsc_result, scratch_dir)

    logger.debug("Running Python sumstats")
    snps_df = run_sumstats(
        sumstats_traits,
        ref_1kg,
        maf_filter=sumstats_config.maf_filter,
        info_filter=sumstats_config.info_filter,
        ambig=sumstats_config.ambig,
    )
    logger.debug(f"sumstats produced {snps_df.height} aligned SNPs")

    return _PythonGWASInputs(
        S_LD=ldsc_result.S,
        V_LD=ldsc_result.V,
        I_LD=ldsc_result.I,
        snps_df=snps_df,
        trait_names=trait_names,
    )


def _build_subtraction_frames(inputs: _PythonGWASInputs) -> SubtractionFrames:
    """Run the Cholesky kernel and package the F/R per-factor result tables."""
    snps_pd = inputs.snps_df.to_pandas()
    # Columns are ordered by trait: col 0 = T1 (composite), col 1 = T2 (reference),
    # matching the kernel's beta_SNP convention.
    beta_SNP = np.column_stack(
        [snps_pd[f"beta.{name}"].to_numpy(dtype=float) for name in inputs.trait_names]
    )
    SE_SNP = np.column_stack(
        [snps_pd[f"se.{name}"].to_numpy(dtype=float) for name in inputs.trait_names]
    )
    maf = snps_pd["MAF"].to_numpy(dtype=float)
    varSNP = 2.0 * maf * (1.0 - maf)

    result = fit_gwas_by_subtraction(
        S_LD=inputs.S_LD,
        V_LD=inputs.V_LD,
        I_LD=inputs.I_LD,
        beta_SNP=beta_SNP,
        SE_SNP=SE_SNP,
        varSNP=varSNP,
    )

    f_df = make_result_df(
        snps_pd,
        est=result.beta_F,
        se_c=result.se_c_F,
        z=result.z_F,
        p=result.p_F,
        n_eff=result.n_eff_F,
        fail=result.fail,
        lhs="F",
    )
    r_df = make_result_df(
        snps_pd,
        est=result.beta_R,
        se_c=result.se_c_R,
        z=result.z_R,
        p=result.p_R,
        n_eff=result.n_eff_R,
        fail=result.fail,
        lhs="R",
    )
    return SubtractionFrames(f_df=f_df, r_df=r_df)
