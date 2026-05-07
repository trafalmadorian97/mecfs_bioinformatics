"""
Run a basic GenomicSEM workflow: munge, multivariable LDSC, and a user-specified
factor model expressed in lavaan syntax.

See: https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Genome%E2%80%90wide-Models

GenomicSEM is an R library, so it is accessed through Python via rpy2.
"""

import gc
import tempfile
from contextlib import contextmanager
from pathlib import Path, PurePath
from typing import Iterator, Literal, Sequence

import pandas as pd
import polars as pl
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
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    MULTI_TRAIT,
    BinaryPhenotypeSampleInfo,
    PhenotypeInfo,
    QuantPhenotype,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_P_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

logger = structlog.get_logger()


# Column names expected by GenomicSEM::munge by default.
MUNGE_SNP_COL = "SNP"
MUNGE_A1_COL = "A1"
MUNGE_A2_COL = "A2"
MUNGE_EFFECT_COL = "effect"
MUNGE_SE_COL = "SE"
MUNGE_P_COL = "P"
MUNGE_N_COL = "N"
MUNGE_MAF_COL = "MAF"

# Output structure produced by execute()
MUNGED_SUBDIR = "munged"
LDSC_LOG_PREFIX = "genomic_sem"  # GenomicSEM::ldsc appends "_ldsc.log" to this
LDSC_LOG_FILENAME = f"{LDSC_LOG_PREFIX}_ldsc.log"
LDSC_S_FILENAME = "ldsc_S.csv"
LDSC_V_FILENAME = "ldsc_V.csv"
LDSC_S_STAND_FILENAME = "ldsc_S_Stand.csv"
LAVAAN_MODEL_FILENAME = "lavaan_model.txt"
MODEL_FIT_FILENAME = "model_fit.csv"
MODEL_RESULTS_FILENAME = "model_results.csv"

GSEstmationType = Literal["DWLS", "ML"]


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
    Configuration for the GenomicSEM workflow.
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
class GenomicSEMTask(Task):
    """
    Run the basic GenomicSEM workflow:

    1. Munge the supplied summary statistics into the format expected by LDSC.
    2. Run multivariable LDSC on the munged sumstats to obtain a genetic
       covariance matrix.
    3. Fit a factor model expressed in the lavaan language to the LDSC output.

    Inputs:
      - ld_ref_task: produces a directory containing reference LD scores in standard form.
      - hapmap_snps_task: produces a file listing the HapMap3 SNPs to keep during munging.
      - sources: list of tasks producing tabular GWAS summary statistics in
        gwaslab format.
      - factor_model: a factor model description in lavaan syntax.

    Output: a directory asset containing munged sumstats, LDSC outputs, and
    factor model results plus logs from intermediate steps.
    """

    meta: Meta
    sources: Sequence[GenomicSEMSumstatsSource]
    ld_ref_task: Task
    hapmap_snps_task: Task
    factor_model: str
    config: GenomicSEMConfig = GenomicSEMConfig()

    def __attrs_post_init__(self):
        assert len(self.sources) >= 2, "GenomicSEM requires at least two traits"
        aliases = [s.alias for s in self.sources]
        assert len(set(aliases)) == len(aliases), (
            f"Trait aliases must be unique. Got: {aliases}"
        )

    @property
    def deps(self) -> list["Task"]:
        result: list[Task] = [self.ld_ref_task, self.hapmap_snps_task]
        for source in self.sources:
            result.append(source.task)
        return result

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        gsem = importr("GenomicSEM")

        ld_asset = fetch(self.ld_ref_task.asset_id)
        assert isinstance(ld_asset, DirectoryAsset)
        # Absolute path so values remain correct after R chdirs during munge.
        ld_path = str(ld_asset.path.resolve())

        hapmap_asset = fetch(self.hapmap_snps_task.asset_id)
        assert isinstance(hapmap_asset, FileAsset)
        hapmap_path = str(hapmap_asset.path.resolve())

        munged_dir = scratch_dir / MUNGED_SUBDIR
        munged_dir.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory() as tmp_dir_str:
            tmp_dir = Path(tmp_dir_str)
            input_files = []
            trait_names = []
            sample_sizes: list[float] = []
            sample_prevs: list[float] = []
            population_prevs: list[float] = []
            for source in self.sources:
                input_path = _write_munge_input(
                    source=source, fetch=fetch, tmp_dir=tmp_dir
                )
                input_files.append(str(input_path))
                trait_names.append(source.alias)
                sample_sizes.append(_get_sample_size(source))
                samp_prev, pop_prev = _get_prevs(source.sample_info)
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
                info_filter=self.config.info_filter,
                maf_filter=self.config.maf_filter,
            )

            munged_paths = [
                str(munged_dir / f"{name}.sumstats.gz") for name in trait_names
            ]
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
                ld_file_basename_prefix=self.config.ld_file_filename_pattern,
                ldsc_log_prefix=str(scratch_dir / LDSC_LOG_PREFIX),
            )

            _save_ldsc_outputs(covstruc=covstruc, scratch_dir=scratch_dir)

            (scratch_dir / LAVAAN_MODEL_FILENAME).write_text(self.factor_model)

            logger.debug("Running GenomicSEM::usermodel")
            model_result = gsem.usermodel(
                covstruc=covstruc,
                estimation=self.config.estimation,
                model=self.factor_model,
                CFIcalc=self.config.cfi_calc,
                std_lv=self.config.std_lv,
                fix_resid=self.config.fix_resid,
            )
            _save_model_outputs(model_result=model_result, scratch_dir=scratch_dir)

        gc.collect()
        ro.r("gc()")
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        sources: Sequence[GenomicSEMSumstatsSource],
        ld_ref_task: Task,
        hapmap_snps_task: Task,
        factor_model: str,
        config: GenomicSEMConfig = GenomicSEMConfig(),
    ) -> "GenomicSEMTask":
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
            factor_model=factor_model,
            config=config,
        )


def _get_sample_size(source: GenomicSEMSumstatsSource) -> float:
    info = source.sample_info
    if isinstance(info, BinaryPhenotypeSampleInfo):
        if info.total_sample_size is not None:
            return float(info.total_sample_size)
        return float("nan")
    if isinstance(info, QuantPhenotype):
        if info.total_sample_size is not None:
            return float(info.total_sample_size)
        return float("nan")
    raise ValueError(f"Unknown sample info type {type(info)}")


def _get_prevs(info: PhenotypeInfo) -> tuple[float, float]:
    if isinstance(info, BinaryPhenotypeSampleInfo):
        return info.sample_prevalence, info.estimated_population_prevalence
    if isinstance(info, QuantPhenotype):
        return float("nan"), float("nan")
    raise ValueError(f"Unknown sample info type {type(info)}")


def _write_munge_input(
    source: GenomicSEMSumstatsSource, fetch: Fetch, tmp_dir: Path
) -> Path:
    """
    Read the source sumstats and write a TSV with the column names munge expects.
    """
    asset = fetch(source.asset_id)
    df = (
        source.pipe.process(scan_dataframe_asset(asset, source.task.meta))
        .collect()
        .to_polars()
    )
    df = _add_sample_size_if_missing(df, sample_info=source.sample_info)

    required = [
        GWASLAB_RSID_COL,
        GWASLAB_EFFECT_ALLELE_COL,
        GWASLAB_NON_EFFECT_ALLELE_COL,
        GWASLAB_BETA_COL,
        GWASLAB_SE_COL,
        GWASLAB_P_COL,
        GWASLAB_SAMPLE_SIZE_COLUMN,
    ]
    missing = [c for c in required if c not in df.columns]
    assert not missing, f"Source '{source.alias}' is missing required columns {missing}"

    select_exprs: list[pl.Expr] = [
        pl.col(GWASLAB_RSID_COL).alias(MUNGE_SNP_COL),
        pl.col(GWASLAB_EFFECT_ALLELE_COL).alias(MUNGE_A1_COL),
        pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).alias(MUNGE_A2_COL),
        pl.col(GWASLAB_BETA_COL).alias(MUNGE_EFFECT_COL),
        pl.col(GWASLAB_SE_COL).alias(MUNGE_SE_COL),
        pl.col(GWASLAB_P_COL).alias(MUNGE_P_COL),
        pl.col(GWASLAB_SAMPLE_SIZE_COLUMN).alias(MUNGE_N_COL),
    ]
    if GWASLAB_EFFECT_ALLELE_FREQ_COL in df.columns:
        select_exprs.append(pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL).alias(MUNGE_MAF_COL))
    munge_df = df.select(select_exprs)

    output_path = tmp_dir / f"{source.alias}.sumstats.txt"
    munge_df.write_csv(output_path, separator="\t")
    logger.debug(
        f"Wrote munge input for '{source.alias}' to {output_path} "
        f"({len(munge_df)} variants)"
    )
    return output_path


def _add_sample_size_if_missing(
    df: pl.DataFrame, sample_info: PhenotypeInfo
) -> pl.DataFrame:
    if GWASLAB_SAMPLE_SIZE_COLUMN in df.columns:
        return df
    if isinstance(sample_info, BinaryPhenotypeSampleInfo):
        if sample_info.total_sample_size is None:
            raise ValueError(
                "BinaryPhenotypeSampleInfo.total_sample_size must be set when "
                "the source dataframe does not include an N column"
            )
        n = sample_info.total_sample_size
    elif isinstance(sample_info, QuantPhenotype):
        if sample_info.total_sample_size is None:
            raise ValueError(
                "QuantPhenotype.total_sample_size must be set when the source "
                "dataframe does not include an N column"
            )
        n = sample_info.total_sample_size
    else:
        raise ValueError(f"Unknown sample info type {type(sample_info)}")
    return df.with_columns(pl.lit(n).alias(GWASLAB_SAMPLE_SIZE_COLUMN))


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
    cwd_before = str(base.getwd()[0])  # ty: ignore[not-subscriptable]
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


@contextmanager
def _ld_dir_with_genomic_sem_naming(
    ld_path: str, basename_prefix: str
) -> Iterator[str]:
    """
    GenomicSEM::ldsc constructs LD score paths as "<ld>/<chr>.l2.ldscore.gz",
    so when the on-disk files have a basename prefix (e.g. "LDscore.") we
    build a tempdir of symlinks with the prefix stripped and yield that.
    """
    if not basename_prefix:
        yield ld_path
        return
    src_dir = Path(ld_path)
    with tempfile.TemporaryDirectory() as tmp_str:
        tmp = Path(tmp_str)
        for src in src_dir.iterdir():
            if not src.name.startswith(basename_prefix):
                continue
            stripped = src.name[len(basename_prefix) :]
            (tmp / stripped).symlink_to(src.resolve())
        yield str(tmp)


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

    GSEM: The Genomic SEM R module

    According to the genomic SEM wiki, he LDSC function returns an object with the following attributes:


        LDSCoutput$S is the covariance matrix (on the liability scale for case/control designs).
        LDSCoutput$V which is the sampling covariance matrix in the format expected by lavaan.
        LDSCoutput$I is the matrix of LDSC intercepts and cross-trait (i.e., bivariate) intercepts.
        LDSCoutput$N contains the sample sizes (N) for the heritabilities and sqrt(N1N2) for the co-heritabilities. These are the sample sizes provided in the munging process.
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
