"""
Task to apply two sample mendelian randomization to GWAS data, together with associated axillary functions.
"""

import typing
from pathlib import Path, PurePath

import pandas as pd
import rpy2.robjects as ro  # type: ignore
import structlog
from attrs import frozen
from rpy2.robjects import pandas2ri  # type: ignore
from rpy2.robjects.conversion import localconverter  # type: ignore
from rpy2.robjects.packages import (  # type: ignore
    InstalledPackage,
    InstalledSTPackage,
    importr,
)

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_N_CASE_COL,
    GWASLAB_N_CONTROL_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_P_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SE_COL,
)

logger = structlog.get_logger()

TSM_RSID_COL = "SNP"
TSM_BETA_COL = "beta"
TSM_SE_COL = "se"
TSM_EFFECT_ALLELE_COL = "effect_allele"
TSM_OTHER_ALLELE_COL = "other_allele"
TSM_EAF_COL = "eaf"
TSM_CHR_COL = "chr"
TSM_POS_COL = "position"
TSM_PHENOTYPE = "Phenotype"
TSM_P_VALUE = "pval"
TSM_N_CASE = "ncase"
TSM_N_CONTROL = "ncontrol"
TSM_UNITS_COL = "units"
TSM_SAMPLE_SIZE_COL = "samplesize"


from rpy2.robjects.vectors import DataFrame as RDataFrame  # type: ignore

IgnoreOrRaise = typing.Literal["ignore", "raise"]


@frozen
class MRInputColSpec:
    rsid_col: str
    beta_col: str
    se_col: str
    ea_col: str
    nea_col: str | None = None
    eaf_col: str | None = None
    phenotype_col: str | None = None
    chrom_col: str | None = None
    pos_col: str | None = None
    pval_col: str | None = None
    n_case_col: str | None = None
    n_control_col: str | None = None
    sample_size_col: str | None = None

    errors: IgnoreOrRaise = "raise"

    def make_renamer(self) -> dict[str, str]:
        base = {
            self.rsid_col: TSM_RSID_COL,
            self.beta_col: TSM_BETA_COL,
            self.se_col: TSM_SE_COL,
            self.ea_col: TSM_EFFECT_ALLELE_COL,
        }
        if self.nea_col is not None:
            base[self.nea_col] = TSM_OTHER_ALLELE_COL
        if self.eaf_col is not None:
            base[self.eaf_col] = TSM_EAF_COL
        if self.chrom_col is not None:
            base[self.chrom_col] = TSM_CHR_COL
        if self.pos_col is not None:
            base[self.pos_col] = TSM_POS_COL
        if self.phenotype_col is not None:
            base[self.phenotype_col] = TSM_PHENOTYPE
        if self.pval_col is not None:
            base[self.pval_col] = TSM_P_VALUE
        if self.n_case_col is not None:
            base[self.n_case_col] = TSM_N_CASE
        if self.n_control_col is not None:
            base[self.n_control_col] = TSM_N_CONTROL
        if self.sample_size_col is not None:
            base[self.sample_size_col] = TSM_SAMPLE_SIZE_COL
        return base


GWASLAB_MR_INPUT_COL_SPEC = MRInputColSpec(
    rsid_col=GWASLAB_RSID_COL,
    beta_col=GWASLAB_BETA_COL,
    se_col=GWASLAB_SE_COL,
    ea_col=GWASLAB_EFFECT_ALLELE_COL,
    nea_col=GWASLAB_NON_EFFECT_ALLELE_COL,
    eaf_col=GWASLAB_EFFECT_ALLELE_FREQ_COL,
    phenotype_col=None,
    pos_col=GWASLAB_POS_COL,
    chrom_col=GWASLAB_CHROM_COL,
    n_case_col=GWASLAB_N_CASE_COL,
    n_control_col=GWASLAB_N_CONTROL_COL,
    pval_col=GWASLAB_P_COL,
    errors="ignore",
)

SUN_ET_AL_MR_INPUT_COL_SPEC_hg37 = MRInputColSpec(
    rsid_col="rsID",
    beta_col="BETA",
    se_col="SE",
    ea_col="A1",
    nea_col="A0",
    eaf_col="A1FREQ",
    phenotype_col="protein_exposure_id",
    chrom_col="CHROM_hg37",
    pos_col="GENPOS_hg37",
    pval_col=GWASLAB_P_COL,
)


RPackageType = typing.Union[InstalledSTPackage, InstalledPackage]


TSM_OUTPUT_METHOD_COL = "method"
TSM_OUTPUT_NSNP_COL = "nsnp"
TSM_OUTPUT_B_COL = "b"
TSM_OUTPUT_SE_COL = "se"
TSM_OUTPUT_P_COL = "pval"
TSM_OUTPUT_EXPOSURE_COL = "exposure"

TSM_OUTPUT_STEIGER_DIR_COL = "steiger_dir"
TSM_OUTPUT_STEIGER_P_COL = "steiger_pval"


MAIN_RESULT_DF_PATH = "mr_result.csv"

REPORT_SUBDIR_PATH = PurePath("reports")
STEIGER_RESULT_PATH = PurePath("steiger_result.csv")


@frozen
class ClumpOptions:
    pass


@frozen
class MRReportOptions:
    pass


@frozen
class SteigerFilteringOptions:
    drop_failures: bool
    p_value_thresh: float | None = None


@frozen
class TwoSampleMRConfig:
    clump_exposure_data: ClumpOptions | None
    report_options: MRReportOptions | None = None
    pre_filter_outcome_variants: bool = False
    steiger_filter: SteigerFilteringOptions | None = None


# TODO: Finish this two-sample-MR-task
# TODO: - Add A directionality test.
# Add F test option.
# Add cochranes Q test
# Add report output


NEEDED_COLS = [TSM_RSID_COL, TSM_BETA_COL, TSM_SE_COL, TSM_EFFECT_ALLELE_COL]


@frozen
class TwoSampleMRTask(Task):
    """
    Task to run mendelian randomization using the R package TwoSampleMR.
    This R package is accessed through Python via rpy2.

    Note that some of the calls to the TSMR library below (like clumping) require access
    to the OpenGWAS database.  This in turn requires an access token.
    You can get a token here: https://api.opengwas.io/
    Add to your .Renviron file the following line:
    OPENGWAS_JWT=<token>

    """

    _meta: Meta
    outcome_data_task: Task
    exposure_data_task: Task
    config: TwoSampleMRConfig
    exposure_col_spec: MRInputColSpec | None = None
    outcome_col_spec: MRInputColSpec | None = None
    exposure_pipe: DataProcessingPipe = IdentityPipe()
    outcome_pipe: DataProcessingPipe = IdentityPipe()
    mr_method_list: list[str] | None = None

    @property
    def exposure_id(self) -> AssetId:
        return self.exposure_data_task.asset_id

    @property
    def exposure_meta(self) -> Meta:
        return self.exposure_data_task.meta

    @property
    def outcome_id(self) -> AssetId:
        return self.outcome_data_task.asset_id

    @property
    def outcome_meta(self) -> Meta:
        return self.outcome_data_task.meta

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.exposure_data_task, self.outcome_data_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        tsmr = importr("TwoSampleMR")
        conv = ro.default_converter + pandas2ri.converter

        exposure_asset = fetch(self.exposure_id)
        outcome_asset = fetch(self.outcome_id)

        logger.debug("loading exposure and outcome dataframes...")
        exposure_df = (
            self.exposure_pipe.process(
                scan_dataframe_asset(exposure_asset, meta=self.exposure_meta)
            )
            .collect()
            .to_pandas()
        )
        outcome_df = (
            self.outcome_pipe.process(
                scan_dataframe_asset(outcome_asset, meta=self.outcome_meta)
            )
            .collect()
            .to_pandas()
        )

        logger.debug("Renaming columns...")

        if self.exposure_col_spec is not None:
            exposure_df = exposure_df.rename(
                columns=self.exposure_col_spec.make_renamer(),
                errors=self.exposure_col_spec.errors,
            )

        if self.outcome_col_spec is not None:
            outcome_df = outcome_df.rename(
                columns=self.outcome_col_spec.make_renamer(),
                errors=self.outcome_col_spec.errors,
            )

        assert set(NEEDED_COLS) <= set(exposure_df.columns)
        assert set(NEEDED_COLS) <= set(outcome_df.columns)

        outcome_df = pre_filter_outcome_variants(
            exposure_df=exposure_df, outcome_df=outcome_df, config=self.config
        )

        exposure_rdf, outcome_rdf = convert_outcome_and_exposure_to_r(
            exposure_df=exposure_df, outcome_df=outcome_df
        )

        f_exposure_rdf, f_outcome_rdf = format_data_no_conversion(
            exposure_rdf, outcome_rdf, tsmr=tsmr
        )

        f_exposure_rdf = optionally_clump_exposure_data_no_conversion(
            formatted_exposure=f_exposure_rdf,
            clump_options=self.config.clump_exposure_data,
            tsmr=tsmr,
        )

        harmonized = harmonize_data_no_conversion(
            formatted_exposure=f_exposure_rdf,
            formatted_outcome=f_outcome_rdf,
            tsmr=tsmr,
        )

        harmonized = steiger_filtering_write_output(
            harmonized=harmonized,
            options=self.config.steiger_filter,
            scratch_dir=scratch_dir,
            tsmr=tsmr,
        )

        # with localconverter(conv):
        #     harm_py = ro.conversion.get_conversion().rpy2py(harmonized)
        gen_mr_report(
            harmonized=harmonized,
            options=self.config.report_options,
            target_dir=scratch_dir / REPORT_SUBDIR_PATH,
            tsmr=tsmr,
        )

        result_rdf = run_tsmr_on_harmonized_data_no_conversion(
            harmonized=harmonized, tsmr=tsmr, method_list=self.mr_method_list
        )
        logger.info("Converting R dataframe to pandas dataframe...")
        with localconverter(conv):
            result_df = ro.conversion.get_conversion().rpy2py(result_rdf)

        result_df.to_csv(scratch_dir / MAIN_RESULT_DF_PATH, index=False)
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        outcome_data_task: Task,
        exposure_data_task: Task,
        config: TwoSampleMRConfig,
        exposure_pipe: DataProcessingPipe,
        outcome_pipe: DataProcessingPipe,
        exposure_col_spec: MRInputColSpec,
        outcome_col_spec: MRInputColSpec,
        method_list: list[str] | None = None,
    ):
        outcome_meta = outcome_data_task.meta
        assert isinstance(outcome_meta, FilteredGWASDataMeta)
        meta = ResultDirectoryMeta(
            id=AssetId(asset_id),
            trait=outcome_meta.trait,
            project=outcome_meta.project,
        )
        return cls(
            meta=meta,
            outcome_data_task=outcome_data_task,
            exposure_data_task=exposure_data_task,
            config=config,
            exposure_pipe=exposure_pipe,
            outcome_pipe=outcome_pipe,
            exposure_col_spec=exposure_col_spec,
            outcome_col_spec=outcome_col_spec,
            mr_method_list=method_list,
        )


def steiger_filtering_write_output(
    harmonized: RDataFrame,
    options: SteigerFilteringOptions | None,
    scratch_dir: Path,
    tsmr: RPackageType,
) -> RDataFrame:
    if options is None:
        return harmonized
    logger.debug("Performing steiger filtering")
    filtered = tsmr.steiger_filtering(harmonized)
    conv = ro.default_converter + pandas2ri.converter
    with localconverter(conv):
        py_filtered = ro.conversion.get_conversion().rpy2py(filtered)
        assert not py_filtered[TSM_OUTPUT_STEIGER_P_COL].isnull().any()
        prerows = len(py_filtered)
        logger.debug(f"Pre steiger filtering df has {prerows} rows")
        py_filtered.to_csv(scratch_dir / STEIGER_RESULT_PATH, index=False)
        if options.drop_failures:
            py_filtered = py_filtered.loc[py_filtered[TSM_OUTPUT_STEIGER_DIR_COL]]
            if options.p_value_thresh is not None:
                py_filtered = py_filtered.loc[
                    py_filtered[TSM_OUTPUT_STEIGER_P_COL] <= options.p_value_thresh
                ]
            harmonized = ro.conversion.get_conversion().py2rpy(py_filtered)

        logger.debug(
            f"Post steiger filtering df has {len(py_filtered)} rows, after dropping {prerows - len(py_filtered)}"
        )
        return harmonized


@frozen
class TwoSampleMRResult:
    result: pd.DataFrame


def run_two_sample_mr(
    exposure_df: pd.DataFrame,
    outcome_df: pd.DataFrame,
    config: TwoSampleMRConfig,
) -> TwoSampleMRResult:
    tsmr = importr("TwoSampleMR")
    # base = importr('base')
    formatted_exposure, formatted_outcome = format_data(
        exposure_df, outcome_df, tsmr=tsmr
    )
    return run_tsmr_on_formatted_data(
        formatted_exposure, formatted_outcome, config, tsmr=tsmr
    )


def format_data_no_conversion(
    exposure_rdf: RDataFrame,
    outcome_rdf: RDataFrame,
    tsmr: RPackageType,
) -> tuple[RDataFrame, RDataFrame]:
    logger.debug("formatting exposure...")
    formatted_exposure = tsmr.format_data(
        exposure_rdf,
        type="exposure",
    )

    logger.debug("formatting outcome...")
    formatted_outcome = tsmr.format_data(outcome_rdf, type="outcome")
    return formatted_exposure, formatted_outcome


def convert_outcome_and_exposure_to_r(
    exposure_df: pd.DataFrame,
    outcome_df: pd.DataFrame,
) -> tuple[RDataFrame, RDataFrame]:
    conv = ro.default_converter + pandas2ri.converter

    with localconverter(conv):
        logger.debug("converting exposure to r...")
        exposure_rdf = ro.conversion.get_conversion().py2rpy(exposure_df)
        logger.debug("converting outcome to r...")
        outcome_rdf = ro.conversion.get_conversion().py2rpy(outcome_df)
    return exposure_rdf, outcome_rdf


def format_data(
    exposure_df: pd.DataFrame,
    outcome_df: pd.DataFrame,
    tsmr: RPackageType,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    conv = ro.default_converter + pandas2ri.converter

    with localconverter(conv):
        logger.debug("formatting exposure...")
        formatted_exposure = tsmr.format_data(
            exposure_df,
            type="exposure",
        )

        logger.debug("formatting outcome...")
        formatted_outcome = tsmr.format_data(outcome_df, type="outcome")
        return formatted_exposure, formatted_outcome


def harmonize_data(
    formatted_exposure: pd.DataFrame,
    formatted_outcome: pd.DataFrame,
    tsmr: RPackageType,
) -> pd.DataFrame:
    conv = ro.default_converter + pandas2ri.converter
    with localconverter(conv):
        logger.debug("Harmonizing outcome and exposure...")
        harmonized = tsmr.harmonise_data(formatted_exposure, formatted_outcome)
        return harmonized


def harmonize_data_no_conversion(
    formatted_exposure: RDataFrame,
    formatted_outcome: RDataFrame,
    tsmr: RPackageType,
) -> RDataFrame:
    logger.debug("Harmonizing outcome and exposure...")
    harmonized = tsmr.harmonise_data(formatted_exposure, formatted_outcome)
    return harmonized


def optionally_clump_exposure_data(
    formatted_exposure: pd.DataFrame,
    clump_options: ClumpOptions | None,
    tsmr: RPackageType,
):
    if clump_options is None:
        return formatted_exposure
    logger.debug(f"shape of exposure df before clumping:{formatted_exposure.shape}")
    conv = ro.default_converter + pandas2ri.converter
    with localconverter(conv):
        formatted_exposure = tsmr.clump_data(
            formatted_exposure,
        )
    logger.debug(f"shape of exposure df after clumping:{formatted_exposure.shape}")
    return formatted_exposure


def optionally_clump_exposure_data_no_conversion(
    formatted_exposure: RDataFrame,
    clump_options: ClumpOptions | None,
    tsmr: RPackageType,
) -> RDataFrame:
    if clump_options is None:
        return formatted_exposure
    logger.debug("clumping")
    formatted_exposure = tsmr.clump_data(
        formatted_exposure,
    )

    logger.debug(" done clumping")

    return formatted_exposure


def run_tsmr_on_formatted_data(
    formatted_exposure: pd.DataFrame,
    formatted_outcome: pd.DataFrame,
    config: TwoSampleMRConfig,
    tsmr: RPackageType,
) -> TwoSampleMRResult:
    conv = ro.default_converter + pandas2ri.converter

    with localconverter(conv):
        if config.clump_exposure_data:
            logger.debug("Clumping exposure data...")
            formatted_exposure = tsmr.clump_data(formatted_exposure)
        logger.debug("Harmonizing outcome and exposure...")
        harmonized = tsmr.harmonise_data(formatted_exposure, formatted_outcome)
        return run_tsmr_on_harmonized_data(harmonized=harmonized, tsmr=tsmr)


def run_tsmr_on_harmonized_data(
    harmonized: pd.DataFrame,
    tsmr: RPackageType,
) -> TwoSampleMRResult:
    conv = ro.default_converter + pandas2ri.converter
    with localconverter(conv):
        logger.debug("performing Mendelian randomization...")
        output = tsmr.mr(harmonized)
        return TwoSampleMRResult(output)


def run_tsmr_on_harmonized_data_no_conversion(
    harmonized: RDataFrame,
    tsmr: RPackageType,
    method_list: list[str] | None = None,
) -> RDataFrame:
    logger.debug("performing Mendelian randomization...")
    if method_list is None:
        method_dict = {}
    else:
        method_dict = {"method_list": method_list}
    output = tsmr.mr(harmonized, **method_dict)
    return output


def gen_mr_report(
    harmonized: RDataFrame,
    options: MRReportOptions | None,
    target_dir: Path,
    tsmr: RPackageType,
):
    if options is None:
        return
    tsmr.mr_report(harmonized, str(target_dir))


def pre_filter_outcome_variants(
    exposure_df: pd.DataFrame, outcome_df: pd.DataFrame, config: TwoSampleMRConfig
) -> pd.DataFrame:
    if not config.pre_filter_outcome_variants:
        return outcome_df
    logger.debug("filtering outcome variants...")
    return outcome_df.loc[outcome_df[TSM_RSID_COL].isin(exposure_df[TSM_RSID_COL])]
