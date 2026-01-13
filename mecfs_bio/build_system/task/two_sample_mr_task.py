import typing

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

from mecfs_bio.build_system.task.gwaslab.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
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


TWO_SAMPLE_MR_RENAMER = {
    TSM_RSID_COL: GWASLAB_RSID_COL,
    TSM_BETA_COL: GWASLAB_BETA_COL,
    TSM_SE_COL: GWASLAB_SE_COL,
    TSM_EFFECT_ALLELE_COL: GWASLAB_EFFECT_ALLELE_COL,
    TSM_OTHER_ALLELE_COL: GWASLAB_NON_EFFECT_ALLELE_COL,
}
RPackageType = typing.Union[InstalledSTPackage, InstalledPackage]


TSM_OUTPUT_METHOD_COL = "method"
TSM_OUTPUT_NSNP_COL = "nsnp"
TSM_OUTPUT_B_COL = "b"
TSM_OUTPUT_SE_COL = "se"
TSM_OUTPUT_P_COL = "pval"
TSM_OUTPUT_EXPOSURE_COL = "exposure"


@frozen
class TwoSampleMRConfig:
    clump_exposure_data: bool


# TODO: Finish this two-sample-MR-task
#
#
# @frozen
# class TwoSampleMRTask(Task):
#     _meta: Meta
#     outcome_data_task: Task
#     exposure_data_task: Task
#     config: TwoSampleMRConfig
#
#     @property
#     def exposure_id(self)->AssetId :
#         return self.exposure_data_task.asset_id
#
#     @property
#     def exposure_meta(self)-> Meta:
#         return self.exposure_data_task.meta
#
#     @property
#     def outcome_id(self)->AssetId :
#         return self.outcome_data_task.asset_id
#
#     @property
#     def outcome_meta(self)-> Meta:
#         return self.outcome_data_task.meta
#
#     @property
#     def meta(self) -> Meta:
#         return self._meta
#
#
#     @property
#     def deps(self) -> list["Task"]:
#         return [self.exposure_data_task, self.outcome_data_task]
#
#
#
#
#     def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
#
#         exposure_asset = fetch(self.exposure_id)
#         outcome_asset = fetch(self.outcome_id)
#
#         exposure_df =scan_dataframe_asset(exposure_asset, meta=self.exposure_meta).collect().to_pandas()
#         outcome_df = scan_dataframe_asset(outcome_asset, meta=self.outcome_meta).collect().to_pandas()
#         exposure_df = exposure_df.rename(columns=TWO_SAMPLE_MR_RENAMER)
#         outcome_df = outcome_df.rename(columns=TWO_SAMPLE_MR_RENAMER)
#
#


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
    conv = ro.default_converter + pandas2ri.converter

    with localconverter(conv):
        formatted_exposure = tsmr.format_data(
            exposure_df,
            type="exposure",
        )
        formatted_outcome = tsmr.format_data(outcome_df, type="outcome")
        if config.clump_exposure_data:
            logger.debug("Clumping exposure data...")
            formatted_exposure = tsmr.clump_data(formatted_exposure)
        logger.debug("Harmonizing outcome and exposure...")
        harmonized = tsmr.harmonise_data(formatted_exposure, formatted_outcome)
        logger.debug("performing Mendelian randomization...")
        output = tsmr.mr(harmonized)
    return TwoSampleMRResult(output)
