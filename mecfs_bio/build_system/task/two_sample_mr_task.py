import typing
from pathlib import Path

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
from mecfs_bio.build_system.meta import asset_id
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SE_COL, GWASLAB_EFFECT_ALLELE_FREQ_COL, GWASLAB_POS_COL, GWASLAB_CHROM_COL,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.wf.base_wf import WF

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


IgnoreOrRaise = typing.Literal["ignore", "raise"]

@frozen
class MRInputColSpec:
    rsid_col: str
    beta_col: str
    se_col: str
    ea_col:str
    nea_col:str|None=None
    eaf_col:str|None=None
    phenotype_col: str|None=None
    chrom_col: str|None=None
    pos_col: str|None=None
    errors:IgnoreOrRaise = "raise"
    def make_renamer(self)-> dict[str,str]:
        base = {
            self.rsid_col:TSM_RSID_COL,
            self.beta_col:TSM_BETA_COL,
            self.se_col:TSM_SE_COL,
            self.ea_col:TSM_EFFECT_ALLELE_COL
        }
        if self.nea_col is not None:
            base[self.nea_col]=TSM_OTHER_ALLELE_COL
        if self.eaf_col is not None:
            base[self.eaf_col]= TSM_EAF_COL
        if self.chrom_col is not None:
            base[self.chrom_col]= TSM_CHR_COL
        if self.pos_col is not None:
            base[self.pos_col]= TSM_POS_COL
        if self.phenotype_col is not None:
            base[self.phenotype_col]= TSM_PHENOTYPE
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
    errors="ignore"
)

SUN_ET_AL_MR_INPUT_COL_SPEC_hg37 = MRInputColSpec(
    rsid_col="rsID",
    beta_col="BETA",
    se_col="SE",
    ea_col="A1",
    nea_col="A0",
    eaf_col="A1FREQ",

)




# TWO_SAMPLE_MR_RENAMER = {
#     TSM_RSID_COL: GWASLAB_RSID_COL,
#     TSM_BETA_COL: GWASLAB_BETA_COL,
#     TSM_SE_COL: GWASLAB_SE_COL,
#     TSM_EFFECT_ALLELE_COL: GWASLAB_EFFECT_ALLELE_COL,
#     TSM_OTHER_ALLELE_COL: GWASLAB_NON_EFFECT_ALLELE_COL,
# }
RPackageType = typing.Union[InstalledSTPackage, InstalledPackage]


TSM_OUTPUT_METHOD_COL = "method"
TSM_OUTPUT_NSNP_COL = "nsnp"
TSM_OUTPUT_B_COL = "b"
TSM_OUTPUT_SE_COL = "se"
TSM_OUTPUT_P_COL = "pval"
TSM_OUTPUT_EXPOSURE_COL = "exposure"


MAIN_RESULT_DF_PATH = "mr_result.csv"
REPORT_SUBDIR_PATH="reports"


@frozen
class MRReportOptions:
    pass

@frozen
class TwoSampleMRConfig:
    clump_exposure_data: bool
    report_options: MRReportOptions|None = None


# TODO: Finish this two-sample-MR-task
#TODO: - Add A directionality test.
# Add F test option.
# Add cochranes Q test
# Add report output




@frozen
class TwoSampleMRTask(Task):
    _meta: Meta
    outcome_data_task: Task
    exposure_data_task: Task
    config: TwoSampleMRConfig
    exposure_pipe: DataProcessingPipe
    outcome_pipe: DataProcessingPipe
    exposure_col_spec: MRInputColSpec
    outcome_col_spec: MRInputColSpec

    @property
    def exposure_id(self)->AssetId :
        return self.exposure_data_task.asset_id

    @property
    def exposure_meta(self)-> Meta:
        return self.exposure_data_task.meta

    @property
    def outcome_id(self)->AssetId :
        return self.outcome_data_task.asset_id

    @property
    def outcome_meta(self)-> Meta:
        return self.outcome_data_task.meta

    @property
    def meta(self) -> Meta:
        return self._meta


    @property
    def deps(self) -> list["Task"]:
        return [self.exposure_data_task, self.outcome_data_task]




    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        tsmr = importr("TwoSampleMR")

        exposure_asset = fetch(self.exposure_id)
        outcome_asset = fetch(self.outcome_id)

        exposure_df = self.exposure_pipe.process(scan_dataframe_asset(exposure_asset, meta=self.exposure_meta)).collect().to_pandas()
        outcome_df =  self.outcome_pipe.process(scan_dataframe_asset(outcome_asset, meta=self.outcome_meta)).collect().to_pandas()

        exposure_df = exposure_df.rename(columns=self.exposure_col_spec.make_renamer(), errors=self.exposure_col_spec.errors)
        outcome_df = outcome_df.rename(columns=self.outcome_col_spec.make_renamer(), errors=self.outcome_col_spec.errors)


        f_exposure, f_outcome =format_data(exposure_df, outcome_df,tsmr=tsmr)

        harmonized = harmonize_data(formatted_exposure=f_exposure, formatted_outcome=outcome_df,tsmr=tsmr)

        result =run_tsmr_on_harmonized_data(
            harmonized=harmonized,
            tsmr=tsmr,
        )



        result.result.to_csv(scratch_dir / MAIN_RESULT_DF_PATH)
        return DirectoryAsset(
            scratch_dir
        )





    @classmethod
    def create(cls,
               asset_id:str,
               outcome_data_task: Task,
            exposure_data_task: Task,
            config: TwoSampleMRConfig,
            exposure_pipe: DataProcessingPipe,
            outcome_pipe: DataProcessingPipe,
        ):
        outcome_meta = outcome_data_task.meta
        assert isinstance(outcome_meta, FilteredGWASDataMeta)
        meta =ResultDirectoryMeta(
            asset_id=AssetId(asset_id),
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
        )




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
    formatted_exposure, formatted_outcome = format_data(exposure_df, outcome_df, tsmr=tsmr)
    return run_tsmr_on_formatted_data(formatted_exposure, formatted_outcome, config, tsmr=tsmr)


def format_data(
        exposure_df: pd.DataFrame,
        outcome_df: pd.DataFrame,
        tsmr: RPackageType,
    )-> typing.Tuple[pd.DataFrame, pd.DataFrame]:
    conv = ro.default_converter + pandas2ri.converter

    with localconverter(conv):
        formatted_exposure = tsmr.format_data(
            exposure_df,
            type="exposure",
        )
        formatted_outcome = tsmr.format_data(outcome_df, type="outcome")
        return formatted_exposure, formatted_outcome


def harmonize_data(
        formatted_exposure: pd.DataFrame,
        formatted_outcome: pd.DataFrame,
        tsmr: RPackageType,
)->pd.DataFrame:
    conv = ro.default_converter + pandas2ri.converter
    with localconverter(conv):
        logger.debug("Harmonizing outcome and exposure...")
        harmonized = tsmr.harmonise_data(formatted_exposure, formatted_outcome)
        return harmonized




def run_tsmr_on_formatted_data(
        formatted_exposure: pd.DataFrame,
                               formatted_outcome: pd.DataFrame,

                               config: TwoSampleMRConfig,
                               tsmr: RPackageType) -> TwoSampleMRResult:
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
        )   -> TwoSampleMRResult:
    conv = ro.default_converter + pandas2ri.converter
    with localconverter(conv):
        logger.debug("performing Mendelian randomization...")
        output = tsmr.mr(harmonized)
        return TwoSampleMRResult(output)


def gen_mr_report(
harmonized: pd.DataFrame,
       options: MRReportOptions| None ,
        target_dir: Path,
        tsmr: RPackageType,
):
    if harmonized is None:
        return
    conv = ro.default_converter + pandas2ri.converter
    with localconverter(conv):
        tsmr.mr_report(harmonized, str(target_dir/REPORT_SUBDIR_PATH))

