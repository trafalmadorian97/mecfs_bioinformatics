"""
Task to use pQTLs derived from the Pharma Proteomics Project as exposures in a Mendelian randomization analysis of the Han et al. asthma GWAS.
"""

import attrs

from mecfs_bio.assets.gwas.asthma.han_et_al_2022.analysis.asthma_standard_analysis import (
    HAN_ASTHMA_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.reference_data.pqtls.processed.sun_et_al_2023_pqtls_combined_extracted import (
    SUN_ET_AL_2023_COMBINED_PQTLS_EXTRACTED,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.compute_p_pipe import (
    ComputePIfNeededPipe,
    ComputePPipe,
)
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.pipes.concat_str_pipe import ConcatStrPipe
from mecfs_bio.build_system.task.pipes.filter_rows_by_value import FilterRowsByValue
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.build_system.task.two_sample_mr_task import (
    GWASLAB_MR_INPUT_COL_SPEC,
    TSM_UNITS_COL,
    SteigerFilteringOptions,
    SUN_ET_AL_MR_INPUT_COL_SPEC_hg37,
    TwoSampleMRConfig,
    TwoSampleMRTask,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_MLOG10P_COL,
    GWASLAB_N_CASE_COL,
    GWASLAB_N_CONTROL_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
)
from mecfs_bio.constants.sun_et_al_pqtl_constants import (
    SUN_ASSAY_TARGET,
    SUN_CIS,
    SUN_CIS_TRANS_COL,
    SUN_MLOG10_P_COL,
    SUN_TARGET_UNIPROT,
)

PROTEIN_EXPOSURE_COL = "protein_exposure_id"

HAN_2022_ASTHMA_TSMR = TwoSampleMRTask.create(
    asset_id="two_sample_mr_asthma_han",
    outcome_data_task=HAN_ASTHMA_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
    exposure_data_task=SUN_ET_AL_2023_COMBINED_PQTLS_EXTRACTED,
    config=TwoSampleMRConfig(
        clump_exposure_data=None,
        pre_filter_outcome_variants=True,
        steiger_filter=SteigerFilteringOptions(True),
        # report_options=MRReportOptions(),
    ),
    exposure_pipe=CompositePipe(
        [
            ConcatStrPipe(
                [SUN_ASSAY_TARGET, SUN_TARGET_UNIPROT],
                sep="_",
                new_col_name=PROTEIN_EXPOSURE_COL,
            ),
            # FilterRowsByValue(
            #     target_column="protein_exposure_id", valid_values=["IL1R1_P14778"]
            # ),
            #
            # FilterRowsByValue(
            #     target_column=PROTEIN_EXPOSURE_COL,
            #     valid_values=["TLR1_Q15399", "IL1R1_P14778"],
            # ),
            FilterRowsByValue(target_column=SUN_CIS_TRANS_COL, valid_values=[SUN_CIS]),
            RenameColPipe(SUN_MLOG10_P_COL, GWASLAB_MLOG10P_COL),
            ComputePPipe(),
            SetColToConstantPipe(
                GWASLAB_SAMPLE_SIZE_COLUMN, constant=54_219
            ),  # source: sun et al. abstract
        ]
    ),
    exposure_col_spec=attrs.evolve(
        SUN_ET_AL_MR_INPUT_COL_SPEC_hg37, sample_size_col=GWASLAB_SAMPLE_SIZE_COLUMN
    ),
    outcome_col_spec=GWASLAB_MR_INPUT_COL_SPEC,
    outcome_pipe=CompositePipe(
        [
            ComputeBetaPipe(),
            ComputeSEPipe(),
            SetColToConstantPipe(col_name=GWASLAB_N_CASE_COL, constant=64_538),
            SetColToConstantPipe(
                col_name=GWASLAB_N_CONTROL_COL, constant=329_321
            ),  # source: https://www.nature.com/articles/s41467-020-15649-3.pdf
            ComputePIfNeededPipe(),
            SetColToConstantPipe(TSM_UNITS_COL, constant="log odds"),
        ]
    ),
)
