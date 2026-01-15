from mecfs_bio.assets.gwas.asthma.han_et_al_2022.analysis.asthma_standard_analysis import (
    HAN_ASTHMA_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.reference_data.pqtls.processed.sun_et_al_2023_pqtls_combined_extracted import (
    SUN_ET_AL_2023_COMBINED_PQTLS_EXTRACTED,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.pipes.concat_str_pipe import ConcatStrPipe
from mecfs_bio.build_system.task.pipes.filter_rows_by_value import FilterRowsByValue
from mecfs_bio.build_system.task.two_sample_mr_task import (
    GWASLAB_MR_INPUT_COL_SPEC,
    SUN_ET_AL_MR_INPUT_COL_SPEC_hg37,
    TwoSampleMRConfig,
    TwoSampleMRTask,
)

HAN_2022_ASTHMA_TSMR = TwoSampleMRTask.create(
    asset_id="two_sample_mr_asthma_han",
    outcome_data_task=HAN_ASTHMA_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
    exposure_data_task=SUN_ET_AL_2023_COMBINED_PQTLS_EXTRACTED,
    config=TwoSampleMRConfig(
        clump_exposure_data=None, pre_filter_outcome_variants=True
    ),
    exposure_pipe=CompositePipe(
        [
            ConcatStrPipe(
                ["Assay Target", "Target UniProt"],
                sep="_",
                new_col_name="protein_exposure_id",
            ),
            FilterRowsByValue(
                target_column="protein_exposure_id", valid_values=["IL1R1_P14778"]
            ),
        ]
    ),
    exposure_col_spec=SUN_ET_AL_MR_INPUT_COL_SPEC_hg37,
    outcome_col_spec=GWASLAB_MR_INPUT_COL_SPEC,
    outcome_pipe=CompositePipe([ComputeBetaPipe(), ComputeSEPipe()]),
)
