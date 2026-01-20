from pathlib import PurePath

import attrs

from mecfs_bio.assets.reference_data.pqtls.processed.sun_et_al_2023_pqtls_combined_extracted import \
    SUN_ET_AL_2023_COMBINED_PQTLS_EXTRACTED
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec, DataFrameTextFormat
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.copy_file_from_directory_task import CopyFileFromDirectoryTask
from mecfs_bio.build_system.task.fdr_multiple_testing_table_task import MultipleTestingTableTask
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaIfNeededPipe
from mecfs_bio.build_system.task.pipes.compute_p_pipe import ComputePIfNeededPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEIfNeededPipe
from mecfs_bio.build_system.task.pipes.concat_str_pipe import ConcatStrPipe
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.filter_rows_by_value import FilterRowsByValue
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.build_system.task.pipes.str_split_exact_col import SplitExactColPipe
from mecfs_bio.build_system.task.plot_mr_effect_measure_task import PlotMREffectMeasure, EffectMeasurePlotConfig
from mecfs_bio.build_system.task.two_sample_mr_task import TwoSampleMRTask, TwoSampleMRConfig, SteigerFilteringOptions, \
    SUN_ET_AL_MR_INPUT_COL_SPEC_hg37, GWASLAB_MR_INPUT_COL_SPEC, TSM_UNITS_COL, MAIN_RESULT_DF_PATH, TSM_OUTPUT_P_COL, \
    TSM_OUTPUT_B_COL, TSM_OUTPUT_SE_COL, TSM_OUTPUT_EXPOSURE_COL
from mecfs_bio.constants.gwaslab_constants import GWASLAB_SAMPLE_SIZE_COLUMN, GWASLAB_N_CASE_COL, GWASLAB_N_CONTROL_COL
from mecfs_bio.constants.sun_et_al_pqtl_constants import SUN_ASSAY_TARGET, SUN_TARGET_UNIPROT, SUN_CIS_TRANS_COL, \
    SUN_CIS

PROTEIN_EXPOSURE_COL = "protein_exposure_id"


@attrs.frozen
class BinaryOutcomeConfig:
    n_case: int
    n_control: int


@attrs.frozen
class QuantOutcomeConfig:
    sample_size: int

@attrs.frozen
class CisMRTasks:
    mr_task:Task
    extract_task: Task
    multiple_testing_task:Task
    plot_task: Task

def mr_cis_pqtl_asset_generator(
        gwas_dataframe_with_rsids: Task,
        base_name:str,
        outcome_config: BinaryOutcomeConfig| QuantOutcomeConfig
):
    outcome_pipes:list[DataProcessingPipe] =[
        ComputeBetaIfNeededPipe(),
        ComputeSEIfNeededPipe(),
        ComputePIfNeededPipe(),
    ]
    if isinstance(outcome_config, BinaryOutcomeConfig):
        outcome_pipes.extend(
            [
                SetColToConstantPipe(col_name=GWASLAB_N_CASE_COL, constant=outcome_config.n_case),
                SetColToConstantPipe(
                    col_name=GWASLAB_N_CONTROL_COL, constant=outcome_config.n_control
                ),
                SetColToConstantPipe(TSM_UNITS_COL, constant="log odds"),

            ]
        )
    elif isinstance(outcome_config, QuantOutcomeConfig):
        outcome_pipes.append(
            SetColToConstantPipe(
                GWASLAB_SAMPLE_SIZE_COLUMN, constant=outcome_config.sample_size
            )
        )

    two_sample_task = TwoSampleMRTask.create(
        asset_id=base_name+"_two_sample_mr_cis_pqtl",
        outcome_data_task=gwas_dataframe_with_rsids,
        exposure_data_task=SUN_ET_AL_2023_COMBINED_PQTLS_EXTRACTED,
        config=TwoSampleMRConfig(
            clump_exposure_data=None,
            pre_filter_outcome_variants=True,
            steiger_filter=SteigerFilteringOptions(True),
        ),
        exposure_pipe=CompositePipe(
            [
                ConcatStrPipe(
                    [SUN_ASSAY_TARGET, SUN_TARGET_UNIPROT],
                    sep="__",
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
                ComputePIfNeededPipe(),
                SetColToConstantPipe(
                    GWASLAB_SAMPLE_SIZE_COLUMN, constant=54_219
                ),  # This is the UKBB PPP sample size. source: sun et al. abstract
            ]
        ),
        exposure_col_spec=attrs.evolve(
            SUN_ET_AL_MR_INPUT_COL_SPEC_hg37, sample_size_col=GWASLAB_SAMPLE_SIZE_COLUMN
        ),
        outcome_col_spec=GWASLAB_MR_INPUT_COL_SPEC,
        outcome_pipe=CompositePipe(
            outcome_pipes
        ),
    )
    extract_mr_result_task=CopyFileFromDirectoryTask.create_result_table(
        asset_id=base_name+"_two_sample_mr_cis_pqtl_extract_result",
        source_directory_task=two_sample_task,
        path_inside_directory=PurePath(MAIN_RESULT_DF_PATH),
        extension=".csv",
        read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
    )
    multiple_testing_task=MultipleTestingTableTask.create_from_result_table_task(
        alpha=0.01,
        method="bonferroni",
        asset_id=base_name+"_two_sample_mr_cis_pqtl_multiple_testing",
        p_value_column=TSM_OUTPUT_P_COL,
        source_task=extract_mr_result_task,
        apply_filter=True,
    )

    plot_mr_result_task=PlotMREffectMeasure.create(
        asset_id=base_name+"_cis_pqtl_mr_effect_measure_plot",
        source_df_task=multiple_testing_task,
        config=EffectMeasurePlotConfig(
            y_label_col=SUN_ASSAY_TARGET,
            y_label="Protein",
            effect_size_col=TSM_OUTPUT_B_COL,
            effect_size_label="b",
            se_col=TSM_OUTPUT_SE_COL,
            ref_line_center=0,
            figsize=(14 * 1.3, 6 * 1.3),
            t_adjust=0.02,
        ),
        pre_pipe=SplitExactColPipe(
            TSM_OUTPUT_EXPOSURE_COL,
            split_by="__",
            new_col_names=(SUN_ASSAY_TARGET, SUN_TARGET_UNIPROT),
        ),
    )
    return CisMRTasks(
        mr_task=two_sample_task,
        extract_task=extract_mr_result_task,
        multiple_testing_task=multiple_testing_task,
        plot_task=plot_mr_result_task,
    )


