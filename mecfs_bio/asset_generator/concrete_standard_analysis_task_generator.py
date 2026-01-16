"""
Task generators that apply a collection of standard analysis techniques to GWAS summary statistics using standard reference data.
"""

from attrs import frozen

from mecfs_bio.asset_generator.annovar_37_basic_rsid_assignment import (
    RSIDAssignmentTaskGroup,
    annovar_37_basic_rsid_assignment,
)
from mecfs_bio.asset_generator.concrete_magma_asset_generator import (
    concrete_magma_assets_generate,
)
from mecfs_bio.asset_generator.concrete_sldsc_generator import (
    standard_sldsc_task_generator,
)
from mecfs_bio.asset_generator.hba_magma_asset_generator import (
    HBAMagmaTasks,
    generate_human_brain_atlas_magma_tasks,
)
from mecfs_bio.asset_generator.labeled_lead_variants_asset_generator import (
    LabelLeadVariantsTasks,
    generate_tasks_labeled_lead_variants,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.combine_gene_lists_task import SrcGeneList
from mecfs_bio.build_system.task.gwaslab.gwaslab_constants import (
    GWASLAB_SAMPLE_SIZE_COLUMN,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
    ValidGwaslabFormat,
)
from mecfs_bio.build_system.task.magma.magma_plot_brain_atlas_result_with_stepwise_labels import (
    HBAIndepPlotOptions,
)
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import PlotSettings
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.build_system.task_generator.magma_task_generator import (
    MagmaTaskGeneratorFromRaw,
)
from mecfs_bio.build_system.task_generator.master_gene_list_task_generator import (
    MasterGeneListTasks,
    generate_master_gene_list_tasks,
)
from mecfs_bio.build_system.task_generator.sldsc_task_generator import (
    SLDSCTaskGenerator,
)


@frozen
class StandardAnalysisTaskGroup:
    """
    Collection of standard analysis tasks for GWAS summary statistics.
    """

    sldsc_tasks: SLDSCTaskGenerator
    magma_tasks: MagmaTaskGeneratorFromRaw
    labeled_lead_variant_tasks: LabelLeadVariantsTasks
    master_gene_list_tasks: MasterGeneListTasks | None
    hba_magma_tasks: HBAMagmaTasks | None = None

    def get_terminal_tasks(self) -> list[Task]:
        result = (
            list(self.sldsc_tasks.get_terminal_tasks())
            + self.magma_tasks.inner.terminal_tasks()
            + self.labeled_lead_variant_tasks.terminal_tasks()
        )
        if self.master_gene_list_tasks is not None:
            result = result + self.master_gene_list_tasks.terminal_tasks()
        if self.hba_magma_tasks is not None:
            result.extend(self.hba_magma_tasks.terminal_tasks())
        return result


def concrete_standard_analysis_generator_assume_already_has_rsid(
    base_name: str,
    raw_gwas_data_task: Task,
    fmt: ValidGwaslabFormat,
    sample_size: int,
    sample_size_for_sldsc: int | None = None,
    pre_pipe: DataProcessingPipe = IdentityPipe(),
    pre_sldsc_pipe: DataProcessingPipe = IdentityPipe(),
    include_master_gene_lists: bool = True,
    include_hba_magma_tasks: bool = False,
    include_independent_cluster_plot_in_hba: bool = False,
    hba_plot_settings: PlotSettings = PlotSettings(),
    gtex_magma_number_of_bars: int = 20,
    hba_indep_plot_options: HBAIndepPlotOptions = HBAIndepPlotOptions(),
) -> StandardAnalysisTaskGroup:
    """
    Generate standard MAGMA and S-LDSC analysis tasks for given GWAS data,
    assuming that GWAS data already contains rsids
    """

    labeled_lead_variant_task_group = generate_tasks_labeled_lead_variants(
        base_name=base_name,
        raw_gwas_data_task=raw_gwas_data_task,
        fmt=fmt,
        pre_pipe=pre_pipe,
    )
    magma_tasks = concrete_magma_assets_generate(
        base_name=base_name,
        raw_gwas_data_task=raw_gwas_data_task,
        fmt=fmt,
        sample_size=sample_size,
        pre_pipe=pre_pipe,
        number_of_bars=gtex_magma_number_of_bars,
    )
    sldsc_ss = (
        sample_size_for_sldsc if sample_size_for_sldsc is not None else sample_size
    )
    pre_sldsc_pipe = CompositePipe(
        [
            SetColToConstantPipe(
                col_name=GWASLAB_SAMPLE_SIZE_COLUMN, constant=sldsc_ss
            ),
            pre_sldsc_pipe,
        ]
    )

    sldsc_tasks = standard_sldsc_task_generator(
        sumstats_task=magma_tasks.sumstats_task,
        base_name=base_name,
        pre_pipe=pre_sldsc_pipe,
    )
    if include_master_gene_lists:
        master_gene_list_tasks = generate_master_gene_list_tasks(
            base_name=base_name,
            gene_sources=[
                SrcGeneList(
                    task=magma_tasks.inner.filtered_gene_analysis_task,
                    name="MAGMA",
                    ensemble_id_column="GENE",
                ),
                SrcGeneList(
                    task=labeled_lead_variant_task_group.labeled_lead_variants_task,
                    name="GWASLAB",
                    ensemble_id_column="Gene stable ID",
                ),
            ],
        )
    else:
        master_gene_list_tasks = None
    if include_hba_magma_tasks:
        hba_magma = generate_human_brain_atlas_magma_tasks(
            base_name=base_name,
            gwas_parquet_with_rsids_task=magma_tasks.parquet_file_task,
            sample_size=sample_size,
            plot_settings=hba_plot_settings,
            include_independent_cluster_plot=include_independent_cluster_plot_in_hba,
            hba_indep_plot_options=hba_indep_plot_options,
        )
    else:
        hba_magma = None

    return StandardAnalysisTaskGroup(
        sldsc_tasks=sldsc_tasks,
        magma_tasks=magma_tasks,
        labeled_lead_variant_tasks=labeled_lead_variant_task_group,
        master_gene_list_tasks=master_gene_list_tasks,
        hba_magma_tasks=hba_magma,
    )


@frozen
class StandardAnalysisTaskGroupAddRSIDS:
    """
    Collection to tasks to assign rsids to GWAS data, then apply standard analysis techniques.
    """

    tasks: StandardAnalysisTaskGroup
    initial_sumstats_task: Task
    assign_rsids_task_group: RSIDAssignmentTaskGroup

    def terminal_tasks(self) -> list[Task]:
        return list(self.tasks.get_terminal_tasks())


def concrete_standard_analysis_generator_no_rsid(
    base_name: str,
    raw_gwas_data_task: Task,
    fmt: ValidGwaslabFormat,
    sample_size: int,
    sample_size_for_sldsc: int | None = None,
    pre_pipe: DataProcessingPipe = IdentityPipe(),
    pre_sldsc_pipe: DataProcessingPipe = IdentityPipe(),
    include_master_gene_lists: bool = True,
) -> StandardAnalysisTaskGroupAddRSIDS:
    """

    Generate standard MAGMA and S-LDSC analysis tasks for given GWAS data,
    Assume that the GWAS data does not contain rsids,and so these need to be assigned.

    """
    sumstats_37_task = GWASLabCreateSumstatsTask(
        df_source_task=raw_gwas_data_task,
        asset_id=AssetId(base_name + "_initial_sumstats_37"),
        basic_check=True,
        genome_build="infer",
        liftover_to="19",
        fmt=fmt,
        pre_pipe=pre_pipe,
    )
    rsids_assigned_task_group = annovar_37_basic_rsid_assignment(
        sumstats_task=sumstats_37_task,
        base_name=base_name,
        use_gwaslab_rsids_convention=True,
    )
    standard_tasks = concrete_standard_analysis_generator_assume_already_has_rsid(
        base_name=base_name,
        raw_gwas_data_task=rsids_assigned_task_group.join_task,
        fmt="gwaslab",
        sample_size=sample_size,
        pre_pipe=pre_pipe,
        pre_sldsc_pipe=pre_sldsc_pipe,
        sample_size_for_sldsc=sample_size_for_sldsc,
        include_master_gene_lists=include_master_gene_lists,
    )
    return StandardAnalysisTaskGroupAddRSIDS(
        tasks=standard_tasks,
        initial_sumstats_task=sumstats_37_task,
        assign_rsids_task_group=rsids_assigned_task_group,
    )
