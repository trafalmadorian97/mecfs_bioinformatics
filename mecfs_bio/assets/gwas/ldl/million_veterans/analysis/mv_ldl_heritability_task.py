from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.million_veterans_ldl_eur_magma_task_generator import (
    MILLION_VETERANS_EUR_LDL_MAGMA_TASKS,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
)
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    QuantPhenotype,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_snp_heritability_by_ldsc_task import (
    SNPHeritabilityByLDSCTask,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe
from mecfs_bio.build_system.task.pipes.format_numbers_pipe import FormatFloatNumbersPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.pipes.rename_col_by_position_pipe import (
    RenameColByPositionPipe,
)
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe
from mecfs_bio.build_system.task.pipes.transpose_pipe import TransposePipe

MV_LDL_HERITABILITY_TASK = SNPHeritabilityByLDSCTask.create(
    asset_id="million_veterans_eur_ldl_ldsc_heritability",
    source_sumstats_task=MILLION_VETERANS_EUR_LDL_MAGMA_TASKS.sumstats_task,
    ld_ref_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
    phenotype_info=QuantPhenotype(),
    build="19",
    set_sample_size=404741,  # source: data file,
    pipe=IdentityPipe(),
)

MV_LDL_LDSC_RESULTS_MARKDOWN = (
    ConvertDataFrameToMarkdownTask.create_from_result_table_task(
        asset_id="million_veterans_ldl_ldsc_heritability_markdown",
        source_task=MV_LDL_HERITABILITY_TASK,
        pipe=CompositePipe(
            [
                DropColPipe(["Catagories"]),
                TransposePipe(),
                RenameColByPositionPipe(col_position=0, col_new_name="Parameter"),
                RenameColByPositionPipe(col_position=1, col_new_name="Value"),
                SelectColPipe(["Parameter", "Value"]),
                FormatFloatNumbersPipe(col="Value", format_str=".4g"),
            ]
        ),
    )
)
