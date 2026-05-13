from mecfs_bio.assets.gwas.me_cfs.decode_me.auxiliary.prevalance_info import (
    DECODE_ME_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_with_annovar_37_rsids_sumstats import (
    DECODEME_ME_SUMSTATS_37_WITH_ANNOVAR_RSID,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
)
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
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

DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC = SNPHeritabilityByLDSCTask.create(
    asset_id="decode_me_gwas_1_heritability_by_ldsc",
    # pipe=pre_sldsc_pipe,
    phenotype_info=DECODE_ME_PREVALENCE_INFO,
    build="19",
    set_sample_size=275488,
    source_sumstats_task=DECODEME_ME_SUMSTATS_37_WITH_ANNOVAR_RSID,
    ld_ref_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
    pipe=IdentityPipe(),
)

DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC_MD = (
    ConvertDataFrameToMarkdownTask.create_from_result_table_task(
        asset_id="decode_me_gwas_1_ldsc_heritability_markdown",
        source_task=DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC,
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
