from mecfs_bio.asset_generator.hba_magma_asset_generator import (
    generate_human_brain_atlas_magma_tasks,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_37_harmonized_assign_rsid_via_snp150_annovar import (
    LIU_ET_AL_2023_ASSIGN_RSID_VIA_SNP150_ANNOVAR,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_constants import GWASLAB_RSID_COL
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import PlotSettings
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe

IBD_HBA_MAGMA_TASKS = generate_human_brain_atlas_magma_tasks(
    base_name="liu_et_al_ibd",
    gwas_parquet_with_rsids_task=LIU_ET_AL_2023_ASSIGN_RSID_VIA_SNP150_ANNOVAR,
    sample_size=59957,
    plot_settings=PlotSettings("plotly_white"),
    pipes=[
        RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
    ],
)
