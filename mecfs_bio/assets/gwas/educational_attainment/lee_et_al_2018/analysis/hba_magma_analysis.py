"""
Task generator to use the Human Brain Atlas reference dataset to analyze the educational attainment GWAS via MAGMA.
"""

from mecfs_bio.asset_generator.hba_magma_asset_generator import (
    generate_human_brain_atlas_magma_tasks,
)
from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.processed_gwas_data.lee_et_al_magma_task_generator import (
    LEE_ET_AL_2018_COMBINED_MAGMA_TASKS,
)
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import PlotSettings

LEE_ET_AL_2018_HBA_MAGMA_TASKS_EDU = generate_human_brain_atlas_magma_tasks(
    base_name="lee_et_al_2018_edu",
    gwas_parquet_with_rsids_task=LEE_ET_AL_2018_COMBINED_MAGMA_TASKS.parquet_file_task,
    sample_size=257841,
    plot_settings=PlotSettings("plotly_white"),
    include_independent_cluster_plot=True,
)
