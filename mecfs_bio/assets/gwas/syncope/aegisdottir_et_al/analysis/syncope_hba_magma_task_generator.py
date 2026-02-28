from mecfs_bio.asset_generator.hba_magma_asset_generator import (
    generate_human_brain_atlas_magma_tasks,
)
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.processed.sumstats_liftover_dump_to_parquet import (
    AEGISDOTTIR_SUMSTATS_DUMP_TO_PARQUET,
)
from mecfs_bio.build_system.task.magma.magma_plot_brain_atlas_result_with_stepwise_labels import (
    HBAIndepPlotOptions,
)
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import PlotSettings

HBA_MAGMA_AEGISDOTTIR_SYCNOPE_GWAS = generate_human_brain_atlas_magma_tasks(
    base_name="aegisdottir_syncope",
    gwas_parquet_with_rsids_task=AEGISDOTTIR_SUMSTATS_DUMP_TO_PARQUET,
    sample_size=946_861,  # see: https://www.cardioaragon.com/wp-content/uploads/Genetic-variants-associated-with-syncope.EHeartJ.2023.pdf
    plot_settings=PlotSettings("plotly_white"),
    include_independent_cluster_plot=True,
    hba_indep_plot_options=HBAIndepPlotOptions(annotation_text_size=10),
)
