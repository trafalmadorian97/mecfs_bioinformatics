"""
Task generator to analyze decodeME GWAS 1 using the human brian atlas reference dataset
"""

from mecfs_bio.asset_generator.hba_magma_asset_generator import (
    generate_human_brain_atlas_magma_tasks,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS,
)
from mecfs_bio.build_system.task.magma.magma_plot_brain_atlas_result_with_stepwise_labels import (
    HBAIndepPlotOptions,
)
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import PlotSettings
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL

DECODE_ME_HBA_MAGMA_TASKS = generate_human_brain_atlas_magma_tasks(
    base_name="decode_me_hba_magma_tasks",
    gwas_parquet_with_rsids_task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS.join_task,
    sample_size=275488,  # this is the total sample size, which the MAGMA manual seems to imply is correct.  Could also try effective sample size.
    plot_settings=PlotSettings("plotly_white"),
    include_independent_cluster_plot=True,
    pipes=[RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL)],
    hba_indep_plot_options=HBAIndepPlotOptions(annotation_text_size=13),
)
