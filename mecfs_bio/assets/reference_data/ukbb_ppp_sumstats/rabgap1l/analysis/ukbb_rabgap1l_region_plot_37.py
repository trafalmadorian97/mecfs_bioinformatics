from mecfs_bio.assets.reference_data.ukbb_ppp_sumstats.rabgap1l.processed.ukbb_rabgap1l_sumstats_37 import (
    UKBBPPP_RABGAP1l_SUMSTATS_37,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_region_plot_task_arbitrary_locus import (
    GwasLabRegionPlotTargetLocusTask,
)

UKBBPPP_RABGAP1L_RABGAP1L_REGION_PLOT_37 = GwasLabRegionPlotTargetLocusTask.create(
    asset_id="rabgap1l_region_plot_rabgap1l_locus",
    sumstats_task=UKBBPPP_RABGAP1l_SUMSTATS_37,
    vcf_name_for_ld="1kg_eur_hg19",
    chrom=1,
    pos=174104743,
)
