from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_annovar_37_sumstats import (
    DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_region_plot_task_arbitrary_locus import (
    GwasLabRegionPlotTargetLocusTask,
)

DECODE_ME_BTN1A1_REGION_PLOT_37 = GwasLabRegionPlotTargetLocusTask.create(
    asset_id="decode_me_region_plot_btn1a1_locus",
    sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
    vcf_name_for_ld="1kg_eur_hg19",
    chrom=6,
    pos=26465384,
)
