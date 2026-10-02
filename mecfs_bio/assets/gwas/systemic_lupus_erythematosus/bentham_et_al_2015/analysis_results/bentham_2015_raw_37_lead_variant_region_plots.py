from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.analysis_results.bentham_2015_raw_37_lead_variants import (
    BENTHAM_2015_RAW_LEAD_VARIANTS,
)
from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.processed_gwas_data.bentham_2015_raw_sumstats_minimal_processing import (
    BENTHAM_2015_SUMSTATS_MINIMAL_PROCESSING_FROM_RAW_BUILD_37,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_region_plots_task import (
    GwasLabRegionPlotsFromLeadVariantsTask,
)

BENTHAM_2105_RAW_LEAD_VARIANT_REGION_PLOTS = GwasLabRegionPlotsFromLeadVariantsTask(
    lead_variants_task=BENTHAM_2015_RAW_LEAD_VARIANTS,
    sumstats_task=BENTHAM_2015_SUMSTATS_MINIMAL_PROCESSING_FROM_RAW_BUILD_37,
    vcf_name_for_ld=None,
    short_id="bentham_2015_raw_37_lead_variant_region_plots",
    plot_top=5,
)
