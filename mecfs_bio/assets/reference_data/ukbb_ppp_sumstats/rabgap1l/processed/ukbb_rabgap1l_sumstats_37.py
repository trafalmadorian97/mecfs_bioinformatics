from mecfs_bio.assets.reference_data.ukbb_ppp_sumstats.rabgap1l.processed.stack_ukbb_rabgap1l import (
    STACK_UKBBPPP_RABGAP1L,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
)

UKBBPPP_RABGAP1l_SUMSTATS_37 = GWASLabCreateSumstatsTask(
    df_source_task=STACK_UKBBPPP_RABGAP1L,
    asset_id=AssetId("sumstats_37_ukbb_rabgap1l_gwaslab"),
    basic_check=True,
    genome_build="infer",
    liftover_to="19",
    fmt="regenie",
    harmonize_options=None,
)
