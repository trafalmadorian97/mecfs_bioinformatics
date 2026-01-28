from mecfs_bio.assets.reference_data.ukbb_ppp_sumstats.rabgap1l.processed.stack_ukbb_rabgap1l import (
    STACK_UKBBPPP_RABGAP1L,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
    GWASLabVCFRef,
    HarmonizationOptions,
)

UKBBPPP_RABGAP1l_SUMSTATS_37_HARMONIZED = GWASLabCreateSumstatsTask(
    df_source_task=STACK_UKBBPPP_RABGAP1L,
    asset_id=AssetId("sumstats_37_ukbb_rabgap1l"),
    basic_check=True,
    genome_build="infer",
    liftover_to="19",
    fmt="regenie",
    harmonize_options=HarmonizationOptions(
        ref_infer=GWASLabVCFRef(name="1kg_eur_hg19", ref_alt_freq="AF"),
        ref_seq="ucsc_genome_hg19",
        check_ref_files=True,
        drop_missing_from_ref=True,
        cores=4,
    ),
)
