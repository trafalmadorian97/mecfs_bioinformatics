from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap3_ppp_membership_list import (
    HAPMAP_3_MEMBERSHIP_LIST,
)
from mecfs_bio.assets.reference_data.ukbb_ppp_sumstats.rabgap1l.processed.stack_ukbb_rabgap1l import (
    STACK_UKBBPPP_RABGAP1L,
)
from mecfs_bio.build_system.task.ppp_database.construct_ppp_variant_index_task import (
    ConstructPppVariantIndexTask,
)

HAPMAP_3_PPP_DATABASE_INDEX = ConstructPppVariantIndexTask.create(
    template_protein_task=STACK_UKBBPPP_RABGAP1L,
    membership_task=HAPMAP_3_MEMBERSHIP_LIST,
    asset_id="hapmap3_ppp_database_index",
)
