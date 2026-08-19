"""
The shared CSF pQTL variant index in HapMap3 mode.

Templated off one aptamer's variants intersected with the HapMap3 membership list;
its row order is the alignment slot for every per-aptamer slim file. The HapMap3
membership list is reused unchanged from the PPP database (it is a normalized
gwaslab snplist, independent of which proteomics dataset consumes it), so no CSF
copy is built.
"""

from mecfs_bio.assets.gwas.csf_pqtl.raw.csf_template_aptamer import (
    CSF_TEMPLATE_APTAMER,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap3_ppp_membership_list import (
    HAPMAP_3_MEMBERSHIP_LIST,
)
from mecfs_bio.build_system.task.csf_database.construct_csf_variant_index_task import (
    ConstructCsfVariantIndexTask,
)

HAPMAP_3_CSF_DATABASE_INDEX = ConstructCsfVariantIndexTask.create(
    template_aptamer_task=CSF_TEMPLATE_APTAMER,
    membership_task=HAPMAP_3_MEMBERSHIP_LIST,
    asset_id="hapmap3_csf_database_index",
)
