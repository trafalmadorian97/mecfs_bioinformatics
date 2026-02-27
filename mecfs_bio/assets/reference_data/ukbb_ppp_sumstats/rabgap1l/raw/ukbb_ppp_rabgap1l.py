from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.get_file_from_synapse_task import (
    GetFileFromSynapseTask,
)

UKBB_PPP_RABGAP1L = GetFileFromSynapseTask(
    ReferenceFileMeta(
        group="ukbb_ppp",
        sub_group="sun_et_al_2023",
        sub_folder=PurePath("raw"),
        id="rabgap1l_sumstats",
        filename="RABGAP1L_Q5R372_OID20447_v1_Inflammation",
        extension=".tar",
    ),
    synid="syn51470978",
    expected_filename="RABGAP1L_Q5R372_OID20447_v1_Inflammation.tar",
)
