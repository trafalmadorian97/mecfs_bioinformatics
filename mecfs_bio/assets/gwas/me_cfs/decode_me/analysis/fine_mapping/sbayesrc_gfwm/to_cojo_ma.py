from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_p_pipe import ComputePPipe
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.build_system.task.sbayesrc.sumstats_to_cojo_ma_task import (
    SumstatsToCojoMaTask,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL

DECODE_ME_BUILD_37_ANNOVAR_RSID_COJO_MA = SumstatsToCojoMaTask.create(
    sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.join_task,
    id="decode_me_build_37_annovar_rsid_cojo_ma",
    pipe=CompositePipe(
        [RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL), ComputePPipe()]
    ),
)
