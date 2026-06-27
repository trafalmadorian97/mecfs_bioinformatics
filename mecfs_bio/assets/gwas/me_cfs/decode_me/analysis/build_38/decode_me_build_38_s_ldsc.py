"""
Generate a chain of tasks to apply S-LDSC to DecodeME GWAS 1 in build 38
"""

from mecfs_bio.asset_generator.concrete_sldsc_generator import (
    standard_sldsc_task_generator,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.build_38.decode_me_gwas_1_assign_rsids_build_38 import \
    DECODE_ME_GWAS_1_ASSIGN_RSIDS_BUILD_38
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_annovar_37_sumstats import (
    DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
)
from mecfs_bio.build_system.task.pipes.compute_p_pipe import ComputePIfNeededPipe

DECODE_ME_S_LDSC_BUILD_38 = standard_sldsc_task_generator(
    sumstats_task=DECODE_ME_GWAS_1_ASSIGN_RSIDS_BUILD_38,
    base_name="build_38_decode_me_gwas_1",
    pre_pipe=ComputePIfNeededPipe(),
)
