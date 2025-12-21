"""
Generate a chain of tasks to apply S-LDSC to DecodeME GWAS 1
"""

from mecfs_bio.asset_generator.concrete_sldsc_generator import (
    standard_sldsc_task_generator,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_annovar_37_sumstats import (
    DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
)
from mecfs_bio.build_system.task.pipes.compute_p_pipe import ComputePIfNeededPipe

DECODE_ME_S_LDSC = standard_sldsc_task_generator(
    sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
    base_name="decode_me_gwas_1",
    pre_pipe=ComputePIfNeededPipe(),
)
