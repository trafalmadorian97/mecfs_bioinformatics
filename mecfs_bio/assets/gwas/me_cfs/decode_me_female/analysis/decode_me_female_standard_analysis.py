from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_no_rsid,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_female.auxiliary.decode_me_female_phenotype_info import (
    DECODE_ME_FEMALE_PHENOTYPE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_female.processed.filtered_snps_decode_me_female import (
    DECODE_ME_FEMALE_FILTER_SNPS_TASK,
)
from mecfs_bio.build_system.task.pipes.compute_p_pipe import (
    ComputePIfNeededPipe,
)

DECODE_ME_FEMALE_STANDARD_ANALYSIS = concrete_standard_analysis_generator_no_rsid(
    base_name="decode_me_female",
    raw_gwas_data_task=DECODE_ME_FEMALE_FILTER_SNPS_TASK,
    fmt="regenie",
    sample_size=12_833 + 218_949,
    pre_pipe_after_rsid_assignment=ComputePIfNeededPipe(),
    phenotype_info_for_ldsc=DECODE_ME_FEMALE_PHENOTYPE_INFO,
)
