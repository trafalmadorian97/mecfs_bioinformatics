from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_no_rsid,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_infectious_onset.auxiliary.decode_me_infectious_onset_phenotype_info import (
    DECODE_ME_INFECTIOUS_ONSET_PHENOTYPE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_infectious_onset.processed.filtered_snps_decode_me_infectious_onset import (
    DECODE_ME_INFECTIOUS_ONSET_FILTER_SNPS_TASK,
)
from mecfs_bio.build_system.task.pipes.compute_p_pipe import (
    ComputePIfNeededPipe,
)

DECODE_ME_INFECTIOUS_ONSET_STANDARD_ANALYSIS = (
    concrete_standard_analysis_generator_no_rsid(
        base_name="decode_me_infectious_onset",
        raw_gwas_data_task=DECODE_ME_INFECTIOUS_ONSET_FILTER_SNPS_TASK,
        fmt="regenie",
        sample_size=9_738 + 259_909,
        pre_pipe_after_rsid_assignment=ComputePIfNeededPipe(),
        phenotype_info_for_ldsc=DECODE_ME_INFECTIOUS_ONSET_PHENOTYPE_INFO,
    )
)
