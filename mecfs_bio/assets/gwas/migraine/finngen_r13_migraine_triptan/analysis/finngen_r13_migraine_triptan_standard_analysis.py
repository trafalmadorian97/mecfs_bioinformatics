from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    ManhattanPlotSettings,
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.migraine.finngen_r13.analysis.finngen_r13_column_specifiers import (
    FINNGEN_R13_COLUMN_SPECIFIERS,
    FINNGEN_R13_RSID_COL,
)
from mecfs_bio.assets.gwas.migraine.finngen_r13_migraine_triptan.auxiliary.prevalence_info import (
    FINNGEN_R13_MIGRAINE_TRIPTAN_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.migraine.finngen_r13_migraine_triptan.processed.finngen_r13_migraine_triptan_common_snps import (
    FINNGEN_R13_MIGRAINE_TRIPTAN_COMMON_SNPS_TASK,
)
from mecfs_bio.build_system.task.pipes.drop_null_pipe import DropNullsPipe

FINNGEN_R13_MIGRAINE_TRIPTAN_STANDARD_ANALYSIS = concrete_standard_analysis_generator_assume_already_has_rsid(
    base_name="finngen_r13_migraine_triptan",
    raw_gwas_data_task=FINNGEN_R13_MIGRAINE_TRIPTAN_COMMON_SNPS_TASK,
    fmt=FINNGEN_R13_COLUMN_SPECIFIERS,
    sample_size=FINNGEN_R13_MIGRAINE_TRIPTAN_PREVALENCE_INFO.total_sample_size,
    pre_pipe=DropNullsPipe(subset=[FINNGEN_R13_RSID_COL]),
    phenotype_info_for_ldsc=FINNGEN_R13_MIGRAINE_TRIPTAN_PREVALENCE_INFO,
    manhattan_settings=ManhattanPlotSettings(),
    include_h_magma_tasks=True,
)
