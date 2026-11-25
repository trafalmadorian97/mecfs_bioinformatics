from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.ankylosing_spondylitis.multi_study.fixed_effects_meta_analysis.processed.ank_spond_meta_analysis import (
    ANK_SPOND_FIXED_EFFECTS_META_ANALYSIS,
)
from mecfs_bio.build_system.sample_size_spec import PerVariantSampleSize
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_mlog10p_pipe import (
    ComputeMlog10pIfNeededPipe,
)
from mecfs_bio.build_system.task.pipes.compute_p_from_beta_se import (
    ComputePFromBetaSEPipeIfNeeded,
)

ANK_SPOND_FIXED_EFFECTS_STANDARD_ANALYSIS = concrete_standard_analysis_generator_assume_already_has_rsid(
    base_name="fixed_effects_ank_spond_meta_analysis",
    raw_gwas_data_task=ANK_SPOND_FIXED_EFFECTS_META_ANALYSIS,
    fmt="gwaslab",
    sample_size=PerVariantSampleSize(),  # per-variant sum of contributing studies' effective sample sizes
    pre_pipe=CompositePipe(
        [ComputePFromBetaSEPipeIfNeeded(), ComputeMlog10pIfNeededPipe()]
    ),
    include_hba_magma_tasks=True,
    include_independent_cluster_plot_in_hba=True,
)
