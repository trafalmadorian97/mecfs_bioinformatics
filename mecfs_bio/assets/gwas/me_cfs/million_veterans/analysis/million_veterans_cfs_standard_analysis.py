from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    ManhattanPlotSettings,
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.me_cfs.million_veterans.auxiliary.pheotype_info import (
    MILLION_VETERAN_CFS_PHENOTYPE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.million_veterans.raw.million_veterans_cfs_download import (
    MILLION_VETERANS_CFS_RAW,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaIfNeededPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe

MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP = (
    concrete_standard_analysis_generator_assume_already_has_rsid(
        base_name="million_veterans_cfs",
        raw_gwas_data_task=MILLION_VETERANS_CFS_RAW,
        fmt=GWASLabColumnSpecifiers(
            rsid="rsid",
            chrom="chromosome",
            pos="base_pair_location",
            ea="effect_allele",
            nea="other_allele",
            OR="odds_ratio",
            eaf="effect_allele_frequency",
            p="p_value",
            n="n",
        ),
        sample_size=443093,  # from sumstats file
        manhattan_settings=ManhattanPlotSettings(),
        phenotype_info_for_ldsc=MILLION_VETERAN_CFS_PHENOTYPE_INFO,
        pre_sldsc_pipe=CompositePipe([ComputeBetaIfNeededPipe(), ComputeSEPipe()]),
    )
)
