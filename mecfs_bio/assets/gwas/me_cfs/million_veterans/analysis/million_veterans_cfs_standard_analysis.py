from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import \
    concrete_standard_analysis_generator_assume_already_has_rsid, ManhattanPlotSettings
from mecfs_bio.assets.gwas.me_cfs.million_veterans.auxillary.pheotype_info import MILLION_VETERAN_CFS_PHENOTYPE_INFO
from mecfs_bio.assets.gwas.me_cfs.million_veterans.raw.million_veterans_cfs_download import MILLION_VETERANS_CFS_RAW
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import GWASLabColumnSpecifiers

MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP= concrete_standard_analysis_generator_assume_already_has_rsid(
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
    sample_size=443093, # from sumstats file
    manhattan_settings=ManhattanPlotSettings(),
    phenotype_info_for_ldsc=MILLION_VETERAN_CFS_PHENOTYPE_INFO,

)
