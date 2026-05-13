from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.me_cfs.neale_lab.auxiliary.neale_lab_cfs_phenotype_info import (
    NEALE_LAB_CFS_PHENOTYPE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.neale_lab.processed.neale_lab_cfs_gwas_joined import (
    NEALE_LAB_CFS_JOINED,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)

NEALE_LAB_CFS_STANDARD_ANALYSIS_TASK_GROUP = (
    concrete_standard_analysis_generator_assume_already_has_rsid(
        base_name="neale_lab_cfs",
        raw_gwas_data_task=NEALE_LAB_CFS_JOINED,
        fmt=GWASLabColumnSpecifiers(
            snpid="variant",
            chrom="chr",
            pos="pos",
            ea="alt",
            nea="ref",
            beta="beta",
            rsid="rsid",
            se="se",
            p="pval",
        ),
        sample_size=359482 + 1659,  # source: Neale Lab phenotype file
        phenotype_info_for_ldsc=NEALE_LAB_CFS_PHENOTYPE_INFO,
    )
)
