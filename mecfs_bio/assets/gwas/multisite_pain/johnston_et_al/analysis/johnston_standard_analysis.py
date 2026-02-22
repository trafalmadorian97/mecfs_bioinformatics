"""
Apply standard MAGMA and S-LDSC analysis to Johnston et al.'s GWAS of multisite pain
"""

from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.processed.extracted_johnston_data import (
    JOHNSTON_ET_AL_PAIN_EXTRACTED,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import PlotSettings

JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS = (
    concrete_standard_analysis_generator_assume_already_has_rsid(
        base_name="johnston_et_al_pain",
        raw_gwas_data_task=JOHNSTON_ET_AL_PAIN_EXTRACTED,
        sample_size=387649,  # from gwas catalog
        include_hba_magma_tasks=True,
        fmt=GWASLabColumnSpecifiers(
            rsid="SNP",
            chrom="CHR",
            pos="BP",
            ea="ALLELE1",
            nea="ALLELE0",
            beta="BETA",
            se="SE",
            p="P_BOLT_LMM_INF",
            eaf="A1FREQ",
        ),
        hba_plot_settings=PlotSettings("plotly_white"),
        include_independent_cluster_plot_in_hba=True,
    )
)
