"""
Task generator to perform standard analysis (S-LDSC and MAGMA using standard reference datasets)
on asthma GWAS from Han et al.
"""

from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.asthma.han_et_al_2022.raw.raw_asthma_data import (
    HAN_ET_AL_ASTHMA_RAW,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe

HAN_ASTHMA_STANDARD_ANALYSIS = (
    concrete_standard_analysis_generator_assume_already_has_rsid(
        base_name="han_asthma",
        raw_gwas_data_task=HAN_ET_AL_ASTHMA_RAW,
        fmt=GWASLabColumnSpecifiers(
            rsid="SNP",
            chrom="CHR",
            pos="BP",
            ea="EA",
            nea="NEA",
            eaf="EAF",
            info="INFO",
            OR="OR",
            or_95l="OR_95L",
            or_95u="OR_95U",
            p="P",
            snpid=None,
            se=None,  # need to compute
            #
        ),
        sample_size=393859,  # from summary statistics file
        include_master_gene_lists=False,
        pre_sldsc_pipe=CompositePipe([ComputeBetaPipe(), ComputeSEPipe()]),
        include_hba_magma_tasks=True,
        include_independent_cluster_plot_in_hba=True,
    )
)
