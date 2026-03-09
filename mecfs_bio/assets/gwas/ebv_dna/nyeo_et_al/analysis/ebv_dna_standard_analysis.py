"""
Task generator to perform standard analysis on Nyeo et al's GWAS of EBV DNA.
"""

from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_no_rsid,
)
from mecfs_bio.assets.gwas.ebv_dna.nyeo_et_al.raw.ebv_dna import NYEO_EBV_DNA_SUMSTATS
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_p_pipe import ComputePPipe

EBV_DNA_STANDARD_ANALYSIS = concrete_standard_analysis_generator_no_rsid(
    base_name="nyeo_ebv_dna",
    raw_gwas_data_task=NYEO_EBV_DNA_SUMSTATS,
    fmt=GWASLabColumnSpecifiers(
        chrom="#chrom",
        pos="pos",
        nea="ref",
        ea="alt",
        mlog10p="neg_log_pvalue",
        beta="beta",
        se="stderr_beta",
        eaf="alt_allele_freq",
    ),
    sample_size=490_560,
    pre_pipe_after_rsid_assignment=CompositePipe([ComputePPipe()]),
    drop_palindromic_in_harmonized=False,
    include_hba_magma_tasks=True,
    include_independent_cluster_plot_in_hba=True,
)
