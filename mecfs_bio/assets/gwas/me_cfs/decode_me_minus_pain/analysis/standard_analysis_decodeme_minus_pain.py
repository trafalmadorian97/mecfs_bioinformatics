import narwhals.dtypes

from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    ManhattanPlotSettings,
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_minus_pain.raw.decode_me_minus_pain_extract import (
    DECODE_ME_MINUS_PAIN_EXTRACT,
)
from mecfs_bio.build_system.sample_size_spec import PerVariantSampleSize
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    QuantPhenotype,
)
from mecfs_bio.build_system.task.pipes.cast_pipe import CastPipe
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_config import (
    LDSC_BP_COL,
    LDSC_CHR_COL,
    MUNGE_A1_COL,
    MUNGE_A2_COL,
    MUNGE_MAF_COL,
    MUNGE_SNP_COL,
    SUBTRACTION_EST_COL,
    SUBTRACTION_N_EFF_COL,
    SUBTRACTION_P_COL,
    SUBTRACTION_SE_COL,
)

GENOMIC_SEM_COL_SPECIFIERS = GWASLabColumnSpecifiers(
    rsid=MUNGE_SNP_COL,
    chrom=LDSC_CHR_COL,
    pos=LDSC_BP_COL,
    # est (the per-SNP factor effect) is aligned to the reference A1 allele,
    # so A1 is the effect allele and A2 the non-effect allele.
    ea=MUNGE_A1_COL,
    nea=MUNGE_A2_COL,
    beta=SUBTRACTION_EST_COL,
    se=SUBTRACTION_SE_COL,
    p=SUBTRACTION_P_COL,
    n=SUBTRACTION_N_EFF_COL,
    # MAF here is the reference (1000G) minor-allele frequency, not an
    # effect-allele frequency, so it maps to `maf` rather than `eaf`.
    maf=MUNGE_MAF_COL,
)
DECODE_ME_MINUS_PAIN_STANDARD_ANALYSIS = concrete_standard_analysis_generator_assume_already_has_rsid(
    base_name="decode_me_minus_johnston",
    raw_gwas_data_task=DECODE_ME_MINUS_PAIN_EXTRACT,
    fmt=GENOMIC_SEM_COL_SPECIFIERS,
    sample_size=PerVariantSampleSize(),
    phenotype_info_for_ldsc=QuantPhenotype(),  # Applying GenomicSEM makes everything act like a quant phenotype
    manhattan_settings=ManhattanPlotSettings(),
    pre_pipe=CastPipe(
        target_column=SUBTRACTION_N_EFF_COL,
        type=narwhals.dtypes.Int64(),
        new_col_name=SUBTRACTION_N_EFF_COL,
    ),
)
