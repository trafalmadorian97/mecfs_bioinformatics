"""GWAS-by-subtraction of Johnston et al. multisite chronic pain from DecodeME
GWAS-1, run entirely in Python (no rpy2).

DecodeME is the composite trait (T1): the trait we decompose. Johnston pain is
the reference trait (T2), a pure indicator of the common factor F shared by the
two. The remainder factor R is therefore the ME/CFS-specific signal left after
the pain-shared component is subtracted out.

DecodeME GWAS-1 is a case/control study (logistic), while Johnston's phenotype
is the 0-7 multisite-chronic-pain count analysed with BOLT-LMM as a continuous
trait (OLS).
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.auxiliary.prevalance_info import (
    DECODE_ME_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import (
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.auxiliary.sample_info import (
    JOHNSTON_ET_AL_SAMPLE_INFO,
)
from mecfs_bio.assets.reference_data.genomic_sem_reference.genomes_1k import (
    GENOMES1K_REFERENCE_FOR_GENOMIC_SEM,
)
from mecfs_bio.assets.reference_data.genomic_sem_reference.hapmap3_snpist import (
    HAPMAP3_SNPLIST_FOR_GENOMIC_SEM,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_p_pipe import ComputePPipe
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_config import (
    GenomicSEMGWASSumstatsSource,
    GenomicSEMSumstatsSource,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_gwas_by_subtraction_full_python_task import (
    GenomicSEMGWASBySubtractionFullPythonTask,
    GWASBySubtractionFullPythonConfig,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL

DECODE_ME_MINUS_PAIN_SUBTRACTION = GenomicSEMGWASBySubtractionFullPythonTask.create(
    asset_id="decode_me_minus_pain_subtraction_full_python",
    # Composite trait (T1): DecodeME, decomposed into the pain-shared common
    # factor F and the ME/CFS-specific remainder factor R.
    composite_trait_source=GenomicSEMGWASSumstatsSource(
        source=GenomicSEMSumstatsSource(
            task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS.join_task,
            alias="DECODE_ME",
            sample_info=DECODE_ME_PREVALENCE_INFO,
            pipe=CompositePipe(
                [
                    RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
                    ComputePPipe(),
                ]
            ),
        ),
        gwas_method="logistic",
    ),
    # Reference trait (T2): Johnston multisite chronic pain, a pure indicator
    # of the common factor F. The 0-7 pain-site count was analysed with
    # BOLT-LMM as a continuous trait, so it is an OLS source.
    reference_trait_source=GenomicSEMGWASSumstatsSource(
        source=GenomicSEMSumstatsSource(
            task=JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
            alias="Multsite_Pain",
            sample_info=JOHNSTON_ET_AL_SAMPLE_INFO,
        ),
        gwas_method="ols",
    ),
    ld_ref_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
    hapmap_snps_task=HAPMAP3_SNPLIST_FOR_GENOMIC_SEM,
    sumstats_ref_task=GENOMES1K_REFERENCE_FOR_GENOMIC_SEM,
    # The extracted LD reference files are named LDscore.<chr>.l2.ldscore.gz;
    # the task strips this prefix via a symlink farm before calling ldsc.
    config=GWASBySubtractionFullPythonConfig(ld_file_filename_pattern="LDscore."),
)
