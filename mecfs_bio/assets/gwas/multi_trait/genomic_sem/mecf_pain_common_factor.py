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
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_common_factor_gwas_task import (
    GenomicSEMCommonFactorGWASTask,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_task import (
    GenomicSEMConfig,
    GenomicSEMSumstatsSource,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_user_gwas_task import (
    GenomicSEMGWASRunConfig,
    GenomicSEMGWASSumstatsSource,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL

MECFS_PAIN_COMMON_FACTOR = GenomicSEMCommonFactorGWASTask.create(
    asset_id="decode_me_pain_common_factor_gwas",
    sources=[
        GenomicSEMGWASSumstatsSource(
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
        GenomicSEMGWASSumstatsSource(
            source=GenomicSEMSumstatsSource(
                task=JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
                alias="Multsite_Pain",
                sample_info=JOHNSTON_ET_AL_SAMPLE_INFO,
            ),
            gwas_method="ols",
        ),
    ],
    ld_ref_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
    hapmap_snps_task=HAPMAP3_SNPLIST_FOR_GENOMIC_SEM,
    sumstats_ref_task=GENOMES1K_REFERENCE_FOR_GENOMIC_SEM,
    # The extracted LD reference files are named LDscore.<chr>.l2.ldscore.gz;
    # the task strips this prefix via a symlink farm before calling ldsc.
    munge_config=GenomicSEMConfig(ld_file_filename_pattern="LDscore."),
    # commonfactorGWAS forks one R worker per core (default = nproc-1). Each
    # worker holds ~5GB of state with this dataset, so on a 15GB WSL2 box the
    # default explodes the memory budget and the kernel OOM-kills workers
    # mid-run (surfacing as "error reading from connection"). Cap cores so
    # cores * 5GB stays under the system's available memory.
    run_config=GenomicSEMGWASRunConfig(cores=2),
)
