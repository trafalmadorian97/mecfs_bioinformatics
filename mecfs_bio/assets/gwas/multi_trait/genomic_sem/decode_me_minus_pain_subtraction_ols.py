"""GWAS-by-subtraction of Johnston multisite pain from DecodeME GWAS-1.

DecodeME is the composite trait (T1); Johnston pain is the reference trait (T2),
a pure indicator of the common factor F. The remainder factor R is the
ME/CFS-specific signal left after the pain-shared component is subtracted out.

Both traits are standardised on a single linear scale and no liability
conversion is applied -- this is the only coherent treatment (see
`genomic_sem_gwas_by_subtraction_full_python_task` for why DecodeME being a
case/control trait is *not* given binary/logistic special-casing: doing so would
leave R non-orthogonal to pain). The `_ols` suffix is historical; it is now the
standard subtraction.
"""

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
    GenomicSEMSumstatsSource,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_gwas_by_subtraction_full_python_task import (
    GenomicSEMGWASBySubtractionFullPythonTask,
    GWASBySubtractionFullPythonConfig,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL

# DecodeME GWAS-1 total sample size (15,579 cases + 259,909 controls), used to
# linearly standardise the per-SNP effects.
_DECODE_ME_GWAS_1_TOTAL_N = 259_909 + 15_579

DECODE_ME_MINUS_PAIN_SUBTRACTION_OLS = GenomicSEMGWASBySubtractionFullPythonTask.create(
    asset_id="decode_me_minus_pain_subtraction_full_python_ols",
    # Composite trait (T1): DecodeME.
    composite_trait_source=GenomicSEMSumstatsSource(
        task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS.join_task,
        alias="DECODE_ME",
        sample_size=_DECODE_ME_GWAS_1_TOTAL_N,
        pipe=CompositePipe(
            [
                RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
                ComputePPipe(),
            ]
        ),
    ),
    # Reference trait (T2): Johnston multisite chronic pain.
    reference_trait_source=GenomicSEMSumstatsSource(
        task=JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
        alias="Multsite_Pain",
        sample_size=JOHNSTON_ET_AL_SAMPLE_INFO.total_sample_size,
    ),
    ld_ref_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
    hapmap_snps_task=HAPMAP3_SNPLIST_FOR_GENOMIC_SEM,
    sumstats_ref_task=GENOMES1K_REFERENCE_FOR_GENOMIC_SEM,
    config=GWASBySubtractionFullPythonConfig(ld_file_filename_pattern="LDscore."),
)
