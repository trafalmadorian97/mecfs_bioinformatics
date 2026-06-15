"""GWAS-by-subtraction of Johnston multisite pain from DecodeME GWAS-1, with
DecodeME treated as a LINEAR (OLS) trait instead of logistic.

Why this variant exists
-----------------------
The standard task (`decode_me_minus_pain_subtraction.py`) marks DecodeME as
`logistic` and supplies case/control prevalences. That puts the two halves of
the subtraction on inconsistent scales:

- the subtraction coefficient c = S12/S22 comes from ldsc(), where DecodeME's
  genetic covariance is converted to the population *liability* scale; but
- the per-SNP betas come from sumstats(), where `se.logit` standardises DecodeME
  to the unit-variance *logistic latent* scale (beta = logOR / sqrt(logOR^2*varSNP
  + pi^2/3)).

Because the projection coefficient and the betas it is applied to live on
different scales, the remainder factor R is only approximately orthogonal to
pain -- ct-LDSC shows rg(R, pain) ~ 0.23 instead of 0 (see
experiments/claude/residual_rg_*.py).

This variant removes the mismatch by keeping DecodeME on a single linear scale
end to end:

- `gwas_method="ols"`: the OLS branch of sumstats() standardises the effect as
  beta = Z / sqrt(N*varSNP) -- a purely linear map with no pi^2/3 denominator
  (unlike `logistic` and even `linear_prob`, both of which apply that nonlinear
  term).
- `sample_info=QuantPhenotype(...)`: no prevalences, so ldsc() does NOT apply the
  liability conversion and S stays on the observed linear scale, consistent with
  the OLS betas.

With both halves on the same observed-linear scale, c = S12/S22 is the genetic
regression slope of the very betas being subtracted, so R is genuinely
orthogonal to pain. Genetic-correlation / Z-scores / p-values of R are
scale-invariant, so this does not distort which variants are significant for the
ME/CFS-specific factor; it only fixes the orthogonalisation.
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
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    QuantPhenotype,
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

# DecodeME GWAS-1 total sample size (15,579 cases + 259,909 controls). Carried
# through as the linear-trait N; no prevalences, so ldsc keeps S on the observed
# (non-liability) scale.
_DECODE_ME_GWAS_1_TOTAL_N = 259_909 + 15_579

DECODE_ME_MINUS_PAIN_SUBTRACTION_OLS = GenomicSEMGWASBySubtractionFullPythonTask.create(
    asset_id="decode_me_minus_pain_subtraction_full_python_ols",
    # Composite trait (T1): DecodeME, treated as a linear (OLS) trait so its
    # standardised betas share a scale with the genetic-covariance matrix.
    composite_trait_source=GenomicSEMGWASSumstatsSource(
        source=GenomicSEMSumstatsSource(
            task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS.join_task,
            alias="DECODE_ME",
            # QuantPhenotype (not BinaryPhenotypeSampleInfo): supplies N but
            # no prevalences, so ldsc does not liability-convert DecodeME.
            sample_info=QuantPhenotype(total_sample_size=_DECODE_ME_GWAS_1_TOTAL_N),
            pipe=CompositePipe(
                [
                    RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
                    ComputePPipe(),
                ]
            ),
        ),
        gwas_method="ols",
    ),
    # Reference trait (T2): Johnston multisite chronic pain (OLS), unchanged.
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
    config=GWASBySubtractionFullPythonConfig(ld_file_filename_pattern="LDscore."),
)
