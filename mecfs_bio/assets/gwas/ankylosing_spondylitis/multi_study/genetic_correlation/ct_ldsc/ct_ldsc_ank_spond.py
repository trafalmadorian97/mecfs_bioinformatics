from mecfs_bio.asset_generator.genetic_correlation_asset_generator import (
    genetic_corr_by_ct_ldsc_asset_generator,
)
from mecfs_bio.assets.gwas.ankylosing_spondylitis.finngen.auxiliary.finngen_ank_spond_phenotype_info import (
    FINNGEN_ANK_SPOND_PHENOTYPE_INFO,
)
from mecfs_bio.assets.gwas.ankylosing_spondylitis.finngen.processed.finngen_ank_spond_sumstats_harmonized import (
    FINNGEN_ANK_SPOND_SUMSTATS_HARMONIZED,
)
from mecfs_bio.assets.gwas.ankylosing_spondylitis.million_veterans.auxiliary.million_veterans_ank_spond_phenotype_info import (
    MILLION_VETERANS_ANK_SPOND_PHENOTYPE_INFO,
)
from mecfs_bio.assets.gwas.ankylosing_spondylitis.million_veterans.processed.mv_eur_ank_spond_sumstats_harmonized import (
    MILLION_VETERANS_ANK_SPOND_SUMSTATS_HARMONIZED,
)
from mecfs_bio.assets.gwas.ankylosing_spondylitis.ukbb.auxiliary.ukbb_ank_spond_phenotype_info import (
    UKBB_ANK_SPOND_PHENOTYPE_INFO,
)
from mecfs_bio.assets.gwas.ankylosing_spondylitis.ukbb.processed.ukbb_ank_spond_sumstats_harmonized import (
    UKBB_ANK_SPOND_SUMSTATS_HARMONIZED,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    SumstatsSource,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaIfNeededPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEIfNeededPipe
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

pipe = CompositePipe(
    [
        UniquePipe(  # Most the duplications in this gwas appear to be due to anomalies caused by conversion of multi-allelic to bi-allelic variants
            by=[
                GWASLAB_CHROM_COL,
                GWASLAB_POS_COL,
                GWASLAB_EFFECT_ALLELE_COL,
                GWASLAB_NON_EFFECT_ALLELE_COL,
            ],
            keep="none",
            order_by=[GWASLAB_POS_COL],
        ),
        ComputeBetaIfNeededPipe(),
        ComputeSEIfNeededPipe(),
    ]
)

ANK_SPOND_CT_LDSC_ASSET_GENERATOR = genetic_corr_by_ct_ldsc_asset_generator(
    base_name="ank_spond_ct_ldsc",
    sources=[
        SumstatsSource(
            UKBB_ANK_SPOND_SUMSTATS_HARMONIZED,
            alias="UKBB_Ank_Spond",
            sample_info=UKBB_ANK_SPOND_PHENOTYPE_INFO,
            pipe=pipe,
        ),
        SumstatsSource(
            MILLION_VETERANS_ANK_SPOND_SUMSTATS_HARMONIZED,
            alias="Million_Veterans_Ank_Spond",
            sample_info=MILLION_VETERANS_ANK_SPOND_PHENOTYPE_INFO,
            pipe=pipe,
        ),
        SumstatsSource(
            FINNGEN_ANK_SPOND_SUMSTATS_HARMONIZED,
            alias="FinnGen_Ank_Spond",
            sample_info=FINNGEN_ANK_SPOND_PHENOTYPE_INFO,
            pipe=pipe,
        ),
    ],
)
