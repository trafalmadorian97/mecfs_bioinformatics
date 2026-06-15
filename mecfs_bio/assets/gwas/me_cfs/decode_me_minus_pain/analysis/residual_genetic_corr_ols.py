"""Genetic correlation of the OLS-variant remainder factor with DecodeME and
multisite pain.

Mirrors `residual_genetic_corr.py` but uses the linear-scale subtraction output
(`..._ols`). If the scale-consistency fix works, rg(remainder, pain) should be
~0 here, versus ~0.23 for the logistic variant.
"""

from mecfs_bio.asset_generator.genetic_correlation_asset_generator import (
    genetic_corr_by_ct_ldsc_asset_generator,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.auxiliary.prevalance_info import (
    DECODE_ME_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_with_annovar_37_rsids_sumstats import (
    DECODEME_ME_SUMSTATS_37_WITH_ANNOVAR_RSID,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_minus_pain.analysis.standard_analysis_decodeme_minus_pain_ols import (
    DECODE_ME_MINUS_PAIN_OLS_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import (
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    QuantPhenotype,
    SumstatsSource,
)
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_SAMPLE_SIZE_COLUMN

DECODE_ME_MINUS_PAIN_OLS_GENETIC_CORR_GENERATOR = (
    genetic_corr_by_ct_ldsc_asset_generator(
        base_name="decode_me_minus_pain_ols_genetic_corr",
        sources=[
            SumstatsSource(
                DECODE_ME_MINUS_PAIN_OLS_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
                alias="DecodeME_Minus_Pain_OLS",
                sample_info=QuantPhenotype(),
            ),
            SumstatsSource(
                DECODEME_ME_SUMSTATS_37_WITH_ANNOVAR_RSID,
                alias="DecodeME",
                pipe=SetColToConstantPipe(
                    GWASLAB_SAMPLE_SIZE_COLUMN,
                    constant=275488,
                ),  # true total sample size. From preprint
                sample_info=DECODE_ME_PREVALENCE_INFO,
            ),
            SumstatsSource(
                JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
                alias="Multisite_pain",
                pipe=SetColToConstantPipe(
                    GWASLAB_SAMPLE_SIZE_COLUMN, constant=387649
                ),  # True total sample size. From Gwas catalog
                sample_info=QuantPhenotype(),
            ),
        ],
    )
)
