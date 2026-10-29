"""
Asset generator to apply standard analysis to Willer et al.'s GWAS of LDL.

The analysis uses the GWAS Catalog's harmonised alleles (hm_effect_allele, hm_other_allele,
hm_beta). The file's only allele frequency, eaf_ref, is the 1000 Genomes EUR frequency of the
ORIGINAL effect_allele, which the Catalog swapped for about half the variants. So eaf_ref is
re-expressed as the frequency of hm_effect_allele before use. Whether the Catalog swapped a
variant is read from the sign of its beta: eaf_ref is kept where hm_beta == beta, replaced by
1 - eaf_ref where hm_beta == -beta, and null otherwise (including beta == 0 and rows the
Catalog could not harmonise).

Palindromic SNVs (A/T, C/G) get a null frequency. Their eaf_ref is a placeholder, not a
frequency: it is null for about half of them and exactly 1.0 for 93% of the rest, and it
agrees with neither allele's 1000 Genomes EUR frequency.
"""

import narwhals

from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.ldl.willer_et_al.raw.raw_willer_ldl_data import (
    WILLER_LDL_EUR_DATA_RAW,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    QuantPhenotype,
)
from mecfs_bio.build_system.task.pipes.expr_pipe import ExprPipe

_HM_EAF_COL = "hm_oriented_eaf_ref"

_BETA = narwhals.col("beta")
_HM_BETA = narwhals.col("hm_beta")
_EAF_REF = narwhals.col("eaf_ref")
_HM_EA = narwhals.col("hm_effect_allele")
_HM_NEA = narwhals.col("hm_other_allele")
_IS_PALINDROMIC = (
    ((_HM_EA == "A") & (_HM_NEA == "T"))
    | ((_HM_EA == "T") & (_HM_NEA == "A"))
    | ((_HM_EA == "C") & (_HM_NEA == "G"))
    | ((_HM_EA == "G") & (_HM_NEA == "C"))
)
HM_ORIENTED_EAF_PIPE = ExprPipe(
    narwhals.when(~_IS_PALINDROMIC & (_BETA != 0) & (_HM_BETA == _BETA))
    .then(_EAF_REF)
    .otherwise(
        narwhals.when(~_IS_PALINDROMIC & (_BETA != 0) & (_HM_BETA == -_BETA)).then(
            1 - _EAF_REF
        )
    )
    .alias(_HM_EAF_COL)
)

WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS = (
    concrete_standard_analysis_generator_assume_already_has_rsid(
        base_name="willer_et_al_2023_ldl_eur",
        raw_gwas_data_task=WILLER_LDL_EUR_DATA_RAW,
        fmt=GWASLabColumnSpecifiers(
            rsid="hm_rsid",
            chrom="hm_chrom",
            pos="hm_pos",
            nea="hm_other_allele",
            ea="hm_effect_allele",
            beta="hm_beta",
            se="standard_error",
            p="p_value",
            eaf=_HM_EAF_COL,
        ),
        sample_size=89138,
        phenotype_info_for_ldsc=QuantPhenotype(),
        pre_pipe=HM_ORIENTED_EAF_PIPE,
    )
)
