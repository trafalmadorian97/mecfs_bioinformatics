"""
Standard analysis of the GIANT + UK Biobank meta-analysis GWAS of body mass index from Pulit et al. (2019).

Variants are first restricted to MAF >= 0.01.

The SNP column has the form rsid:allele1:allele2 (or a bare rsid for some variants), so the rsid
is taken from its prefix.  The rsid is extracted lazily and the SNP string column is not passed
to gwaslab, which keeps peak memory within a 16GB machine for this large file.

Sample size varies by variant (median 484,680; up to 806,834 where GIANT and UK Biobank
overlap), so the per-variant N column is used.
"""

import narwhals

from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.bmi.pulit_et_al_2018.processed.pulit_bmi_common_snps import (
    PULIT_ET_AL_2018_BMI_COMMON_SNPS,
)
from mecfs_bio.build_system.sample_size_spec import PerVariantSampleSize
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    QuantPhenotype,
)
from mecfs_bio.build_system.task.pipes.cast_pipe import CastPipe
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.drop_null_pipe import DropNullsPipe
from mecfs_bio.build_system.task.pipes.str_replace_pipe import StrReplacePipe

PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS = (
    concrete_standard_analysis_generator_assume_already_has_rsid(
        base_name="pulit_et_al_2018_bmi",
        raw_gwas_data_task=PULIT_ET_AL_2018_BMI_COMMON_SNPS,
        fmt=GWASLabColumnSpecifiers(
            rsid="RSID",
            chrom="CHR",
            pos="POS",
            ea="Tested_Allele",
            nea="Other_Allele",
            eaf="Freq_Tested_Allele",
            beta="BETA",
            se="SE",
            p="P",
            info="INFO",
            n="N",
            OR=None,
        ),
        sample_size=PerVariantSampleSize(),
        pre_pipe=CompositePipe(
            [
                DropNullsPipe(subset=["CHR", "POS", "N"]),
                CastPipe(target_column="N", type=narwhals.Int64(), new_col_name="N"),
                StrReplacePipe(
                    target_column="SNP",
                    new_column_name="RSID",
                    replace_what=":.*$",
                    replace_with="",
                ),
            ]
        ),
        phenotype_info_for_ldsc=QuantPhenotype(),
    )
)
