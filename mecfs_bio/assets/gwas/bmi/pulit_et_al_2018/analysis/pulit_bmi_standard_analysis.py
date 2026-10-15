"""
Standard analysis of the GIANT + UK Biobank meta-analysis GWAS of body mass index from Pulit et al. (2019).

The SNP column has the form rsid:allele1:allele2 (or a bare rsid for some variants), so the rsid is
split out of it.  Sample size varies by variant (median 484,680; up to 806,834 where GIANT and
UK Biobank overlap), so the per-variant N column is used.
"""

import narwhals

from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.bmi.pulit_et_al_2018.raw.pulit_bmi_raw import (
    PULIT_ET_AL_2018_BMI_RAW,
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
from mecfs_bio.build_system.task.pipes.str_split_exact_col import SplitExactColPipe

PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS = (
    concrete_standard_analysis_generator_assume_already_has_rsid(
        base_name="pulit_et_al_2018_bmi",
        raw_gwas_data_task=PULIT_ET_AL_2018_BMI_RAW,
        fmt=GWASLabColumnSpecifiers(
            rsid="RSID",
            snpid="SNP",
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
                SplitExactColPipe(
                    col_to_split="SNP",
                    split_by=":",
                    new_col_names=("RSID", "al1", "al2"),
                ),
            ]
        ),
        phenotype_info_for_ldsc=QuantPhenotype(),
    )
)
