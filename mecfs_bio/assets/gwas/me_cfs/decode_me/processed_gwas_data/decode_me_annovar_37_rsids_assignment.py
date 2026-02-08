"""
Generate a chain of tasks to assign rsids to the decodeME genetic variants using dbSNP 150 reference data as pre-processed by annovar.
Note that that unlike the approach based on ucsc hg19 sbSNP, this approach does allow assignment of rsids to non-SNP genetic variants.

There are two versions of the SNP-assigned data.  One discards strand-ambiguous palindromic variants, and one keeps them.
The version that keeps ambiguous variants is suitable for downstream tasks like MAGMA, which ignore variant direction of effect.
On the other hand, it should not be used for downstream tasks like Mendelian Randomization, which crucially depend on direction of effect.
"""

from mecfs_bio.asset_generator.annovar_37_basic_rsid_assignment import (
    annovar_37_basic_rsid_assignment,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_sumstats_liftover_to_37 import (
    DECODE_ME_GWAS_1_SUMSTATS_LIFTOVER_TO_37,
)

DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED = annovar_37_basic_rsid_assignment(
    sumstats_task=DECODE_ME_GWAS_1_SUMSTATS_LIFTOVER_TO_37,
    base_name="decode_me_gwas_1",
    drop_palindromic_ambiguous=True,
)


DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS = (
    annovar_37_basic_rsid_assignment(
        sumstats_task=DECODE_ME_GWAS_1_SUMSTATS_LIFTOVER_TO_37,
        base_name="decode_me_gwas_1_keep_ambiguous",
        drop_palindromic_ambiguous=False,
    )
)
