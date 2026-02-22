"""
Task generator to apply S-LDSC to the European results from the Inflammatory Bowel Disease study of Liu et al, using
rsids assigned via SNP150
"""

from mecfs_bio.asset_generator.concrete_sldsc_generator import (
    standard_sldsc_task_generator,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_37_harmonized_with_rsid_sumstats import (
    LIU_ET_AL_SUMSTATS_WITH_RSID_FROM_SNP150,
)
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_SAMPLE_SIZE_COLUMN,
)

LIU_ET_AL_S_LSDC_FROM_SNP_150 = standard_sldsc_task_generator(
    sumstats_task=LIU_ET_AL_SUMSTATS_WITH_RSID_FROM_SNP150,
    base_name="liu_et_al_ibd_from_snp150",
    pre_pipe=SetColToConstantPipe(
        col_name=GWASLAB_SAMPLE_SIZE_COLUMN,
        constant=59957,  # source: https://pmc.ncbi.nlm.nih.gov/articles/PMC10290755/pdf/nihms-1903722.pdf
    ),
)
