"""
Constants for interacting with gwaslab (https://github.com/Cloufield/gwaslab/tree/main/src/gwaslab).
"""

from collections.abc import Mapping
from pathlib import PurePath
from typing import Literal

GWASLAB_CHROM_COL = "CHR"
GWASLAB_CHISQ_COL = "CHISQ"
GWASLAB_POS_COL = "POS"
GWASLAB_N_CASE_COL = "N_CASE"
GWASLAB_N_CONTROL_COL = "N_CONTROL"

GWASLAB_EFFECT_ALLELE_COL = "EA"
GWASLAB_NON_EFFECT_ALLELE_COL = "NEA"
GWASLAB_INFO_SCORE_COL = "INFO"
GWASLAB_ODDS_RATIO_COL = "OR"
GWASLAB_BETA_COL = "BETA"
GWASLAB_P_COL = "P"
GWASLAB_MLOG10P_COL = "MLOG10P"
GWASLAB_EFFECT_ALLELE_FREQ_COL = "EAF"
GWASLAB_SAMPLE_SIZE_COLUMN = "N"
GWASLAB_SE_COL = "SE"
GWASLAB_STATUS_COL = "STATUS"
GWASLAB_RSID_COL = "rsID"
GWASLAB_SNPID_COL = "SNPID"

GwaslabKnownFormat = Literal["gwaslab", "regenie", "gwascatalog"]

GWASLabVCFRefFile = Literal["1kg_eur_hg38", "1kg_eur_hg19"]

GWASLAB_EUR_1K_GENOMES_NAME_38: GWASLabVCFRefFile = "1kg_eur_hg38"
GWASLAB_HUMAN_GENOME_NAME_38 = "ucsc_genome_hg38"

GWASLAB_EUR_1K_GENOMES_NAME_19 = "1kg_eur_hg19"
GWASLAB_HUMAN_GENOME_NAME_19 = "ucsc_genome_hg19"

# gwaslab download-ref key for the 1000G dbSNP151 hg38 autosomal variant table
# (columns SNPID, rsID, CHR, POS, NEA, EA); used to annotate rsIDs.
GWASLAB_1KG_DBSNP151_HG38_AUTO_NAME = "1kg_dbsnp151_hg38_auto"

# HapMap3 hg38 snplist bundled inside the gwaslab package (columns rsid, A1, A2,
# #CHROM, POS), relative to the gwaslab package directory. PurePath (not str) since
# it is a relative path to be resolved against an as-yet-unspecified base directory.
GWASLAB_HAPMAP3_HG38_SNPLIST_RELPATH = PurePath(
    "data/hapmap3_SNPs/hapmap3_db151_hg38.snplist.gz"
)

GWASLAB_EFFECTIVE_SAMPLE_SIZE = "N_EFF"

# gwaslab codes the sex chromosomes and the mitochondrion numerically, X as 23, Y as
# 24 and MT as 25 (see get_chr_to_number in gwaslab.bd.bd_common_data), and we follow
# it so that CHR stays an integer column. Reference VCFs and the UKB-PPP per-chromosome
# filenames spell them with letters instead, so these translate at that boundary.
# Note UKB-PPP does both at once: its file is named chrX while the CHROM inside it is
# already 23.
GWASLAB_CHROM_CODE_FOR_NAME: Mapping[str, int] = {"X": 23, "Y": 24, "MT": 25}
GWASLAB_CHROM_NAME_FOR_CODE: Mapping[int, str] = {
    code: name for name, code in GWASLAB_CHROM_CODE_FOR_NAME.items()
}

# Further gwaslab standard column names, from the "gwaslab" entry of
# https://github.com/Cloufield/formatbook/blob/main/formatbook.json
GWASLAB_Z_COL = "Z"
GWASLAB_T_STATISTIC_COL = "T"
GWASLAB_F_STATISTIC_COL = "F"
GWASLAB_HAZARD_RATIO_COL = "HR"
GWASLAB_ODDS_RATIO_95L_COL = "OR_95L"
GWASLAB_ODDS_RATIO_95U_COL = "OR_95U"
GWASLAB_HAZARD_RATIO_95L_COL = "HR_95L"
GWASLAB_HAZARD_RATIO_95U_COL = "HR_95U"
GWASLAB_P_HET_COL = "P_HET"
GWASLAB_I2_COL = "I2"
GWASLAB_SNPR2_COL = "SNPR2"
GWASLAB_DOF_COL = "DOF"
GWASLAB_MAF_COL = "MAF"
