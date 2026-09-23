"""DecodeME-specific per-SNP heritability prior (PolyFun Approach 2).

An L2-regularized S-LDSC over the baseline-LF annotations, fit on the DecodeME
build-37 sumstats themselves (genome-reference harmonized, so NEA is the hg19 REF
allele and joins exactly with the REF-oriented baseline-LF tables). Used as an
alternative to the precomputed 15-trait PolyFun prior, to check whether the
generic prior mismatches DecodeME.
"""

from pathlib import PurePath

from mecfs_bio.asset_generator.polyfun_explain_fine_mapping_asset_generator import (
    PolyfunPriorSource,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.assets.reference_data.polyfun.annotations.baseline_lf_annotations import (
    BASELINE_LF_ANNOTATION_MATRIX,
)
from mecfs_bio.assets.reference_data.polyfun.annotations.baseline_lf_ldscores import (
    BASELINE_LF_ANNOTATION_LDSCORE_MEMBERS,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.task.annotation_weights.l2_sldsc_snpvar_task import (
    SNPVAR_COL,
    SNPVAR_PARQUET_FILENAME,
    TAU_EVEN_WEIGHTS_FILENAME,
    TAU_ODD_WEIGHTS_FILENAME,
    L2RegularizedSldscSnpvarTask,
)
from mecfs_bio.build_system.task.copy_file_from_directory_task import (
    CopyFileFromDirectoryTask,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

# 4 / (1/cases + 1/controls), as in the DecodeME fine-mapping modules.
DECODE_ME_EFFECTIVE_SAMPLE_SIZE = int(4 / (1 / 15_579 + 1 / 259_909))

DECODE_ME_L2_SLDSC_SNPVAR = L2RegularizedSldscSnpvarTask.create(
    asset_id="decode_me_gwas_1_l2_sldsc_snpvar",
    sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.join_task,
    effective_sample_size=DECODE_ME_EFFECTIVE_SAMPLE_SIZE,
    annotation_ldscore_members_task=BASELINE_LF_ANNOTATION_LDSCORE_MEMBERS,
    annotation_matrix_task=BASELINE_LF_ANNOTATION_MATRIX,
)

DECODE_ME_L2_SLDSC_SNPVAR_TABLE = CopyFileFromDirectoryTask.create_result_table(
    asset_id="decode_me_gwas_1_l2_sldsc_snpvar_table",
    source_directory_task=DECODE_ME_L2_SLDSC_SNPVAR,
    path_inside_directory=PurePath(SNPVAR_PARQUET_FILENAME),
    extension=".parquet",
    read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
)

# tau fit on the odd chromosomes (it scores, so explains, the even chromosomes).
DECODE_ME_L2_SLDSC_TAU_ODD_WEIGHTS = CopyFileFromDirectoryTask.create_result_table(
    asset_id="decode_me_gwas_1_l2_sldsc_tau_odd_weights",
    source_directory_task=DECODE_ME_L2_SLDSC_SNPVAR,
    path_inside_directory=PurePath(TAU_ODD_WEIGHTS_FILENAME),
    extension=".parquet",
    read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
)

# tau fit on the even chromosomes (it scores, so explains, the odd chromosomes).
DECODE_ME_L2_SLDSC_TAU_EVEN_WEIGHTS = CopyFileFromDirectoryTask.create_result_table(
    asset_id="decode_me_gwas_1_l2_sldsc_tau_even_weights",
    source_directory_task=DECODE_ME_L2_SLDSC_SNPVAR,
    path_inside_directory=PurePath(TAU_EVEN_WEIGHTS_FILENAME),
    extension=".parquet",
    read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
)

DECODE_ME_L2_SLDSC_PRIOR_SOURCE = PolyfunPriorSource(
    prior_task=DECODE_ME_L2_SLDSC_SNPVAR_TABLE,
    weight_col=SNPVAR_COL,
    odd_chrom_explanation_weights_task=DECODE_ME_L2_SLDSC_TAU_EVEN_WEIGHTS,
    even_chrom_explanation_weights_task=DECODE_ME_L2_SLDSC_TAU_ODD_WEIGHTS,
    chr_col=GWASLAB_CHROM_COL,
    pos_col=GWASLAB_POS_COL,
    nea_col=GWASLAB_NON_EFFECT_ALLELE_COL,
    ea_col=GWASLAB_EFFECT_ALLELE_COL,
)
