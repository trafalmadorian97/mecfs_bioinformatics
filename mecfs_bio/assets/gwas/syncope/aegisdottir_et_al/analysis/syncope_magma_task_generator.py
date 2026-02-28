"""
Generate MAGMA tasks for the Syncope GWAS study of AegisDottir et al.

Looks like the data is already filtered to include only SNPs with imputation quality of at least 0.9
So no further filtering is required

"""

from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.processed.sumstats_liftover_dump_to_parquet import \
    AEGISDOTTIR_SUMSTATS_DUMP_TO_PARQUET
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.raw.raw_syncope_data import (
    AEGISDOTTIR_ET_AL_RAW_SYNCOPE_DATA,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.raw.magma_ensembl_gene_location_reference_data_build_37 import (
    MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
)
from mecfs_bio.assets.reference_data.magma_ld_reference.magma_eur_build_37_1k_genomes_ref_extracted import (
    MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
)
from mecfs_bio.assets.reference_data.rna_seq_data.gtex.gtex_v10_median_tissue_expression_rna_seq_prep_for_magma import \
    GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task_generator.magma_task_generator import (
    MagmaTaskGeneratorFromRaw,
)

AEGISDOTTIR_COMBINED_MAGMA_TASKS = MagmaTaskGeneratorFromRaw.create(
    raw_gwas_data_task=AEGISDOTTIR_SUMSTATS_DUMP_TO_PARQUET,
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    gene_loc_file_task=MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
    magma_ld_ref_task=MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
    tissue_expression_gene_set_task=GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA,
    base_name="aegisdottir_et_al_2023",
    sample_size=946_861,
    fmt="gwaslab",
    post_pipe=ComputeBetaPipe(),
)
