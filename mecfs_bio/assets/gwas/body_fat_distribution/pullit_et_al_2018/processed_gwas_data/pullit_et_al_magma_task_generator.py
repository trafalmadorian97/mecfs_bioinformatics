"""
Task generator to apply the main steps of MAGMA analysis is Pullit et al's
Body fat distribution GWAS
"""

from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.body_fat_distribution.pullit_et_al_2018.raw_gwas_data.pullit_et_al_raw_gwas_data import (
    PULLIT_ET_AL_2018_RAW,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.raw.magma_ensembl_gene_location_reference_data_build_37 import (
    MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
)
from mecfs_bio.assets.reference_data.magma_ld_reference.magma_eur_build_37_1k_genomes_ref_extracted import (
    MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
)
from mecfs_bio.assets.reference_data.rna_seq_data.gtex_v10_median_tissue_expression_rna_seq_prep_for_magma import (
    GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.drop_null_pipe import DropNullsPipe
from mecfs_bio.build_system.task.pipes.str_split_col import SplitColPipe
from mecfs_bio.build_system.task_generator.magma_task_generator import (
    MagmaTaskGeneratorFromRaw,
)

PULLIT_ET_AL_2018_COMBINED_MAGMA_TASKS = MagmaTaskGeneratorFromRaw.create(
    raw_gwas_data_task=PULLIT_ET_AL_2018_RAW,
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    gene_loc_file_task=MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
    magma_ld_ref_task=MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
    tissue_expression_gene_set_task=GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA,
    base_name="pullit_et_al_2018",
    sample_size=484_564,  # source: data file
    fmt=GWASLabColumnSpecifiers(
        # rsid="SNP",
        chrom="CHR",
        pos="POS",
        ea="Tested_Allele",
        nea="Other_Allele",
        eaf="Freq_Tested_Allele",
        beta="BETA",
        p="P",
        snpid="SNP",
        OR=None,
        se="SE",
        info="INFO",
        rsid="RSID",
    ),
    pre_pipe=CompositePipe(
        [
            DropNullsPipe(subset=["CHR", "POS"]),
            SplitColPipe(
                col_to_split="SNP", split_by=":", new_col_names=("RSID", "al1", "al2")
            ),
        ]
    ),
)
