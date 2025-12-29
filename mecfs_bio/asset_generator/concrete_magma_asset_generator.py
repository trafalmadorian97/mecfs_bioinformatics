from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.reference_data.ensembl_biomart.gene_thesaurus import (
    GENE_THESAURUS,
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
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    ValidGwaslabFormat,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task_generator.magma_task_generator import (
    MagmaTaskGeneratorFromRaw,
)


def concrete_magma_assets_generate(
    base_name: str,
    raw_gwas_data_task: Task,
    fmt: ValidGwaslabFormat,
    sample_size: int,
    pre_pipe: DataProcessingPipe = IdentityPipe(),
) -> MagmaTaskGeneratorFromRaw:
    return MagmaTaskGeneratorFromRaw.create(
        raw_gwas_data_task=raw_gwas_data_task,
        base_name=base_name,
        sample_size=sample_size,
        magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
        gene_loc_file_task=MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
        magma_ld_ref_task=MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
        tissue_expression_gene_set_task=GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA,
        fmt=fmt,
        pre_pipe=pre_pipe,
        gene_thesaurus_task=GENE_THESAURUS,
    )
