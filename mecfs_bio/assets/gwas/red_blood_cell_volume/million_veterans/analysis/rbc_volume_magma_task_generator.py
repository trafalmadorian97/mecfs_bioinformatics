"""
Generate the tasks required to apply MAGMA to the red blood cell volume GWAS
"""

from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.red_blood_cell_volume.million_veterans.raw_gwas_data.raw_red_blood_volume import (
    MILLION_VETERAN_RED_BLOOD_VOLUME_EUR_DATA_RAW_HARMONIZED,
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
from mecfs_bio.assets.reference_data.rna_seq_data.gtex.gtex_v10_median_tissue_expression_rna_seq_prep_for_magma import (
    GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.pipes.filter_rows_by_info_score import (
    FilterRowsByInfoScorePipe,
)
from mecfs_bio.build_system.task_generator.magma_task_generator import (
    MagmaTaskGeneratorFromRaw,
)

MILLION_VETERANS_EUR_RBC_VOLUME_MAGMA_TASKS = MagmaTaskGeneratorFromRaw.create(
    raw_gwas_data_task=MILLION_VETERAN_RED_BLOOD_VOLUME_EUR_DATA_RAW_HARMONIZED,
    base_name="million_veterans_eur_rbc_vol",
    sample_size=407294,
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    gene_loc_file_task=MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
    magma_ld_ref_task=MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
    tissue_expression_gene_set_task=GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA,
    fmt=GWASLabColumnSpecifiers(
        rsid="rsid",
        # note: RSIDS do not change with genome build, so this column is valid regardless of any downstream genome build conversion
        OR=None,
        se="standard_error",
        chrom="chromosome",
        pos="base_pair_location",
        ea="effect_allele",
        nea="other_allele",
        eaf="effect_allele_frequency",
        beta="beta",
        p="p_value",
        n="n",
        info="r2",
        snpid=None,
    ),
    pre_pipe=FilterRowsByInfoScorePipe(min_score=0.8, info_col="r2"),
    gene_thesaurus_task=GENE_THESAURUS,
)
