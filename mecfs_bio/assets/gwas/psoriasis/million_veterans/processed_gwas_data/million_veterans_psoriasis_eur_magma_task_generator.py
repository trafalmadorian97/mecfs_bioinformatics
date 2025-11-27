from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.psoriasis.million_veterans.raw_gwas_data.raw_psoriasis_data import (
    MILLION_VETERAN_PSORIASIS_EUR_DATA_RAW,
)
from mecfs_bio.assets.reference_data.db_snp.snp151_build37_parquet import (
    GENOME_ANNOTATION_DATABASE_BUILD_37_PARQUET,
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
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.filter_rows_by_info_score import (
    FilterRowsByInfoScorePipe,
)
from mecfs_bio.build_system.task_generator.magma_task_generator import (
    MagmaTaskGeneratorFromRawCompute37RSIDs,
)

MILLION_VETERANS_EUR_PSORIASIS_COMBINED_MAGMA_TASKS = MagmaTaskGeneratorFromRawCompute37RSIDs.create(
    raw_gwas_data_task=MILLION_VETERAN_PSORIASIS_EUR_DATA_RAW,
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    gene_loc_file_task=MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
    magma_ld_ref_task=MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
    tissue_expression_gene_set_task=GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA,
    base_name="million_veterans_eur_psoriasis",
    sample_size=443794,  # source: https://www.ebi.ac.uk/gwas/studies/GCST90476186
    fmt=GWASLabColumnSpecifiers(
        # rsid="rsid",
        rsid=None,
        chrom="chromosome",
        pos="base_pair_location",
        ea="effect_allele",
        nea="other_allele",
        eaf="effect_allele_frequency",
        p="p_value",
        OR="odds_ratio",
        # se="standard_error",
        se=None,
        n="n",
        ncase="num_cases",
        ncontrol="num_controls",
        info="r2",
        # snpid="variant_id",
        snpid=None,
        or_95l="ci_lower",
        or_95u="ci_upper",
    ),
    pipe=CompositePipe(
        pipes=[
            # FilterRowsByInfoScorePipe(min_score=0.9),
            ComputeBetaPipe()
        ]
    ),
    pre_pipe=FilterRowsByInfoScorePipe(min_score=0.9, info_col="r2"),
    genome_build="infer",  # "38", # source: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476186/GCST90476186.tsv.gz-meta.yaml
    snp151_database_file_task=GENOME_ANNOTATION_DATABASE_BUILD_37_PARQUET,
)
