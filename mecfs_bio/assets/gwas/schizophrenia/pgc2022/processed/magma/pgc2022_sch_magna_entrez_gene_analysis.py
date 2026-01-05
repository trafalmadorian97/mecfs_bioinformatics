from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.magma.pgc2022_sch_entrez_annotations import (
    PGC2022_SCH_MAGMA_ENTREZ_ANNOTATIONS,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.standard_analysis_sc_pgc_2022 import (
    SCH_PGC_2022_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.reference_data.magma_ld_reference.magma_eur_build_37_1k_genomes_ref_extracted import (
    MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
)
from mecfs_bio.build_system.task.magma.magma_gene_analysis_task import (
    MagmaGeneAnalysisTask,
)

PGC2022_SCH_MAGMA_ENTREZ_GENE_ANALYSIS = MagmaGeneAnalysisTask.create(
    asset_id="pgc_2022_sch_magma_entrez_gene_analysis",
    magma_annotation_task=PGC2022_SCH_MAGMA_ENTREZ_ANNOTATIONS,
    magma_p_value_task=SCH_PGC_2022_STANDARD_ANALYSIS.magma_tasks.inner.p_value_task,  # DECODE_ME_GWAS_1_BUILD_37_MAGMA_SNP_P_VALUES,
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    magma_ld_ref_task=MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
    ld_ref_file_stem="g1000_eur",
    sample_size=77258,  # source: GWAS file.  MAGMA documentation seems to imply we should use total sample size, not effective sample size
)
