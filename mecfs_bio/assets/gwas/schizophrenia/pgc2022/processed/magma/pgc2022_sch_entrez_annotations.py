"""
Associate Entrez gene ids with SNPs from the schizophrenia GWAS
Note that we use an asymmetric window of 35 kb/10kb.
This is consistent with what is done in the Duncan paper:
https://github.com/Integrative-Mental-Health-Lab/linking_cell_types_to_brain_phenotypes/blob/675b5c9b58b8762934183a3ca61ae49ad587934a/MAGMA/1.annotationAndGeneAnalysis.sh#L12
"""

from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.standard_analysis_sc_pgc_2022 import (
    SCH_PGC_2022_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.processed.magma_entrez_gene_locations_data_build_37_unzipped import (
    MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_37_EXTRACTED,
)
from mecfs_bio.build_system.task.magma.magma_annotate_task import MagmaAnnotateTask

PGC2022_SCH_MAGMA_ENTREZ_ANNOTATIONS = MagmaAnnotateTask.create(
    asset_id="pgc2022_sch_magma_annotations",
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    snp_loc_file_task=SCH_PGC_2022_STANDARD_ANALYSIS.magma_tasks.inner.snp_loc_task,
    gene_loc_file_task=MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_37_EXTRACTED,
    window=(35, 10),
)
