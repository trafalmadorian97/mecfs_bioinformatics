from pathlib import PurePath

from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.magma.decode_me_gwas_1_build_37_magma_entrez_gene_analysis import (
    DECODE_ME_GWAS_1_MAGMA_ENTREZ_GENE_ANALYSIS,
)
from mecfs_bio.assets.reference_data.gene_set_data.for_magma.from_gsea_msigdb.extracted.gsea_entrez_human_gene_set_2025_v1_extracted import (
    GSEA_HUMAN_GENE_SET_2025_V1_EXTRACTED,
)
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import (
    DirectoryGeneSetSpec,
    MagmaGeneSetAnalysisTask,
)

MAGMA_DECODE_ME_GSEA_HALLMARK_GENE_SET_ANALYSIS = MagmaGeneSetAnalysisTask.create(
    asset_id="decode_me_gwas_1_build_37_magma_entrez_gse_hallmark_gene_set_analysis",
    magma_gene_analysis_task=DECODE_ME_GWAS_1_MAGMA_ENTREZ_GENE_ANALYSIS,
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    gene_set_task=DirectoryGeneSetSpec(
        gene_set_task=GSEA_HUMAN_GENE_SET_2025_V1_EXTRACTED,
        path_in_dir=PurePath(
            "msigdb_v2025.1.Hs_files_to_download_locally/msigdb_v2025.1.Hs_GMTs/h.all.v2025.1.Hs.entrez.gmt"
        ),
    ),
    set_or_covar="set",
    model_params=None,
)
