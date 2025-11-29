"""
Tissue-specific LD scores for use with S-LDSC

Derived from the GTEx and Franke Lab gene expression datasets as described in Finucane, Hilary K., et al.
 "Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types." Nature genetics 50.4 (2018): 621-629.

 The Franke Lab dataset is described as:

 *an aggregation of publicly available microarray gene expression data sets comprising 37,427 samples in human, mouse, and rat*

See instructions here: https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses

Data is officially located at : https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

THOUSAND_GENOME_PARTITIONED_MODEL_MULTI_TISSUE_GENE_EXPR_LD_SCORES_RAW: DownloadFileTask = DownloadFileTask(
    meta=ReferenceFileMeta(
        asset_id="multi_tissue_gene_expression_thousand_genomes_partitioned_ld_scores",
        group="linkage_disequilibrium_scores",
        sub_group="LDSCORE_LDSC_SEG",
        sub_folder=PurePath("raw"),
        filename="LDSCORE_LDSC_SEG_ldscores_Multi_tissue_gene_expr_1000Gv3_ldscores",
        extension=".tar",
    ),
    url="https://www.dropbox.com/scl/fi/uabu3197gchjtuamn0691/LDSCORE_LDSC_SEG_ldscores_Multi_tissue_gene_expr_1000Gv3_ldscores.tar?rlkey=8lanxtrs146nsxrkhyqeiczmu&dl=1",
    md5_hash="fa2750675688474d8fceeab675a14f90",
)
