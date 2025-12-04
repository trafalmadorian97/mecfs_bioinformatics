"""
Tissue-specific LD scores for use with S-LDSC

Derived from chromatin datasets as described in
 "Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types." Nature genetics 50.4 (2018): 621-629.

 The paper says:
*We downloaded narrow peaks from the Roadmap Epigenomics consortium for DNase
hypersensitivity and five activating histone marks: H3K27ac, H3K4me3, H3K4me1,
H3K9ac, and H3K36me3 (see URLs). Each of these six features was present in a subset of
the 88 primary cell types/tissues, for a total of 397 cell-type-/tissue-specific annotations. We
also analyzed peaks called using Homer from EN-TEx, a subgroup of the ENCODE project,
for four activating histone marks: H3K27ac, H3K4m3, H3K4me1, and H3K36me3. Each of
these four features was present in a subset of 27 tissues that were also included in the GTEx
data set, for a total of 93 cell-type-/tissue-specific annotations. For each of these two
datasets, of each of the annotations, we tested for enrichment by adding the annotation to the
baseline model (see Table S1), together with the union of cell-type-specific annotations
within each mark and the average of cell-type-specific annotations within each mark. A
positive regression coefficient for a tissue-/cell-type-specific annotation represents a positive
contribution of the annotation to per-SNP heritability, conditional on the other annotations.
We again computed a P-value to test whether the regression coefficient was positive*

Roughly speaking the underlying logic of this dataset is as follows: In certain tissue or cell types, certain regions of the genome are known to be marked by special histone marks in certain
. These marks indicate that a nearby gene is being transcribed in those tissues or cell types.


Data is officially located at : https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays

I re-hosted on dropbox.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

PARTITIONED_MODEL_MULTI_TISSUE_CHROMATIN_LD_SCORES_RAW: DownloadFileTask = DownloadFileTask(
    meta=ReferenceFileMeta(
        asset_id="multi_tissue_chromatin_partitioned_ld_scores",
        group="linkage_disequilibrium_scores",
        sub_group="LDSCORE_LDSC_SEG",
        sub_folder=PurePath("raw"),
        filename="LDSCORE_LDSC_SEG_ldscores_Multi_tissue_chromatin_1000Gv3_ldscores",
        extension=".tar",
    ),
    url="https://www.dropbox.com/scl/fi/8gdia1hhbjg65l4lvvpnt/LDSCORE_LDSC_SEG_ldscores_Multi_tissue_chromatin_1000Gv3_ldscores.tar?rlkey=94v5r1mjawbd1ha2gzmrg3bur&dl=1",
    md5_hash="8aa379558c9ff8986965b4337b5ce6d2",
)
