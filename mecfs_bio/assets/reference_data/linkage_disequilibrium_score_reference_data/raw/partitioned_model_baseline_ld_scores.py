"""
This is one of the reference datasets from the Stratified Linkage Disequilibrium Score Regression Paper.
It contains per-functional-category LD scores for variants across the genome.
Functional categories


See: Finucane, Hilary K., et al. "Partitioning heritability by functional annotation using genome-wide association summary statistics." Nature genetics 47.11 (2015): 1228-1235.

Excerpt:

Functional category annotations were largely derived from  ROADMAP and ENCODE

In detail (from the paper):
• Coding, 3’-UTR, 5’-UTR promoter, and intron annotations from the RefSeq gene model were
obtained from UCSC17 and post-processed by Gusev et al.13
• Digital genomic footprint and transcription factor binding site annotations were obtained from
ENCODE3 and post-processed by Gusev et al.13
• The combined chromHMM/Segway annotations for six cell lines were obtained from Hoffman et
al.20 CTCF, promoter flanking, transcribed, transcription start site, strong enhancer, and weak
enhancers are a union the six cell lines; repressed is an intersection over the six cell lines.
• DNase I hypersensitive sites (DHSs) are a combination of ENCODE and Roadmap data, postprocessed by Trynka et al.5 We combined the cell-type-specific annotations into two annotations
for inclusion in the full baseline model: a union of all cell types, and a union of only fetal cell types.
DHS and fetal DHS.
• H3K4me1, H3K4me, and H3K9ac were all obtained from Roadmap and post-processed by Trynka
et al.5 For each, we took a union over cell types for the full baseline model, and used the individual
cell types for our cell-type-specific analysis.
• One version of H3K27ac was obtained from Roadmap, post-processed,18 and a second version of
H3K27ac was obtained from the data of Hnisz et al.19 For each, we took a union over cell types
for the full baseline model, and used the individual cell types for our cell-type-specific analysis.
• Super enhancers were also obtained from Hnisz et al,19 and comprise a subset of the H3K27ac
annotation from that paper. We took a union over cell types for the full baseline model
• Regions conserved in mammals were obtained from Lindblad-Toh et al.,21 post-processed by Ward
and Kellis.31
• FANTOM5 enhancers were obtained from Andersson et al.22
• For each of these 24 categories, we added a 500bp window around the category as an additional
category to keep our heritability estimates from being inflated by heritability in flanking regions.
• For each of DHS, H3K4me1, H3K4me3, and H3K9ac, we added a 100bp window around the ChIPseq peak as an additional category.
• We added an additional category containing all SNPs


Since this dataset is hosted in a requester-pays Google Cloud bucket, I have mirrored it in my dropbox.



Data is officially located at : https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

THOUSAND_GENOME_BASE_MODEL_PARTITIONED_LD_SCORES = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="thousand_genomes_baseline_model_partitioned_ld_scores",
        group="linkage_disequilibrium_scores",
        sub_group="LDSCORE_1000G_Phase3_baseline",
        sub_folder=PurePath("raw"),
        filename="LDSCORE_1000G_Phase3_baseline_v1.2_ldscores",
        extension=".tar",
    ),
    url="https://www.dropbox.com/scl/fi/w0syfa29pnhd3heq3slmq/LDSCORE_1000G_Phase3_baseline_v1.2_ldscores.tar?rlkey=vizywdwx79s7d0xg0t5yh2man&dl=1",
    md5_hash="4a8faf1ec665d9cbfb67ea18139d98bf",
)
