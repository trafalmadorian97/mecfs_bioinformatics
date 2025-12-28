"""
Task to combine lists of candidate genes generated from MAGMA and GWASLAB
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_lead_variants_gene_labeled import (
    DECODE_ME_GWAS_1_GWASLAB_FILTERED_GENE_LIST_WITH_GENE_METADATA,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_filtered_gene_list import (
    DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST,
)
from mecfs_bio.build_system.task.combine_gene_lists_task import (
    CombineGeneListsTask,
    SrcGeneList,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat

DECODE_ME_GWAS_1_COMBINED_GENE_LISTS = CombineGeneListsTask.create(
    asset_id="decode_me_gwas_1_combined_gene_list",
    src_gene_lists=[
        SrcGeneList(
            task=DECODE_ME_GWAS_1_GWASLAB_FILTERED_GENE_LIST_WITH_GENE_METADATA,
            name="GWASLAB",
            ensemble_id_column="Gene stable ID",
        ),
        SrcGeneList(
            DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST,
            name="MAGMA",
            ensemble_id_column="GENE",
        ),
    ],
    out_format=ParquetOutFormat(),
)
