"""
MAGMA gene-location file derived from POPs' gene_annot_jun10.txt.

Used to run a POPs-specific MAGMA gene analysis whose gene set is a subset of POPs'
annotation, as POPs requires. See PopsMagmaGeneLocTask.
"""

from pathlib import PurePath

from mecfs_bio.assets.reference_data.pops.source.pops_source_extracted import (
    POPS_SOURCE_EXTRACTED,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.pops.pops_magma_gene_loc_task import (
    PopsMagmaGeneLocTask,
)

POPS_MAGMA_GENE_LOC = PopsMagmaGeneLocTask(
    meta=ReferenceFileMeta(
        group="pops",
        sub_group="magma_gene_loc",
        sub_folder=PurePath("derived"),
        id=AssetId("pops_magma_gene_loc"),
        extension=".gene.loc",
        filename="pops_gene_annot",
    ),
    pops_source_task=POPS_SOURCE_EXTRACTED,
)
