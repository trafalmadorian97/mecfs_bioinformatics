"""
GRCh38 locations of human genes (MAGMA's NCBI38 Entrez gene locations), extracted from the zip
and read with the gene-info column names that the locus gene tracks expect.

The .loc file is whitespace separated with no header: Entrez id, chromosome, start, end, strand,
gene symbol.
"""

from pathlib import PurePath

from mecfs_bio.assets.reference_data.magma_gene_locations.raw.magma_entrez_gene_location_reference_data_build_38 import (
    MAGMA_GENE_LOCATION_REFERENCE_DATA_BUILD_38_RAW,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameWhiteSpaceSepTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.extraction_one_file_from_zip_task import (
    ExtractFromZipTask,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
    GENE_INFO_STRAND_COL,
)

_FILE_TO_EXTRACT = "NCBI38.gene.loc"

MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED = ExtractFromZipTask(
    meta=ReferenceFileMeta(
        group="magma_reference_data",
        sub_group="gene_locations",
        sub_folder=PurePath("extracted"),
        id="magma_gene_location_reference_data_build_38_extracted",
        filename=_FILE_TO_EXTRACT,
        extension=".loc",
        read_spec=DataFrameReadSpec(
            DataFrameWhiteSpaceSepTextFormat(
                comment_code="#",
                col_names=[
                    "entrez_id",
                    GENE_INFO_CHROM_COL,
                    GENE_INFO_START_COL,
                    GENE_INFO_END_COL,
                    GENE_INFO_STRAND_COL,
                    GENE_INFO_NAME_COL,
                ],
            )
        ),
    ),
    source_file_task=MAGMA_GENE_LOCATION_REFERENCE_DATA_BUILD_38_RAW,
    file_to_extract=_FILE_TO_EXTRACT,
)
