from pathlib import PurePath

from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec, DataFrameTextFormat, \
    DataFrameWhiteSpaceSepTextFormat
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.build_system.task.susie_stacked_plot_task import GENE_INFO_CHROM_COL, GENE_INFO_START_COL, \
    GENE_INFO_END_COL, GENE_INFO_STRAND_COL, GENE_INFO_NAME_COL

"""
Official site:https://fuma.ctglab.nl/downloadPage
File labeled "MAGMA gene boundaries Ensembl v110"
No direct link from official site, so hosted on dropbox.

"""
MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="magma_reference_data",
        sub_group="ensembl_gene_locations",
        sub_folder=PurePath("raw"),
        id="magma_ensembl_gene_location_reference_data_build_37_raw",
        extension=".txt",
        filename="ENSGv102.coding.genes",
        read_spec=DataFrameReadSpec(

            DataFrameWhiteSpaceSepTextFormat(
                col_names=["ensembl_name",
                              GENE_INFO_CHROM_COL, GENE_INFO_START_COL,
                              GENE_INFO_END_COL,
                              GENE_INFO_STRAND_COL,
                            GENE_INFO_NAME_COL
                           ],
                comment_code="#"
            )
        )
    ),
    url="https://www.dropbox.com/scl/fi/watz1v757xzjrawkl60se/ENSGv102.coding.genes.txt?rlkey=8lnk7i3nqxfnz7og8fasnn0zg&dl=1",
    md5_hash="8541ce5c7c95d67956672e1d37479ead",
)
