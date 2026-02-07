"""
Download a table containing Ensembl IDs, Entrez IDs, Gene symbols, and gene descriptions
for all known human genes.

Original source:
Ensembl Biomart (https://useast.ensembl.org/info/data/biomart/index.html)
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

GENE_THESAURUS = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="nomeclature",
        sub_group="ensembl_biomart",
        sub_folder=PurePath("raw"),
        filename="biomart_thesaurus",
        extension=".csv",
        id=AssetId("gene_thesaurus"),
        read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
    ),
    url="https://www.dropbox.com/scl/fi/9yamk0bze0f509p9opcnk/biomart_thesaurus.csv?rlkey=r9f4sedug2l3f6b5ak7hs88o2&dl=1",
    md5_hash="cd1774041d88a4f27c9b77e491c23263",
)
