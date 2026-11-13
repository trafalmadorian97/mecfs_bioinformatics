"""
The HGNC complete set: one row per approved human gene, with its current symbol, name, locus
type, Ensembl gene id and the UniProt accessions of its protein products. Pinned to a dated
monthly snapshot.

Source: https://www.genenames.org/download/archive/
"""

from pathlib import PurePath

import polars as pl

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.constants.hgnc_constants import (
    HGNC_ENSEMBL_ID_COL,
    HGNC_LOCUS_TYPE_COL,
    HGNC_NAME_COL,
    HGNC_SYMBOL_COL,
    HGNC_UNIPROT_IDS_COL,
)

HGNC_COMPLETE_SET = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="nomeclature",
        sub_group="hgnc",
        sub_folder=PurePath("raw"),
        filename="hgnc_complete_set_2026-10-06",
        extension=".tsv",
        id=AssetId("hgnc_complete_set"),
        read_spec=DataFrameReadSpec(
            DataFrameTextFormat(
                separator="\t",
                schema_overrides={
                    col: pl.String()
                    for col in (
                        HGNC_SYMBOL_COL,
                        HGNC_NAME_COL,
                        HGNC_LOCUS_TYPE_COL,
                        HGNC_ENSEMBL_ID_COL,
                        HGNC_UNIPROT_IDS_COL,
                    )
                },
            )
        ),
    ),
    url="https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2026-10-06.tsv",
    md5_hash="af43fd56a64ad57c0fcc72d6ef65e516",
)
