"""
The HGNC complete set: one row per approved human gene, with its current symbol, name, locus
type, Ensembl gene id and the UniProt accessions of its protein products. Pinned to a dated
monthly snapshot.

HGNC_UNIPROT_GENES normalizes it into a UniProt-to-gene table: one row per (UniProt accession,
gene), with the gene's symbol, name and Ensembl gene id. Genes without a UniProt accession are
dropped, and so are readthrough genes (a transcript spanning two adjacent genes, such as
PTGES3L-AARSD1, which HGNC lists under the accession of one of the genes it spans). An accession
encoded by duplicated genes (for example P59665, by DEFA1 and DEFA1B) has one row per gene.

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
from mecfs_bio.build_system.task.dataframe_output import ParquetOutFormat
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.build_system.task.pipe_dataframe_task import PipeDataFrameTask
from mecfs_bio.build_system.task.pipes.drop_null_pipe import DropNullsPipe
from mecfs_bio.build_system.task.pipes.drop_rows_by_value import DropRowsByValue
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe
from mecfs_bio.build_system.task.pipes.str_split_col import SplitColPipe
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe
from mecfs_bio.build_system.task.pipes.unnest_pipe import UNNestPipe
from mecfs_bio.constants.hgnc_constants import (
    HGNC_ENSEMBL_ID_COL,
    HGNC_LIST_SEPARATOR,
    HGNC_LOCUS_TYPE_COL,
    HGNC_NAME_COL,
    HGNC_READTHROUGH_LOCUS_TYPE,
    HGNC_SYMBOL_COL,
    HGNC_UNIPROT_ID_COL,
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

HGNC_UNIPROT_GENES = PipeDataFrameTask.create(
    source_task=HGNC_COMPLETE_SET,
    asset_id="hgnc_uniprot_genes",
    out_format=ParquetOutFormat(),
    backend="polars",
    pipes=[
        DropRowsByValue(
            target_column=HGNC_LOCUS_TYPE_COL,
            values_to_drop=[HGNC_READTHROUGH_LOCUS_TYPE],
        ),
        SelectColPipe(
            cols_to_select=[
                HGNC_SYMBOL_COL,
                HGNC_NAME_COL,
                HGNC_ENSEMBL_ID_COL,
                HGNC_UNIPROT_IDS_COL,
            ]
        ),
        SplitColPipe(
            col_to_split=HGNC_UNIPROT_IDS_COL,
            split_by=HGNC_LIST_SEPARATOR,
            new_col_name=HGNC_UNIPROT_ID_COL,
        ),
        UNNestPipe(col_to_unnest=HGNC_UNIPROT_ID_COL),
        SelectColPipe(
            cols_to_select=[
                HGNC_UNIPROT_ID_COL,
                HGNC_SYMBOL_COL,
                HGNC_NAME_COL,
                HGNC_ENSEMBL_ID_COL,
            ]
        ),
        DropNullsPipe(subset=[HGNC_UNIPROT_ID_COL]),
        UniquePipe(by=[HGNC_UNIPROT_ID_COL, HGNC_SYMBOL_COL], keep="any"),
    ],
)
