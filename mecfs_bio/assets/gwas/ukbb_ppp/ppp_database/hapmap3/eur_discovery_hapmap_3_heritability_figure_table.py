"""
The per-protein HapMap3 EUR-discovery heritability table, encoded for display on
the documentation site by the data_table macro.

This is the same table as HAPMAP_3_PPP_HERITABILITY, re-encoded rather than
re-computed. At ~5,900 rows it is far too large for the markdown_table macro,
which renders every row into the page at build time; data_table instead ships
this parquet to the browser and renders it client-side.

Two encoding choices are load-bearing:

- n_snps is cast to Int32. Parquet INT64 decodes to JavaScript BigInt, which
  Tabulator cannot format --- it throws and renders an empty table while the
  data itself decodes fine, so the failure is easy to misread. The counts here
  are around 1.2e6, far inside Int32.
- Floats stay float64 and are not rounded, so the table's CSV download hands
  readers exactly the values the build system computed rather than a rounded
  approximation. The macro rounds for display only. BYTE_STREAM_SPLIT plus
  zstd is what keeps that affordable: it takes the file from ~830KB as plain
  CSV to ~236KB.
"""

import narwhals

from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.eur_discovery_hapmap3_ppp_heritability import (
    HAPMAP_3_PPP_HERITABILITY,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ParquetOutFormat,
    ParquetWriteOptions,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.cast_pipe import CastPipe

HAPMAP_3_PPP_HERITABILITY_FIGURE_TABLE = PipeDataFrameTask.create(
    asset_id="ppp_heritability_hapmap_3_eur_discovery_table",
    source_task=HAPMAP_3_PPP_HERITABILITY,
    pipes=[
        CastPipe(
            target_column="n_snps",
            type=narwhals.Int32(),
            new_col_name="n_snps",
        )
    ],
    out_format=ParquetOutFormat(
        write_options=ParquetWriteOptions(
            compression="zstd",
            # Level 22 costs a few seconds on a table this size and is paid once
            # at build time, against a file every reader of the page downloads.
            compression_level=22,
            byte_stream_split_floats=True,
        )
    ),
    backend="polars",
)
