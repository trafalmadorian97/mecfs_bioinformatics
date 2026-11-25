"""Reference asset: the Eagle hg38 (GRCh38) genetic map and its parsed parquet.

Provides the per-position recombination rate (cM/Mb) for recombination tracks of
GRCh38 locus plots, such as the coloc locus stack plots.

The raw map (genetic_map_hg38_withX.txt.gz) is space-separated with columns
"chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)"; the header names carry
characters we would only rename away, so the download's read-spec renames the
four columns positionally. Parsing is a plain read-plus-sort, so it is expressed
as a PipeDataFrameTask (read via the download's read-spec, sort by CHR/POS, write
parquet) rather than a bespoke task.
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
from mecfs_bio.build_system.task.pipes.sort_pipe import SortPipe
from mecfs_bio.constants.genetic_map_constants import (
    GMAP_CM_COL,
    GMAP_POS_COL,
    GMAP_RATE_COL,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL

_OUTPUT_COLUMNS = [GWASLAB_CHROM_COL, GMAP_POS_COL, GMAP_RATE_COL, GMAP_CM_COL]

GENETIC_MAP_HG38_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="genetic_map",
        sub_group="hg38",
        sub_folder=PurePath("raw"),
        id=AssetId("genetic_map_hg38_eagle"),
        extension=".txt.gz",
        read_spec=DataFrameReadSpec(
            DataFrameTextFormat(
                separator=" ",
                has_header=True,
                column_names=_OUTPUT_COLUMNS,
                schema_overrides={
                    GWASLAB_CHROM_COL: pl.Int64(),
                    GMAP_POS_COL: pl.Int64(),
                    GMAP_RATE_COL: pl.Float64(),
                    GMAP_CM_COL: pl.Float64(),
                },
            )
        ),
    ),
    url="https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz",
    md5_hash="70f0be81e1ad1314ccc88e48091bcba1",
)

GENETIC_MAP_HG38 = PipeDataFrameTask.create(
    source_task=GENETIC_MAP_HG38_RAW,
    asset_id="genetic_map_hg38",
    out_format=ParquetOutFormat(),
    pipes=[SortPipe(by=[GWASLAB_CHROM_COL, GMAP_POS_COL])],
    backend="polars",
)
