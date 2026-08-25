import narwhals
from wcmatch.pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.concat_frames_task import ConcatFramesTask
from mecfs_bio.build_system.task.dataframe_output import ParquetOutFormat
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.expr_pipe import ExprPipe

POLYFUN_H_WEIGHT_COL = "snpvar_bin"
POLYFUN_PRIOR_COL = "prior"


def create_prior_col_pipe(q: int) -> DataProcessingPipe:
    max_weight = narwhals.col(POLYFUN_H_WEIGHT_COL).max()
    floor = max_weight / q
    prior_with_floor_pipe = ExprPipe(
        narwhals.when(narwhals.col(POLYFUN_H_WEIGHT_COL) <= floor)
        .then(floor)
        .otherwise(narwhals.col(POLYFUN_H_WEIGHT_COL))
        .alias(POLYFUN_PRIOR_COL)
    )
    return prior_with_floor_pipe


POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHT_CHR_1_7 = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="polyfun",
        sub_group="precomputed_prior",
        sub_folder=PurePath("raw"),
        id=AssetId("polyfun_precomputed_heritability_weight_chr_1_7"),
        extension=".parquet",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    ),
    url="https://github.com/omerwe/polyfun/raw/refs/heads/master/snpvar_meta.chr1_7.parquet",
    md5_hash="2f6a1509843edb954c4e1a200983683c",
)

POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHT_CHR_8_22 = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="polyfun",
        sub_group="precomputed_prior",
        sub_folder=PurePath("raw"),
        id=AssetId("polyfun_precomputed_heritability_weight_chr_8_22"),
        extension=".parquet",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    ),
    url="https://github.com/omerwe/polyfun/raw/refs/heads/master/snpvar_meta.chr8_22.parquet",
    md5_hash=None,
)

COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS = ConcatFramesTask.create(
    asset_id="polyfun_precomputed_heritability_weight_concat",
    frames_tasks=[
        POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHT_CHR_1_7,
        POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHT_CHR_8_22,
    ],
    out_format=ParquetOutFormat(),
)
