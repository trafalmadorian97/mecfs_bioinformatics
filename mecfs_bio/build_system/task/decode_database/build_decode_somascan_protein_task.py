"""
The stored form of one deCODE SomaScan protein: a directory with

  slim.parquet: BETA, SE, N (float32) of the shared variant index's rows, in index order, null
      for an index row the file lacks; the number of such rows is in its metadata;
  seek_index.gzidx.zst: the GZIDX seek index of the remote file, for window reads;
  posmap.parquet: the position map (see one_pass) with the remote object's key, size and md5
      and the header length in its metadata, so a reader can detect a changed object;
  significant_rows.parquet: every row with p <= 5e-8 on any chromosome (not read yet).

The task downloads the file (about 970 MB) with parallel range requests, makes one rapidgzip
pass over it, writes the four files and deletes the download.
"""

import resource
import time
from pathlib import Path, PurePath

import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.task.dataframe_output import write_parquet_table
from mecfs_bio.build_system.task.decode_database.construct_decode_variant_index_task import (
    download_protein_file,
)
from mecfs_bio.build_system.task.decode_database.decode_somascan_protein import (
    DecodeSomascanProtein,
)
from mecfs_bio.build_system.task.decode_database.line_blocks import MIB
from mecfs_bio.build_system.task.decode_database.one_pass import OnePassResult, one_pass
from mecfs_bio.build_system.task.directory_member_dataframe import (
    DirectoryMemberDataFrame,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.decode_somascan_constants import (
    DECODE_BETA_COL,
    DECODE_INDEX_NAME_COL,
    DECODE_N_COL,
    DECODE_SE_COL,
    DECODE_SIGNIFICANT_ROWS_MIN_MLOG10P,
    DecodeSomascanGroup,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_POS_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

logger = structlog.get_logger()

DECODE_SLIM_FILENAME = "slim.parquet"
DECODE_SEEK_INDEX_FILENAME = "seek_index.gzidx.zst"
DECODE_POSMAP_FILENAME = "posmap.parquet"
DECODE_SIGNIFICANT_ROWS_FILENAME = "significant_rows.parquet"
DECODE_SLIM_MEMBER = DirectoryMemberDataFrame(
    PurePath(DECODE_SLIM_FILENAME), DataFrameReadSpec(DataFrameParquetFormat())
)

# Position-map metadata keys.
POSMAP_SOURCE_KEY = "source_key"
POSMAP_COMPRESSED_SIZE_KEY = "compressed_size"
POSMAP_MD5_KEY = "md5"
POSMAP_HEADER_END_OFFSET_KEY = "header_end_offset"
POSMAP_FORMAT_VERSION_KEY = "posmap_format_version"
GZIDX_FORMAT_VERSION_KEY = "gzidx_format_version"
POSMAP_FORMAT_VERSION = "1"
GZIDX_FORMAT_VERSION = "1"
# Slim metadata key: index rows the file lacks (null in the slim).
SLIM_INDEX_ROWS_MISSING_KEY = "index_rows_missing"
# Significant-rows metadata keys.
SIGNIFICANT_ROWS_SCANNED_KEY = "rows_scanned"
SIGNIFICANT_MIN_MLOG10P_KEY = "min_mlog10p"

SLIM_COLUMNS = {
    DECODE_BETA_COL: GWASLAB_BETA_COL,
    DECODE_SE_COL: GWASLAB_SE_COL,
    DECODE_N_COL: GWASLAB_SAMPLE_SIZE_COLUMN,
}


@frozen(slots=True)
class BuildDecodeSomascanProteinTask(GeneratingTask):
    """
    index_task: the group's ConstructDecodeVariantIndexTask.
    protein: the manifest row of the file to build from.
    parallelization, chunk_size: rapidgzip workers and compressed bytes per seek point.
    posmap_step: uncompressed bytes between position-map rows.
    n_connections: parallel range requests for the download.
    """

    meta: ResultDirectoryMeta
    index_task: Task
    protein: DecodeSomascanProtein
    parallelization: int = 8
    chunk_size: int = 16 * MIB
    posmap_step: int = 1 * MIB
    n_connections: int = 128

    @property
    def deps(self) -> list[Task]:
        return [self.index_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        timer = _StepTimer(self.protein.key)
        index_asset = fetch(self.index_task.asset_id)
        assert isinstance(index_asset, FileAsset)
        index = _read_index(index_asset, self.index_task)
        download = scratch_dir / "protein.txt.gz"
        download_protein_file(self.protein, download, wf, self.n_connections)
        timer.step("download")
        result = one_pass(
            download,
            index,
            parallelization=self.parallelization,
            chunk_size=self.chunk_size,
            posmap_step=self.posmap_step,
        )
        timer.step("one pass")
        logger.info(
            f"{self.protein.key}: {result.index_rows_missing} of {index.height} index rows "
            "missing from the file"
        )
        download.unlink()
        out_dir = scratch_dir / "out"
        out_dir.mkdir()
        _write_outputs(result, self.protein, out_dir)
        timer.step("write")
        return DirectoryAsset(out_dir)

    @classmethod
    def create(
        cls,
        index_task: Task,
        protein: DecodeSomascanProtein,
        asset_id: str,
        index_name: str,
        group: DecodeSomascanGroup,
        parallelization: int = 8,
        chunk_size: int = 16 * MIB,
        posmap_step: int = 1 * MIB,
    ) -> "BuildDecodeSomascanProteinTask":
        return cls(
            meta=ResultDirectoryMeta(
                id=asset_id,
                trait=f"decode_somascan_{group}",
                project=protein.gene_symbol,
                sub_dir=PurePath(f"aligned/{index_name}_index"),
            ),
            index_task=index_task,
            protein=protein,
            parallelization=parallelization,
            chunk_size=chunk_size,
            posmap_step=posmap_step,
        )


class _StepTimer:
    def __init__(self, label: str) -> None:
        self.label = label
        self.start = time.perf_counter()

    def step(self, name: str) -> None:
        now = time.perf_counter()
        peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        logger.info(
            f"{self.label}: {name} took {now - self.start:.1f}s; peak RSS {peak_mb:.0f} MB"
        )
        self.start = now


def _read_index(index_asset: FileAsset, index_task: Task) -> pl.DataFrame:
    return (
        scan_dataframe_asset(
            index_asset, meta=index_task.meta, parquet_backend="polars"
        )
        .to_native()
        .select(GWASLAB_POS_COL, DECODE_INDEX_NAME_COL)
        .collect()
    )


def _write_outputs(
    result: OnePassResult,
    protein: DecodeSomascanProtein,
    out_dir: Path,
) -> None:
    slim = result.slim.select(
        pl.col(source).cast(pl.Float32).alias(name)
        for source, name in SLIM_COLUMNS.items()
    )
    write_parquet_table(
        table=slim.to_arrow().replace_schema_metadata(
            {SLIM_INDEX_ROWS_MISSING_KEY: str(result.index_rows_missing)}
        ),
        out_path=out_dir / DECODE_SLIM_FILENAME,
        compression="zstd",
        compression_level=None,
        byte_stream_split_columns=list(SLIM_COLUMNS.values()),
    )
    (out_dir / DECODE_SEEK_INDEX_FILENAME).write_bytes(result.seek_index)
    _write_with_metadata(
        result.posmap,
        out_dir / DECODE_POSMAP_FILENAME,
        {
            POSMAP_SOURCE_KEY: protein.key,
            POSMAP_COMPRESSED_SIZE_KEY: str(protein.size),
            POSMAP_MD5_KEY: protein.md5,
            POSMAP_HEADER_END_OFFSET_KEY: str(result.header_end_offset),
            POSMAP_FORMAT_VERSION_KEY: POSMAP_FORMAT_VERSION,
            GZIDX_FORMAT_VERSION_KEY: GZIDX_FORMAT_VERSION,
        },
    )
    _write_with_metadata(
        result.significant_rows,
        out_dir / DECODE_SIGNIFICANT_ROWS_FILENAME,
        {
            SIGNIFICANT_ROWS_SCANNED_KEY: str(result.rows_scanned),
            SIGNIFICANT_MIN_MLOG10P_KEY: repr(DECODE_SIGNIFICANT_ROWS_MIN_MLOG10P),
        },
    )


def _write_with_metadata(
    frame: pl.DataFrame, path: Path, metadata: dict[str, str]
) -> None:
    write_parquet_table(
        table=frame.to_arrow().replace_schema_metadata(metadata),
        out_path=path,
        compression="zstd",
        compression_level=None,
        byte_stream_split_columns=[],
    )
