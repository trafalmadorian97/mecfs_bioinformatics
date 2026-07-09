"""
Task to split the monolithic POPs feature collection into a directory of
column-chunk files, streaming in bounded memory.

The full PoPS.features table is a single ~21.5GB tab-separated file (~18k gene rows
by tens of thousands of feature columns). POPs' munge_feature_directory.py reads
each file in its input directory wholesale via pandas, so feeding it the monolithic
file needs far more RAM than a laptop or the CI runner has. This task rewrites the
single file as many smaller files, each carrying the ENSGID column plus a disjoint
slice of feature columns, so the munge step reads one manageable chunk at a time.

The split is a single streaming pass over the source: the header determines the
column groups, then every data row is split and its slices appended to the
corresponding output files. Peak memory is one input row plus the output write
buffers, independent of the size of the whole table.
"""

from pathlib import Path
from typing import TextIO

import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()

# First column of the POPs feature table, shared by every output chunk.
ENSGID_COL = "ENSGID"

# Output chunk filenames. munge_feature_directory.py globs every file in its input
# directory, so the directory must contain only these chunk files.
CHUNK_FILENAME_TEMPLATE = "features_chunk_{:04d}.txt"


@frozen
class PopsSplitFeatureFileTask(Task):
    """Split the monolithic POPs feature file into a directory of column chunks.

    columns_per_file bounds how many feature columns (excluding ENSGID) each output
    chunk holds, which in turn bounds the memory the downstream munge step needs to
    read one chunk.
    """

    meta: Meta
    source_features_task: Task
    columns_per_file: int = 2000

    @property
    def deps(self) -> list["Task"]:
        return [self.source_features_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self.source_features_task.asset_id)
        assert isinstance(source_asset, FileAsset)
        out_dir = scratch_dir / "feature_chunks"
        out_dir.mkdir(parents=True, exist_ok=True)
        num_chunks = split_feature_file(
            source_asset.path, out_dir, self.columns_per_file
        )
        logger.info(
            "Split POPs feature file into column chunks",
            num_chunks=num_chunks,
            columns_per_file=self.columns_per_file,
        )
        return DirectoryAsset(out_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        source_features_task: Task,
        columns_per_file: int = 2000,
    ) -> "PopsSplitFeatureFileTask":
        return cls(
            meta=SimpleDirectoryMeta(id=AssetId(asset_id)),
            source_features_task=source_features_task,
            columns_per_file=columns_per_file,
        )


def split_feature_file(source_path: Path, out_dir: Path, columns_per_file: int) -> int:
    """Stream a tab-separated feature table into column-chunk files under out_dir.

    Each chunk file has the ENSGID column followed by up to columns_per_file feature
    columns. Returns the number of chunk files written. Reads the source in a single
    pass, holding at most one row plus the open output handles in memory.
    """
    assert columns_per_file > 0, "columns_per_file must be positive"
    with source_path.open("r", encoding="utf-8") as source:
        header = source.readline().rstrip("\n").split("\t")
        assert header[0] == ENSGID_COL, (
            f"Expected first column {ENSGID_COL}, got {header[0]!r}"
        )
        num_columns = len(header)
        # Column index ranges (into the split row) for each chunk, skipping ENSGID.
        chunk_ranges = [
            (start, min(start + columns_per_file, num_columns))
            for start in range(1, num_columns, columns_per_file)
        ]
        out_files: list[TextIO] = []
        try:
            for chunk_index, (start, end) in enumerate(chunk_ranges):
                out_path = out_dir / CHUNK_FILENAME_TEMPLATE.format(chunk_index)
                out_file = out_path.open("w", encoding="utf-8")
                out_files.append(out_file)
                out_file.write(ENSGID_COL + "\t" + "\t".join(header[start:end]) + "\n")
            for line in source:
                fields = line.rstrip("\n").split("\t")
                assert len(fields) == num_columns, (
                    f"Row has {len(fields)} fields, expected {num_columns}"
                )
                ensgid = fields[0]
                for out_file, (start, end) in zip(out_files, chunk_ranges):
                    out_file.write(ensgid + "\t" + "\t".join(fields[start:end]) + "\n")
        finally:
            for out_file in out_files:
                out_file.close()
    return len(chunk_ranges)
