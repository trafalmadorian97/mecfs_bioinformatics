import gzip
from pathlib import Path, PurePath

import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.dataframe_output import ParquetOutFormat
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.pipe_dataframe_task import PipeDataFrameTask
from mecfs_bio.build_system.task.pipes.sort_pipe import SortPipe
from mecfs_bio.build_system.wf.base_wf import make_wf


def test_reads_gzip_text_with_positional_rename_and_overrides_then_pipes(
    tmp_path: Path,
) -> None:
    # A gzip-compressed, space-separated source whose header names we replace
    # positionally, with per-column dtype overrides -- read entirely from the
    # source meta's read-spec, then transformed by a pipe and written to parquet.
    # This exercises the text-read path of PipeDataFrameTask that the (deleted)
    # bespoke genetic-map parse task previously covered.
    raw = tmp_path / "source.txt.gz"
    lines = [
        "chr position rate(cM/Mb)",
        "2 100 1.5",
        "1 721290 2.685807669",
        "1 55550 0",
    ]
    raw.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))

    columns = ["CHR", "POS", "rate"]
    source_task = FakeTask(
        ReferenceFileMeta(
            group="genetic_map",
            sub_group="hg19",
            sub_folder=PurePath("raw"),
            id=AssetId("source_raw"),
            extension=".txt.gz",
            read_spec=DataFrameReadSpec(
                DataFrameTextFormat(
                    separator=" ",
                    has_header=True,
                    column_names=columns,
                    schema_overrides={
                        "CHR": pl.Int64(),
                        "POS": pl.Int64(),
                        "rate": pl.Float64(),
                    },
                )
            ),
        )
    )
    task = PipeDataFrameTask.create(
        source_task=source_task,
        asset_id="sorted_out",
        out_format=ParquetOutFormat(),
        pipes=[SortPipe(by=["CHR", "POS"])],
        backend="polars",
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "source_raw":
            return FileAsset(raw)
        raise ValueError("unknown asset id")

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    df = pl.read_parquet(result.path)

    # Header renamed positionally; overrides applied.
    assert df.columns == columns
    assert df.schema["CHR"] == pl.Int64
    assert df.schema["POS"] == pl.Int64
    assert df.schema["rate"] == pl.Float64
    # SortPipe ordered rows by (CHR, POS).
    assert df["CHR"].to_list() == [1, 1, 2]
    assert df["POS"].to_list() == [55550, 721290, 100]
    assert abs(df["rate"][1] - 2.685807669) < 1e-9
