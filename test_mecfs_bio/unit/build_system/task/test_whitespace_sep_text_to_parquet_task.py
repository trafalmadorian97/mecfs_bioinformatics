import gzip
from pathlib import Path

import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
    DataFrameTextFormat,
    DataFrameWhiteSpaceSepTextFormat,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.whitespace_sep_text_to_parquet_task import (
    WhitespaceSepTextToParquetTask,
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF

# Variable-width spacing (like the deCODE alignment padding) plus a column that is "NA"
# in early rows and numeric later, to exercise streaming parse and cross-chunk schema.
_ROWS = [
    "Chr PosB38 Cohorts EAFrq I2",
    "chr1 100 -?????   0.10 NA",
    "chr1 200 ?+???-   0.20 NA",
    "chr2 300 -----?   0.30 NA",
    "chr2 400 ??-???   0.40 42",
    "chr3 500 ?----?   0.50 17",
]


def test_whitespace_sep_text_to_parquet_streams_and_preserves_values(tmp_path: Path):
    src = tmp_path / "gwas.txt.gz"
    with gzip.open(src, "wt") as f:
        f.write("\n".join(_ROWS) + "\n")
    scratch = tmp_path / "scratch"
    scratch.mkdir(parents=True, exist_ok=True)

    source_task = FakeTask(
        meta=GWASSummaryDataFileMeta(
            id=AssetId("raw"),
            trait="rheumatoid_arthritis",
            project="decode_ra_seropositive",
            sub_dir="raw",
            project_path=Path("gwas.txt.gz"),
            read_spec=DataFrameReadSpec(
                format=DataFrameWhiteSpaceSepTextFormat(comment_code="#")
            ),
        )
    )
    # chunk_size=2 forces three chunks, including the NA -> numeric transition.
    task = WhitespaceSepTextToParquetTask.create(
        source_task=source_task, asset_id="raw_parquet", chunk_size=2
    )

    class _Fetch(Fetch):
        def __call__(self, asset_id: AssetId) -> Asset:
            return FileAsset(src)

    result = task.execute(scratch_dir=scratch, fetch=_Fetch(), wf=SimpleWF())
    assert isinstance(result, FileAsset)

    df = pl.read_parquet(result.path)
    assert df.shape == (5, 5)
    assert df["Cohorts"].to_list() == ["-?????", "?+???-", "-----?", "??-???", "?----?"]
    assert df["PosB38"].to_list() == [100, 200, 300, 400, 500]
    # Column is float across all chunks despite starting all-NA, and NAs become null.
    assert df["I2"].dtype == pl.Float64
    assert df["I2"].to_list()[:3] == [None, None, None]
    assert df["I2"].to_list()[3:] == [42.0, 17.0]


def test_create_derives_trait_project_and_parquet_read_spec():
    source_task = FakeTask(
        meta=GWASSummaryDataFileMeta(
            id=AssetId("raw"),
            trait="rheumatoid_arthritis",
            project="decode_ra_seropositive",
            sub_dir="raw",
            project_path=Path("gwas.txt.gz"),
            read_spec=DataFrameReadSpec(
                format=DataFrameWhiteSpaceSepTextFormat(comment_code="#")
            ),
        )
    )
    task = WhitespaceSepTextToParquetTask.create(
        source_task=source_task, asset_id="raw_parquet"
    )
    meta = task.meta
    assert isinstance(meta, FilteredGWASDataMeta)
    assert meta.trait == "rheumatoid_arthritis"
    assert meta.project == "decode_ra_seropositive"
    assert meta.read_spec is not None
    assert isinstance(meta.read_spec.format, DataFrameParquetFormat)


def test_rejects_non_whitespace_source_at_construction():
    # Shift-left: a non-whitespace source must fail when the task is constructed,
    # not later at execute time.
    source_task = FakeTask(
        meta=GWASSummaryDataFileMeta(
            id=AssetId("raw"),
            trait="rheumatoid_arthritis",
            project="decode_ra_seropositive",
            sub_dir="raw",
            project_path=Path("gwas.tsv"),
            read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
        )
    )
    with pytest.raises(AssertionError):
        WhitespaceSepTextToParquetTask.create(
            source_task=source_task, asset_id="raw_parquet"
        )
