from pathlib import Path

import polars as pl
import pyarrow.parquet as pq
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.decode_database.build_decode_variant_info_task import (
    BuildDecodeVariantInfoTask,
    DownloadDecodeVariantInfoTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.decode_somascan_constants import VI_CHR_COL, VI_POS_COL
from test_mecfs_bio.unit.build_system.task.decode_database.decode_fixtures import (
    FakeParallelRangeDownloader,
    FakeS3Access,
    fake_ref,
    synthetic_rows,
    synthetic_variant_info,
    write_synthetic_variant_info,
)

_REF = fake_ref("variant_info_iceland.txt.gz")


def _run(tmp_path: Path, access: FakeS3Access) -> Path:
    """Download the raw text, then build the parquet from it; returns the parquet path."""
    download = DownloadDecodeVariantInfoTask.create("vi_raw", _REF)
    build = BuildDecodeVariantInfoTask.create("vi", download)
    wf = make_wf(
        s3_access=access, parallel_range_downloader=FakeParallelRangeDownloader()
    )
    assets: dict[AssetId, Asset] = {}
    for task in (download, build):
        scratch = tmp_path / f"scratch_{task.asset_id}"
        scratch.mkdir()
        assets[task.asset_id] = task.execute(
            scratch, fetch=lambda asset_id: assets[asset_id], wf=wf
        )
    raw = assets[download.asset_id]
    assert isinstance(raw, FileAsset)
    assert raw.path.read_bytes() == access.objects[_REF].read_bytes()
    out = assets[build.asset_id]
    assert isinstance(out, FileAsset)
    return out.path


def test_variant_info_dedups_sorts_and_checks_etag(tmp_path: Path) -> None:
    rows = synthetic_rows(seed=1)
    source = tmp_path / "vi.txt.gz"
    write_synthetic_variant_info(source, rows)
    out = _run(tmp_path, FakeS3Access(objects={_REF: source}, etag_only=True))

    table = pl.read_parquet(out)
    assert table.height == synthetic_variant_info(rows).unique().height
    assert table.equals(table.sort(VI_CHR_COL, VI_POS_COL, maintain_order=True))
    metadata = pq.ParquetFile(out).metadata
    assert all(
        metadata.row_group(i).num_rows <= 65_536 for i in range(metadata.num_row_groups)
    )


def test_variant_info_etag_mismatch_fails(tmp_path: Path) -> None:
    source = tmp_path / "vi.txt.gz"
    write_synthetic_variant_info(source, synthetic_rows(seed=1))
    access = FakeS3Access(
        objects={_REF: source}, etag_only=True, reported_md5={_REF: "0" * 32}
    )
    with pytest.raises(AssertionError):
        _run(tmp_path, access)
