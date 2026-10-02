from pathlib import Path

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.verifying_trace_info import (
    VerifyingTraceInfo,
)


def test_serialization(tmp_path: Path):
    target_path = tmp_path / "dump_loc.yml"
    info = VerifyingTraceInfo(
        trace_store={
            SimpleFileMeta(AssetId("file1")).asset_id: (
                "fakehash",
                [(SimpleFileMeta(AssetId("file2")).asset_id, "fakehash2")],
            )
        }
    )
    info.serialize(target_path)
    info_recovered = VerifyingTraceInfo.deserialize(target_path)
    assert info_recovered == info
