from pathlib import Path

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.metadata_to_path.simple_meta_to_path import (
    SimpleMetaToPath,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.simple_hasher import (
    SimpleHasher,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.verifying_trace_info import (
    VerifyingTraceInfo,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.verifying_trace_rebuilder_core import (
    VerifyingTraceRebuilder,
)
from mecfs_bio.build_system.scheduler.topological_scheduler import topological
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.tasks.simple_tasks import find_tasks
from mecfs_bio.build_system.wf.base_wf import make_wf


def test_magma_snp_locations_task(tmp_path: Path, magma_snp_locations_task: Task):
    tasks = find_tasks([magma_snp_locations_task])

    wf = make_wf()
    info: VerifyingTraceInfo = VerifyingTraceInfo.empty()

    asset_dir = tmp_path / "asset_dir"
    asset_dir.mkdir(exist_ok=True, parents=True)
    meta_to_path = SimpleMetaToPath(root=asset_dir)

    tracer = SimpleHasher.md5_hasher()
    rebuilder = VerifyingTraceRebuilder(tracer)

    targets = [magma_snp_locations_task.meta.asset_id]

    # Verify that all files are created in the correct location
    store, info = topological(
        rebuilder=rebuilder,
        tasks=tasks,
        targets=targets,
        wf=wf,
        info=info,
        meta_to_path=meta_to_path,
    )
    asset = store[AssetId("magma_snp_locations_task")]
    assert isinstance(asset, FileAsset)
    df = (
        scan_dataframe_asset(asset, magma_snp_locations_task.meta).collect().to_pandas()
    )
    assert len(df) == 1
    assert len(df.columns) == 3
