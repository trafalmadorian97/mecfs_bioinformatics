"""Tests for per-variant sample size support in the MAGMA gene analysis path.

A scalar sample size is passed to MAGMA as ``N=<value>``; a per-variant sample
size is written as an extra column in the (headerless) p-value file and passed as
``ncol=<index>``.
"""

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
from mecfs_bio.build_system.task.magma.magma_snp_location_task import MagmaSNPFileTask
from mecfs_bio.build_system.tasks.simple_tasks import find_tasks
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_P_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
)


def _p_value_task_with_n(assign_rsids_task: Task) -> MagmaSNPFileTask:
    return MagmaSNPFileTask.create_for_magma_snp_p_value_file_compute_if_needed(
        gwas_parquet_with_rsids_task=assign_rsids_task,
        asset_id="magma_pval_with_n",
        sample_size_column=GWASLAB_SAMPLE_SIZE_COLUMN,
    )


def test_pval_file_includes_n_column(tmp_path: Path, assign_rsids_task: Task):
    """When a per-variant sample size column is requested, the p-value file
    carries it as a third column (rsID, P, N)."""
    p_value_task = _p_value_task_with_n(assign_rsids_task)
    tasks = find_tasks([p_value_task])

    asset_dir = tmp_path / "asset_dir"
    asset_dir.mkdir(exist_ok=True, parents=True)
    rebuilder = VerifyingTraceRebuilder(SimpleHasher.md5_hasher())
    store, _info = topological(
        rebuilder=rebuilder,
        tasks=tasks,
        targets=[p_value_task.meta.asset_id],
        wf=make_wf(),
        info=VerifyingTraceInfo.empty(),
        meta_to_path=SimpleMetaToPath(root=asset_dir),
    )

    asset = store[AssetId("magma_pval_with_n")]
    assert isinstance(asset, FileAsset)
    df = scan_dataframe_asset(asset, p_value_task.meta).collect().to_pandas()
    assert list(df.columns) == [
        GWASLAB_RSID_COL,
        GWASLAB_P_COL,
        GWASLAB_SAMPLE_SIZE_COLUMN,
    ]
