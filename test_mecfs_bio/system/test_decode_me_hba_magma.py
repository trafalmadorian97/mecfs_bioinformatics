from pathlib import Path

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_hba_magma_analysis import (
    DECODE_ME_HBA_MAGMA_TASKS,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.magma.magma_forward_stepwise_select_task import (
    RETAINED_CLUSTERS_COLUMN,
)
from test_mecfs_bio.system.util import log_on_error


def test_run_hba_magma(tmp_path: Path):
    """
    Test that we can run MAGMA via the human brain atlas on the DECODE ME data
    """

    info_store = tmp_path / "info_store.yaml"
    asset_root = tmp_path / "asset_store"
    asset_root.mkdir(parents=True, exist_ok=True)

    with log_on_error(info_store):
        test_runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=info_store,
            asset_root=asset_root,
        )
        result = test_runner.run(
            DECODE_ME_HBA_MAGMA_TASKS.terminal_tasks()
            + [DECODE_ME_HBA_MAGMA_TASKS.magma_independent_clusters_csv],
            incremental_save=True,
        )
        assert result is not None
        df = (
            scan_dataframe_asset(
                result[
                    DECODE_ME_HBA_MAGMA_TASKS.magma_independent_clusters_csv.asset_id
                ],
                meta=DECODE_ME_HBA_MAGMA_TASKS.magma_independent_clusters_csv.meta,
            )
            .collect()
            .to_polars()
        )

        assert set(df[RETAINED_CLUSTERS_COLUMN].to_list()) == {
            "Cluster234",
            "Cluster419",
            "Cluster136",
        }
