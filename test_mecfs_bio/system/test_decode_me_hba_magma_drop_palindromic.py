from pathlib import Path

from mecfs_bio.asset_generator.hba_magma_asset_generator import (
    generate_human_brain_atlas_magma_tasks,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_hba_magma_analysis import (
    DECODE_ME_HBA_MAGMA_TASKS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.magma.magma_plot_brain_atlas_result_with_stepwise_labels import (
    HBAIndepPlotOptions,
)
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import PlotSettings
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL
from test_mecfs_bio.system.util import log_on_error

DECODE_ME_HBA_MAGMA_TASKS_DROP_PALINDROMIC = generate_human_brain_atlas_magma_tasks(
    base_name="decode_me_hba_magma_tasks",
    gwas_parquet_with_rsids_task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.join_task,
    sample_size=275488,  # this is the total sample size, which the MAGMA manual seems to imply is correct.  Could also try effective sample size.
    plot_settings=PlotSettings("plotly_white"),
    include_independent_cluster_plot=True,
    pipes=[RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL)],
    hba_indep_plot_options=HBAIndepPlotOptions(annotation_text_size=13),
)


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
            DECODE_ME_HBA_MAGMA_TASKS_DROP_PALINDROMIC.terminal_tasks()
            + [
                DECODE_ME_HBA_MAGMA_TASKS_DROP_PALINDROMIC.magma_independent_clusters_csv
            ],
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
        print("yo")
