from pathlib import Path

import pandas as pd

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.counting_task import CountingTask
from mecfs_bio.build_system.task.external_file_copy_task import ExternalFileCopyTask
from mecfs_bio.build_system.task.filter_snps_task import FilterSNPsTask


def test_filer_snps(tmp_path: Path):
    """
    Test that we can correctly filter a dataframe of SNPs
    Also check that we can use the must_rebuild field of the runner to force tasks to be rebuilt
    """
    external_dir = tmp_path / "external"
    store_path = tmp_path / "store.yaml"
    asset_root = tmp_path / "assets"
    external_dir.mkdir(parents=True, exist_ok=True)
    asset_root.mkdir(parents=True, exist_ok=True)
    df_1 = pd.DataFrame({"ID": [1, 2, 3]})
    df_2 = pd.DataFrame({"ID": [1]})
    path_1 = external_dir / "df_1.parqet"
    path_2 = external_dir / "df_2.parqet"
    df_1.to_parquet(path_1)
    df_2.to_parquet(path_2)
    task_1 = CountingTask(
        ExternalFileCopyTask(
            meta=SimpleFileMeta(
                AssetId("file_1"), read_spec=DataFrameReadSpec(DataFrameParquetFormat())
            ),
            external_path=path_1,
        )
    )

    task_2 = CountingTask(
        ExternalFileCopyTask(
            meta=SimpleFileMeta(
                AssetId("file_2"), read_spec=DataFrameReadSpec(DataFrameParquetFormat())
            ),
            external_path=path_2,
        )
    )
    task_3 = CountingTask(
        FilterSNPsTask(
            raw_gwas_task=task_1,
            snp_list_task=task_2,
            meta=FilteredGWASDataMeta(
                AssetId("file3"), "trait", "proj", "raw", extension=".parquet"
            ),
        )
    )
    runner = SimpleRunner(info_store=store_path, asset_root=asset_root)
    runner.run(targets=[task_3])
    result_df = pd.read_parquet(runner.meta_to_path(task_3.meta))
    assert len(result_df) == 1
    assert task_3.run_count == 1
    assert task_2.run_count == 1
    assert task_1.run_count == 1
    runner.run(targets=[task_3])
    assert task_3.run_count == 1
    assert task_2.run_count == 1
    assert task_1.run_count == 1

    runner.run(targets=[task_3], must_rebuild_transitive=[task_2])
    assert task_3.run_count == 2
    assert task_2.run_count == 2
    assert task_1.run_count == 1
