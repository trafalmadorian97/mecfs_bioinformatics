from pathlib import Path, PurePath

import pandas as pd
import structlog

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.simple_filtered_gwas_data_meta import (
    SimpleFilteredGWASDataMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.filter_snps_by_frequency import (
    FilterSNPsFrequencyTask,
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF

logger = structlog.get_logger()


def test_filter_snps_by_frequency_task(tmp_path: Path):
    df_1_loc = tmp_path / "df_1.csv"
    df_1 = pd.DataFrame(
        {
            "variant_id": [
                "rs10_rare",
                "rs11_common",
                "rs12_common",
                "rs13_common",
                "rs14_rare",
            ],
            "chromosome": [1, 1, 1, 1, 1],
            "base_pair_location": [1000, 1001, 1002, 1003, 1004],
            "effect_allele": ["T", "T", "A", "C", "T"],
            "other_allele": ["A", "C", "G", "G", "C"],
            "effect_allele_frequency": [0.001, 0.01, 0.02, 0.50, 0.995],
        }
    )
    df_1.to_csv(df_1_loc, index=False)
    scratch_loc = tmp_path / "scratch"
    scratch_loc.mkdir(exist_ok=True, parents=True)
    task = FilterSNPsFrequencyTask(
        meta=SimpleFilteredGWASDataMeta(
            trait="fake",
            project="fake",
            sub_dir=PurePath("abc/def"),
            id=AssetId("my_task"),
            read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator=",")),
        ),
        raw_gwas_task=FakeTask(
            SimpleFilteredGWASDataMeta(
                trait="fake",
                project="fake",
                sub_dir=PurePath("abc/def"),
                id=AssetId("Fake"),
                read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator=",")),
            ),
        ),
        allele_freq_col="effect_allele_frequency",
        freq_thresh=0.01,
    )

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(df_1_loc)

    result = task.execute(scratch_dir=scratch_loc, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
    assert result.path.is_file()

    result_df = pd.read_parquet(result.path)
    assert set(result_df["variant_id"]) == set(
        [
            "rs11_common",
            "rs12_common",
            "rs13_common",
        ]
    )
