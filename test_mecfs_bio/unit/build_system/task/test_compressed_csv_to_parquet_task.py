import gzip
import shutil
from pathlib import Path, PurePath

import pandas as pd

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
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.reference.schemas.hg19_sn151_schema import (
    HG19_SNP151_SCHEMA,
)
from mecfs_bio.build_system.task.compressed_csv_to_parquet_task import (
    CompressedCSVToParquetTask,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def _gzip_file(src: Path, target: Path):
    with open(src, "rb") as f:
        with gzip.open(target, "wb") as g:
            shutil.copyfileobj(f, g)


def test_compressed_csv_to_parquet_task(tmp_path: Path):
    gzip_loc = tmp_path / "gzip_loc.tsv.gzip"
    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir(exist_ok=True, parents=True)
    _gzip_file(
        Path("test_mecfs_bio/unit/build_system/task/large_dummy_snp151.txt"), gzip_loc
    )
    dummy_task = FakeTask(
        meta=ReferenceFileMeta(
            group="",
            sub_folder=PurePath(""),
            sub_group="",
            id=AssetId("123"),
            read_spec=DataFrameReadSpec(
                format=DataFrameTextFormat(
                    separator="\t", column_names=HG19_SNP151_SCHEMA
                )
            ),
        )
    )
    task = CompressedCSVToParquetTask.create(
        csv_task=dummy_task,
        asset_id="parquet_Task",
    )

    class MyFetch(Fetch):
        def __call__(self, asset_id: AssetId) -> Asset:
            return FileAsset(gzip_loc)

    result = task.execute(scratch_dir=scratch_dir, fetch=MyFetch(), wf=SimpleWF())
    assert isinstance(result, FileAsset)
    pd.read_parquet(result.path)
