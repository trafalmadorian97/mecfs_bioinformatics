import tempfile
from pathlib import Path
from typing import Iterator

import polars as pl
import pytest

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.reference.schemas.chrom_rename_rules import (
    CHROM_RENAME_RULES,
)
from mecfs_bio.build_system.reference.schemas.hg19_sn151_schema import (
    HG19_SNP151_SCHEMA,
)
from mecfs_bio.build_system.reference.schemas.hg19_snp151_schema_valid_choms import (
    HG19_SNP151_VALID_CHROMS,
)
from mecfs_bio.build_system.task.assign_rsids_via_snp151_task import (
    AssignRSIDSToSNPsViaSNP151Task,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.external_file_copy_task import ExternalFileCopyTask
from mecfs_bio.build_system.task.magma.magma_snp_location_task import MagmaSNPFileTask


@pytest.fixture(scope="package")
def parquet_snp151() -> Iterator[Path]:
    df = pl.scan_csv(
        "test_mecfs_bio/unit/build_system/task/dummy_snp151.txt",
        separator="\t",
        has_header=False,
        with_column_names=lambda x: HG19_SNP151_SCHEMA,
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        out_path = tmpdir_path / "out.parquet"
        df.sink_parquet(out_path)
        yield out_path


@pytest.fixture(scope="package")
def dummy_processed_gwas_parquet() -> Iterator[Path]:
    df = pl.scan_csv(
        "test_mecfs_bio/unit/build_system/task/dummy_processed_gwas.csv",
        separator=",",
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        out_path = tmpdir_path / "out.parquet"
        df.sink_parquet(out_path)
        yield out_path


@pytest.fixture(scope="package")
def raw_snp_data_task(dummy_processed_gwas_parquet: Path) -> Iterator[Task]:
    yield ExternalFileCopyTask(
        meta=SimpleFileMeta(
            AssetId("file_1"),
            read_spec=DataFrameReadSpec(format=DataFrameParquetFormat()),
        ),
        external_path=dummy_processed_gwas_parquet,
    )


@pytest.fixture(scope="package")
def snp_151_database_task(parquet_snp151: Path) -> Iterator[Task]:
    yield ExternalFileCopyTask(
        meta=SimpleFileMeta(
            AssetId("file_2"),
            read_spec=DataFrameReadSpec(format=DataFrameParquetFormat()),
        ),
        external_path=parquet_snp151,
    )


@pytest.fixture(scope="package")
def assign_rsids_task(
    raw_snp_data_task: Task, snp_151_database_task: Task
) -> Iterator[Task]:
    yield AssignRSIDSToSNPsViaSNP151Task.create(
        raw_snp_data_task=raw_snp_data_task,
        snp151_database_file_task=snp_151_database_task,
        asset_id="assign_rsids_via_snp_task",
        valid_chroms=HG19_SNP151_VALID_CHROMS,
        chrom_replace_rules=CHROM_RENAME_RULES,
    )


@pytest.fixture(scope="package")
def magma_snp_locations_task(assign_rsids_task: Task) -> Iterator[Task]:
    yield MagmaSNPFileTask.create_for_magma_snp_pos_file(
        gwas_parquet_with_rsids_task=assign_rsids_task,
        asset_id="magma_snp_locations_task",
    )
