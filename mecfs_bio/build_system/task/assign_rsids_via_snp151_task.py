"""
Assign RSIDs to variants via joining a database file.  Only works for single-nucleotide variations.
"""

from pathlib import Path
from typing import Mapping

import ibis
import narwhals as nw
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameFormat,
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class AssignRSIDSToSNPsViaSNP151Task(Task):
    """
    Assigns RSIDS to the SNP genetic variants in a file of Gwas summary statistics
    Uses SNP151 database file
    Assumes the GWASLAB naming conventions are used in the summary statistics file
    Assumes that both input files are in parquet format

    Note that non-SNP variations (e.g. insertions or deletions) are excluded.
    This operates exclusively on SNPs
    """

    _meta: Meta
    snp151_database_file_task: Task
    raw_snp_data_task: Task
    valid_chroms: list[str]
    chrom_replace_rules: Mapping[str, int]

    def __attrs_post_init__(self):
        db_readspec = self.database_meta.read_spec()
        assert db_readspec is not None
        assert isinstance(db_readspec.format, DataFrameParquetFormat)
        snp_readspec = self.snp_data_meta.read_spec()
        assert snp_readspec is not None
        assert isinstance(snp_readspec.format, DataFrameParquetFormat)

    @property
    def database_meta(self) -> Meta:
        return self.snp151_database_file_task.meta

    @property
    def database_id(self) -> AssetId:
        return self.database_meta.asset_id

    @property
    def snp_data_meta(self) -> Meta:
        return self.raw_snp_data_task.meta

    @property
    def snp_data_id(self) -> AssetId:
        return self.snp_data_meta.asset_id

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.snp151_database_file_task, self.raw_snp_data_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        rename_frame = nw.from_native(
            ibis.memtable(
                {
                    "chrom": list(self.chrom_replace_rules.keys()),
                    GWASLAB_CHROM_COL: list(self.chrom_replace_rules.values()),
                }
            )
        )
        source_data_asset = fetch(self.snp_data_id)
        source_data_lf = scan_dataframe_asset(
            source_data_asset, meta=self.snp_data_meta, parquet_backend="ibis"
        )
        # filter for single nucleotide changes only
        source_data_lf = source_data_lf.filter(
            nw.col(GWASLAB_NON_EFFECT_ALLELE_COL).str.len_chars() == 1
        ).filter(nw.col(GWASLAB_EFFECT_ALLELE_COL).str.len_chars() == 1)
        database_asset = fetch(self.database_id)
        database_lf = scan_dataframe_asset(
            database_asset, meta=self.database_meta, parquet_backend="ibis"
        )
        processed_database_lf = database_lf.filter(nw.col("class") == "single").filter(
            nw.col("chrom").is_in(self.valid_chroms)
        )
        processed_database_lf = processed_database_lf.join(rename_frame, on="chrom")
        processed_database_lf = processed_database_lf.with_columns(
            (nw.col("chromStart_zero_based") + 1).alias(GWASLAB_POS_COL),
            (nw.col("name")).alias(GWASLAB_RSID_COL),
        ).select(
            GWASLAB_CHROM_COL,
            GWASLAB_POS_COL,
            GWASLAB_RSID_COL,
        )
        result: nw.LazyFrame = source_data_lf.join(
            processed_database_lf, on=[GWASLAB_POS_COL, GWASLAB_CHROM_COL]
        )
        out_path = scratch_dir / "snps_with_rsids.parquet"
        result.sink_parquet(
            out_path,
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        snp151_database_file_task: Task,
        raw_snp_data_task: Task,
        asset_id: str,
        valid_chroms: list[str],
        chrom_replace_rules: Mapping[str, int],
    ):
        source_meta = raw_snp_data_task.meta
        meta = create_new_meta(source_meta, asset_id=asset_id)
        return cls(
            meta=meta,
            snp151_database_file_task=snp151_database_file_task,
            raw_snp_data_task=raw_snp_data_task,
            valid_chroms=valid_chroms,
            chrom_replace_rules=chrom_replace_rules,
        )


def create_new_meta(
    source_meta: Meta,
    asset_id: str,
    format: DataFrameFormat = DataFrameParquetFormat(),
    extension=".parquet",
) -> Meta:
    meta: Meta
    if isinstance(source_meta, SimpleFileMeta):
        meta = SimpleFileMeta(
            short_id=AssetId(asset_id),
            read_spec=DataFrameReadSpec(format=format),
        )
    elif isinstance(source_meta, FilteredGWASDataMeta):
        meta = FilteredGWASDataMeta(
            short_id=AssetId(asset_id),
            project=source_meta.project,
            trait=source_meta.trait,
            sub_dir=source_meta.sub_dir,
            read_spec=DataFrameReadSpec(format=format),
            extension=extension,
        )
    else:
        raise ValueError("unknown source meta")
    return meta
