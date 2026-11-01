"""
Docs-facing table of the colocalizations of one trait with proteins.

Keeps the tested (locus, protein) pairs of a TraitProteinColocTask result whose PP.H4 exceeds
h4_threshold, sorted by PP.H4 descending. Columns: GENE, PROTEIN_ID, LOCUS (the locus id,
chr:start-end in GRCh38), LEAD_VARIANT (the trait lead as chr:pos:NEA:EA), CIS_STATUS, PP_H4.
Written as zstd parquet with byte-stream-split floats and Int32 integers, for the data_table
docs macro. An empty table is a valid output.
"""

from pathlib import Path, PurePath

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import (
    float_column_names,
    write_parquet_table,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.protein_coloc_constants import (
    CIS_STATUS_COL,
    DISPLAY_LEAD_VARIANT_COL,
    DISPLAY_LOCUS_COL,
    GENE_COL,
    LOCUS_CHR_COL,
    LOCUS_ID_COL,
    LOCUS_LEAD_EA_COL,
    LOCUS_LEAD_NEA_COL,
    LOCUS_LEAD_POS_COL,
    PP_H4_COL,
    PROTEIN_ID_COL,
    RAW_RESULTS_FILENAME,
    STATUS_COL,
    STATUS_TESTED,
)

DISPLAY_SCHEMA = pl.Schema(
    {
        GENE_COL: pl.String,
        PROTEIN_ID_COL: pl.String,
        DISPLAY_LOCUS_COL: pl.String,
        DISPLAY_LEAD_VARIANT_COL: pl.String,
        CIS_STATUS_COL: pl.String,
        PP_H4_COL: pl.Float64,
    }
)


def read_raw_results(fetch: Fetch, coloc_task: Task) -> pl.DataFrame:
    asset = fetch(coloc_task.asset_id)
    assert isinstance(asset, DirectoryAsset), (
        f"expected {coloc_task.asset_id} to be a DirectoryAsset"
    )
    return pl.read_parquet(asset.path / RAW_RESULTS_FILENAME)


def display_rows(raw: pl.DataFrame, h4_threshold: float) -> pl.DataFrame:
    lead_variant = pl.format(
        "chr{}:{}:{}:{}",
        LOCUS_CHR_COL,
        LOCUS_LEAD_POS_COL,
        LOCUS_LEAD_NEA_COL,
        LOCUS_LEAD_EA_COL,
    )
    return (
        raw.filter(
            pl.col(STATUS_COL) == STATUS_TESTED, pl.col(PP_H4_COL) > h4_threshold
        )
        .sort([PP_H4_COL, PROTEIN_ID_COL], descending=[True, False])
        .select(
            GENE_COL,
            PROTEIN_ID_COL,
            pl.col(LOCUS_ID_COL).alias(DISPLAY_LOCUS_COL),
            lead_variant.alias(DISPLAY_LEAD_VARIANT_COL),
            CIS_STATUS_COL,
            PP_H4_COL,
        )
        .cast(DISPLAY_SCHEMA)
    )


@frozen(slots=True)
class ColocDisplayTableTask(Task):
    meta: Meta
    coloc_task: Task
    h4_threshold: float

    @property
    def deps(self) -> list[Task]:
        return [self.coloc_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        table = display_rows(
            read_raw_results(fetch, self.coloc_task), self.h4_threshold
        )
        out_path = scratch_dir / f"{self.meta.asset_id}.parquet"
        arrow_table = table.to_arrow()
        write_parquet_table(
            table=arrow_table,
            out_path=out_path,
            compression="zstd",
            compression_level=None,
            byte_stream_split_columns=float_column_names(arrow_table),
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls, asset_id: str, coloc_task: Task, h4_threshold: float = 0.8
    ) -> "ColocDisplayTableTask":
        """coloc_task: a TraitProteinColocTask (any Task producing its result directory)."""
        assert 0 < h4_threshold < 1, f"h4_threshold must be in (0, 1): {h4_threshold}"
        source_meta = coloc_task.meta
        assert isinstance(source_meta, ResultDirectoryMeta), (
            f"coloc_task must carry ResultDirectoryMeta, got {type(source_meta).__name__}"
        )
        return cls(
            meta=ResultTableMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                sub_dir=PurePath("analysis"),
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
            coloc_task=coloc_task,
            h4_threshold=h4_threshold,
        )
