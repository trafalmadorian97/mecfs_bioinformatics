"""Parse the Eagle hg19 genetic map into a sorted parquet keyed by (CHR, POS).

The raw map (genetic_map_hg19_withX.txt.gz) is space-separated with columns
"chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)". We keep the per-position
recombination rate (cM/Mb) directly - it is what the explainability plot's
recombination track draws - plus the cumulative genetic-map position (cM). The
polyfun LD panel is hg19/build 37, so this map matches its coordinates.
"""

from pathlib import Path

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL

GMAP_POS_COL = "POS"
GMAP_RATE_COL = "recomb_rate_cm_per_mb"
GMAP_CM_COL = "genetic_map_cm"

_OUTPUT_COLUMNS = [GWASLAB_CHROM_COL, GMAP_POS_COL, GMAP_RATE_COL, GMAP_CM_COL]


@frozen
class ParseHg19GeneticMapTask(Task):
    meta: Meta
    raw_map_task: Task

    @property
    def deps(self) -> list["Task"]:
        return [self.raw_map_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        raw_asset = fetch(self.raw_map_task.asset_id)
        assert isinstance(raw_asset, FileAsset)
        # Read the four space-separated columns positionally (the header names
        # contain characters we would only rename away anyway).
        frame = pl.read_csv(
            raw_asset.path,
            separator=" ",
            has_header=True,
            new_columns=_OUTPUT_COLUMNS,
            schema_overrides={
                GWASLAB_CHROM_COL: pl.Int64,
                GMAP_POS_COL: pl.Int64,
                GMAP_RATE_COL: pl.Float64,
                GMAP_CM_COL: pl.Float64,
            },
        ).sort([GWASLAB_CHROM_COL, GMAP_POS_COL])

        out_path = scratch_dir / str(self.meta.asset_id)
        frame.write_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(cls, asset_id: str, raw_map_task: Task) -> "ParseHg19GeneticMapTask":
        source_meta = raw_map_task.meta
        if not isinstance(source_meta, ReferenceFileMeta):
            raise ValueError(f"Unknown meta for raw genetic map task: {source_meta}")
        meta = ReferenceFileMeta(
            group=source_meta.group,
            sub_group=source_meta.sub_group,
            sub_folder=source_meta.sub_folder,
            id=AssetId(asset_id),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
        return cls(meta=meta, raw_map_task=raw_map_task)
