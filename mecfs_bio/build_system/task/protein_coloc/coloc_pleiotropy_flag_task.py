"""
Flags coloc pairs whose protein lead lies near a highly pleiotropic deCODE SomaScan SMP pQTL
(Eldjarn et al. 2023, ST36), for reading deCODE coloc results: a colocalization driven by a
variant that moves dozens of SomaScan measurements says little about the one protein.

ST36 gives the SMP lead of each pleiotropic pQTL, not every variant in LD with it, so the flag
is proximity, not identity: true when the pair's protein lead (on the locus chromosome) lies
within distance_bp of an ST36 row marked pleiotropic for SMP; false when none does; null when
the pair has no protein lead (untested pairs).
"""

from pathlib import Path, PurePath

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.task.protein_coloc.coloc_display_table_task import (
    coloc_result_dir,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.decode_somascan_constants import (
    ST36_CHROM_COL,
    ST36_PLEIOTROPIC_SMP_COL,
    ST36_SMP_POS_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    LOCUS_CHR_COL,
    PROTEIN_LEAD_POS_COL,
    RAW_RESULTS_FILENAME,
)

_KB = 1_000
_ROW_COL = "__row__"


def near_st36_pleiotropic_smp_col(distance_bp: int) -> str:
    """The flag column's name, which states its distance (e.g. NEAR_ST36_PLEIOTROPIC_SMP_250KB)."""
    assert distance_bp > 0 and distance_bp % _KB == 0, (
        f"distance_bp must be a positive multiple of 1 kb: {distance_bp}"
    )
    return f"NEAR_ST36_PLEIOTROPIC_SMP_{distance_bp // _KB}KB"


def flag_near_pleiotropic(
    results: pl.DataFrame, st36: pl.DataFrame, distance_bp: int
) -> pl.DataFrame:
    """results with the flag column appended (row order kept)."""
    smp = st36.filter(pl.col(ST36_PLEIOTROPIC_SMP_COL)).select(
        pl.col(ST36_CHROM_COL).cast(pl.Int64).alias(LOCUS_CHR_COL),
        pl.col(ST36_SMP_POS_COL).cast(pl.Int64).alias("__st36_pos__"),
    )
    keyed = results.with_row_index(_ROW_COL)
    near = (
        keyed.select(
            _ROW_COL,
            pl.col(LOCUS_CHR_COL).cast(pl.Int64),
            pl.col(PROTEIN_LEAD_POS_COL).cast(pl.Int64),
        )
        .join(smp, on=LOCUS_CHR_COL, how="inner")
        .filter(
            (pl.col(PROTEIN_LEAD_POS_COL) - pl.col("__st36_pos__")).abs() <= distance_bp
        )
        .select(_ROW_COL)
        .unique()
    )
    flag = (
        pl.when(pl.col(PROTEIN_LEAD_POS_COL).is_null())
        .then(None)
        .otherwise(pl.col(_ROW_COL).is_in(near[_ROW_COL].implode()))
    )
    return keyed.with_columns(
        flag.alias(near_st36_pleiotropic_smp_col(distance_bp))
    ).drop(_ROW_COL)


@frozen(slots=True)
class ColocPleiotropyFlagTask(GeneratingTask):
    """
    coloc_task: a TraitProteinColocTask over deCODE SomaScan SMP assays.
    pleiotropic_task: the typed ST36 table (see St36PleiotropicPqtlsPipe).
    distance_bp: how far a protein lead may lie from an ST36 SMP lead and still be flagged.
    """

    meta: ResultTableMeta
    coloc_task: Task
    pleiotropic_task: Task
    distance_bp: int

    @property
    def deps(self) -> list[Task]:
        return [self.coloc_task, self.pleiotropic_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        results = pl.read_parquet(
            coloc_result_dir(fetch, self.coloc_task) / RAW_RESULTS_FILENAME
        )
        st36 = (
            scan_dataframe_asset(
                fetch(self.pleiotropic_task.asset_id), self.pleiotropic_task.meta
            )
            .collect()
            .to_polars()
        )
        flagged = flag_near_pleiotropic(results, st36, self.distance_bp)
        out_path = scratch_dir / f"{self.meta.asset_id}.parquet"
        flagged.write_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        coloc_task: Task,
        pleiotropic_task: Task,
        distance_bp: int = 250_000,
    ) -> "ColocPleiotropyFlagTask":
        near_st36_pleiotropic_smp_col(distance_bp)
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
            pleiotropic_task=pleiotropic_task,
            distance_bp=distance_bp,
        )
