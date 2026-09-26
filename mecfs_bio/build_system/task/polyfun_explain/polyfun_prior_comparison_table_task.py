"""Display table comparing an external and an internal prior's SUSIE runs
against the uniform run at one locus.

One row per variant in any of the three runs' credible sets, carrying each run's
credible-set number and PIP, and each prior's lift m*pi_i (from its contrast
task's prior_lift.parquet). The three runs share one variant set (a prior run
requires its prior to cover every variant), so the PIP and lift columns are
filled on every row; a credible-set column is null only when the variant is not
in that run's credible sets.
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
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    CS_NUMBER_COL,
    DISP_CHR,
    DISP_EA,
    DISP_LIFT,
    DISP_NEA,
    DISP_POS,
    PRIOR_LIFT_FILENAME,
    SECONDARY_POS_COL,
    VARIANT_KEY,
    SecondaryPositionFromSnpid,
    load_cs_numbers,
    load_run_variants,
)
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import PIP_COLUMN
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

CS_PF_EXT_COL = "cs_pf_ext"
CS_PF_INT_COL = "cs_pf_int"
CS_U_COL = "cs_u"
PIP_PF_EXT_COL = "pip_pf_ext"
PIP_PF_INT_COL = "pip_pf_int"
PIP_U_COL = "pip_u"
LIFT_EXT_COL = "lift_ext"
LIFT_INT_COL = "lift_int"

_RUN_COLS = [
    CS_PF_EXT_COL,
    CS_PF_INT_COL,
    CS_U_COL,
    PIP_PF_EXT_COL,
    PIP_PF_INT_COL,
    PIP_U_COL,
    LIFT_EXT_COL,
    LIFT_INT_COL,
]
_FILLED_ON_EVERY_ROW = [
    PIP_PF_EXT_COL,
    PIP_PF_INT_COL,
    PIP_U_COL,
    LIFT_EXT_COL,
    LIFT_INT_COL,
]


@frozen(slots=True)
class PolyfunPriorComparisonTableTask(Task):
    """Write the external-vs-internal prior comparison display table."""

    meta: Meta
    susie_uniform_task: Task
    external_prior_susie_task: Task
    external_prior_contrast_task: Task
    internal_prior_susie_task: Task
    internal_prior_contrast_task: Task
    secondary_position: SecondaryPositionFromSnpid | None = None

    @property
    def deps(self) -> list["Task"]:
        return [
            self.susie_uniform_task,
            self.external_prior_susie_task,
            self.external_prior_contrast_task,
            self.internal_prior_susie_task,
            self.internal_prior_contrast_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        uni_dir = _dir(fetch, self.susie_uniform_task)
        ext_dir = _dir(fetch, self.external_prior_susie_task)
        int_dir = _dir(fetch, self.internal_prior_susie_task)
        cs_frames = {
            CS_U_COL: load_cs_numbers(uni_dir),
            CS_PF_EXT_COL: load_cs_numbers(ext_dir),
            CS_PF_INT_COL: load_cs_numbers(int_dir),
        }
        per_variant_frames = [
            _pip(load_run_variants(uni_dir, self.secondary_position), PIP_U_COL),
            _pip(load_run_variants(ext_dir), PIP_PF_EXT_COL),
            _pip(load_run_variants(int_dir), PIP_PF_INT_COL),
            _lift(fetch, self.external_prior_contrast_task, LIFT_EXT_COL),
            _lift(fetch, self.internal_prior_contrast_task, LIFT_INT_COL),
        ]
        table = _comparison_table(
            cs_frames=cs_frames,
            per_variant_frames=per_variant_frames,
            secondary_pos_col=_secondary_pos_display_col(self.secondary_position),
        )
        out_path = scratch_dir / "prior_comparison_table.parquet"
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
        cls,
        asset_id: str,
        susie_uniform_task: Task,
        external_prior_susie_task: Task,
        external_prior_contrast_task: Task,
        internal_prior_susie_task: Task,
        internal_prior_contrast_task: Task,
        secondary_position: SecondaryPositionFromSnpid | None = None,
    ) -> "PolyfunPriorComparisonTableTask":
        source_meta = susie_uniform_task.meta
        if not isinstance(source_meta, ResultDirectoryMeta):
            raise ValueError(f"Unknown meta for uniform susie task: {source_meta}")
        meta = ResultTableMeta(
            id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            sub_dir=PurePath("analysis"),
        )
        return cls(
            meta=meta,
            susie_uniform_task=susie_uniform_task,
            external_prior_susie_task=external_prior_susie_task,
            external_prior_contrast_task=external_prior_contrast_task,
            internal_prior_susie_task=internal_prior_susie_task,
            internal_prior_contrast_task=internal_prior_contrast_task,
            secondary_position=secondary_position,
        )


def _dir(fetch: Fetch, task: Task) -> Path:
    asset = fetch(task.asset_id)
    assert isinstance(asset, DirectoryAsset)
    return asset.path


def _pip(variants: pl.DataFrame, pip_col: str) -> pl.DataFrame:
    return variants.rename({PIP_COLUMN: pip_col})


def _lift(fetch: Fetch, contrast_task: Task, lift_col: str) -> pl.DataFrame:
    # The lift inherits its prior's weight dtype (float32 for the precomputed
    # PolyFun prior); cast so both lift columns share one dtype.
    return pl.read_parquet(_dir(fetch, contrast_task) / PRIOR_LIFT_FILENAME).select(
        *VARIANT_KEY, pl.col(DISP_LIFT).cast(pl.Float64).alias(lift_col)
    )


def _secondary_pos_display_col(
    secondary_position: SecondaryPositionFromSnpid | None,
) -> str | None:
    if secondary_position is None:
        return None
    return f"pos_{secondary_position.build_label}"


def _comparison_table(
    cs_frames: dict[str, pl.DataFrame],
    per_variant_frames: list[pl.DataFrame],
    secondary_pos_col: str | None,
) -> pl.DataFrame:
    """Rows: the union of the runs' credible-set variants. Each cs_frames entry
    (display column -> load_cs_numbers frame) and each per-variant frame (keyed
    on VARIANT_KEY, covering every locus variant) is left-joined on. Sorted by
    the larger of the two prior PIPs, descending."""
    rows = pl.concat(
        [cs.select(VARIANT_KEY) for cs in cs_frames.values()], how="vertical"
    ).unique()
    for cs_col, cs in cs_frames.items():
        rows = rows.join(
            cs.select(*VARIANT_KEY, pl.col(CS_NUMBER_COL).alias(cs_col)),
            on=VARIANT_KEY,
            how="left",
        )
    for frame in per_variant_frames:
        rows = rows.join(frame, on=VARIANT_KEY, how="left")
    n_missing = {
        col: rows[col].null_count()
        for col in _FILLED_ON_EVERY_ROW
        if rows[col].null_count()
    }
    assert not n_missing, (
        f"Runs do not share one variant set; credible-set variants missing a "
        f"PIP or lift: {n_missing}"
    )
    rename_map = {
        GWASLAB_CHROM_COL: DISP_CHR,
        GWASLAB_POS_COL: DISP_POS,
        GWASLAB_NON_EFFECT_ALLELE_COL: DISP_NEA,
        GWASLAB_EFFECT_ALLELE_COL: DISP_EA,
    }
    position_cols = [DISP_CHR, DISP_POS]
    if secondary_pos_col is not None:
        rename_map[SECONDARY_POS_COL] = secondary_pos_col
        position_cols.append(secondary_pos_col)
    return (
        rows.rename(rename_map)
        .with_columns(pl.col(c).cast(pl.Int32) for c in position_cols)
        .select(*position_cols, DISP_NEA, DISP_EA, *_RUN_COLS)
        .sort(
            pl.max_horizontal(PIP_PF_EXT_COL, PIP_PF_INT_COL),
            descending=True,
        )
    )
