"""Per-variant view of an external and an internal prior's SUSIE runs against the
uniform run at one locus, shared by the prior comparison table and plots.

Every locus variant gets each run's credible-set number and PIP, and each prior's
lift m*pi_i (from its contrast task's prior_lift.parquet). The three runs share
one variant set (a prior run requires its prior to cover every variant), so the
PIP and lift columns are filled on every row; a credible-set column is null only
when the variant is not in that run's credible sets.
"""

from pathlib import Path

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
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
CS_COLS = [CS_PF_EXT_COL, CS_PF_INT_COL, CS_U_COL]

_RUN_COLS = [
    *CS_COLS,
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
class PriorComparisonRuns:
    """The SUSIE runs being compared at one locus under one run config: the
    uniform baseline, and each prior's run with the contrast task explaining it
    against that baseline."""

    susie_uniform_task: Task
    external_prior_susie_task: Task
    external_prior_contrast_task: Task
    internal_prior_susie_task: Task
    internal_prior_contrast_task: Task

    @property
    def tasks(self) -> list[Task]:
        return [
            self.susie_uniform_task,
            self.external_prior_susie_task,
            self.external_prior_contrast_task,
            self.internal_prior_susie_task,
            self.internal_prior_contrast_task,
        ]


def secondary_pos_display_col(
    secondary_position: SecondaryPositionFromSnpid | None,
) -> str | None:
    if secondary_position is None:
        return None
    return f"pos_{secondary_position.build_label}"


def load_prior_comparison_variants(
    fetch: Fetch,
    runs: PriorComparisonRuns,
    secondary_position: SecondaryPositionFromSnpid | None,
) -> pl.DataFrame:
    """One row per locus variant with display columns chr, pos, (the secondary
    position, when configured), nea, ea, then each run's credible-set number, PIP,
    and each prior's lift."""
    uni_dir = _dir(fetch, runs.susie_uniform_task)
    ext_dir = _dir(fetch, runs.external_prior_susie_task)
    int_dir = _dir(fetch, runs.internal_prior_susie_task)
    rows = _pip(load_run_variants(uni_dir, secondary_position), PIP_U_COL)
    for frame in [
        _pip(load_run_variants(ext_dir), PIP_PF_EXT_COL),
        _pip(load_run_variants(int_dir), PIP_PF_INT_COL),
        _lift(fetch, runs.external_prior_contrast_task, LIFT_EXT_COL),
        _lift(fetch, runs.internal_prior_contrast_task, LIFT_INT_COL),
        _cs(uni_dir, CS_U_COL),
        _cs(ext_dir, CS_PF_EXT_COL),
        _cs(int_dir, CS_PF_INT_COL),
    ]:
        rows = rows.join(frame, on=VARIANT_KEY, how="left")
    n_missing = {
        col: rows[col].null_count()
        for col in _FILLED_ON_EVERY_ROW
        if rows[col].null_count()
    }
    assert not n_missing, (
        f"Runs do not share one variant set; variants missing a PIP or lift: "
        f"{n_missing}"
    )
    return _to_display_columns(rows, secondary_pos_display_col(secondary_position))


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


def _cs(run_dir: Path, cs_col: str) -> pl.DataFrame:
    return load_cs_numbers(run_dir).select(
        *VARIANT_KEY, pl.col(CS_NUMBER_COL).alias(cs_col)
    )


def _to_display_columns(
    rows: pl.DataFrame, secondary_pos_col: str | None
) -> pl.DataFrame:
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
    )
