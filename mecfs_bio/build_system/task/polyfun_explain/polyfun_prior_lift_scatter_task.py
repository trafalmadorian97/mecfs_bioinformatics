"""Interactive scatter of an external against an internal prior's lift over every
variant at one locus.

The lift m*pi_i is a variant's prior probability relative to uniform, normalized
to mean 1 over the locus, so both priors share a scale and the dashed y = x line
is where they agree. Both axes are log-scale because lift is multiplicative. A
least-squares line through the log lifts shows the overall relationship: a slope
below 1 means the internal prior spreads its lifts less than the external one.
The title reports that slope and the Spearman correlation of the lifts.

Points are colored by the larger of the two prior-run PIPs, on a scale from 0
to the locus maximum so a locus whose PIPs stay low still shows contrast, and
are drawn in increasing PIP order, so variants either prior run ranks highly sit
on top.
Hovering over a point names the variant and shows its PIPs, credible-set
numbers, and lifts.
"""

from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import polars as pl
from attrs import frozen
from scipy.stats import spearmanr

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.plot_file_meta import GWASPlotFileMeta
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    DISP_CHR,
    DISP_EA,
    DISP_NEA,
    DISP_POS,
    SecondaryPositionFromSnpid,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_prior_comparison_variants import (
    CS_PF_EXT_COL,
    CS_PF_INT_COL,
    CS_U_COL,
    LIFT_EXT_COL,
    LIFT_INT_COL,
    PIP_PF_EXT_COL,
    PIP_PF_INT_COL,
    PIP_U_COL,
    PriorComparisonRuns,
    load_prior_comparison_variants,
    secondary_pos_display_col,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.plotting.save_fig import PlotlyWriteMode

_MAX_PRIOR_PIP_COL = "max_prior_pip"
_VARIANT_LABEL_COL = "variant"


@frozen(slots=True)
class LiftFit:
    """Summary of how the internal prior's lifts track the external prior's."""

    spearman: float
    # Least-squares fit of log10 internal lift on log10 external lift.
    log_slope: float
    log_intercept: float


@frozen(slots=True)
class PolyfunPriorLiftScatterTask(Task):
    """Write the external-vs-internal prior lift scatter as a standalone html."""

    meta: Meta
    runs: PriorComparisonRuns
    secondary_position: SecondaryPositionFromSnpid | None = None
    plotly_js_mode: bool | PlotlyWriteMode = "cdn"

    @property
    def deps(self) -> list["Task"]:
        return self.runs.tasks

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        variants = _with_hover_columns(
            load_prior_comparison_variants(fetch, self.runs, self.secondary_position)
        )
        fit = _fit_lifts(variants)
        fig = _lift_scatter_figure(
            variants,
            fit=fit,
            secondary_pos_col=secondary_pos_display_col(self.secondary_position),
        )
        out_path = scratch_dir / "prior_lift_scatter.html"
        fig.write_html(out_path, include_plotlyjs=self.plotly_js_mode)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        runs: PriorComparisonRuns,
        secondary_position: SecondaryPositionFromSnpid | None = None,
    ) -> "PolyfunPriorLiftScatterTask":
        source_meta = runs.susie_uniform_task.meta
        if not isinstance(source_meta, ResultDirectoryMeta):
            raise ValueError(f"Unknown meta for uniform susie task: {source_meta}")
        meta = GWASPlotFileMeta(
            trait=source_meta.trait,
            project=source_meta.project,
            extension=".html",
            id=AssetId(asset_id),
        )
        return cls(meta=meta, runs=runs, secondary_position=secondary_position)


def _with_hover_columns(variants: pl.DataFrame) -> pl.DataFrame:
    """Add the variant label and the larger prior PIP, sorted by that PIP so the
    highest-PIP points are drawn last (on top)."""
    return variants.with_columns(
        pl.concat_str(
            [
                pl.col(DISP_CHR).cast(pl.String),
                pl.col(DISP_POS).cast(pl.String),
                pl.col(DISP_NEA),
                pl.col(DISP_EA),
            ],
            separator=":",
        ).alias(_VARIANT_LABEL_COL),
        pl.max_horizontal(PIP_PF_EXT_COL, PIP_PF_INT_COL).alias(_MAX_PRIOR_PIP_COL),
    ).sort(_MAX_PRIOR_PIP_COL)


def _fit_lifts(variants: pl.DataFrame) -> LiftFit:
    lift_ext = variants[LIFT_EXT_COL].to_numpy()
    lift_int = variants[LIFT_INT_COL].to_numpy()
    assert (lift_ext > 0).all() and (lift_int > 0).all(), (
        "Lifts must be positive to plot on log axes"
    )
    log_slope, log_intercept = np.polyfit(np.log10(lift_ext), np.log10(lift_int), 1)
    return LiftFit(
        spearman=float(spearmanr(lift_ext, lift_int).statistic),
        log_slope=float(log_slope),
        log_intercept=float(log_intercept),
    )


def _lift_scatter_figure(
    variants: pl.DataFrame, fit: LiftFit, secondary_pos_col: str | None
) -> go.Figure:
    lift_ext = variants[LIFT_EXT_COL].to_numpy()
    lift_int = variants[LIFT_INT_COL].to_numpy()
    lo = float(min(lift_ext.min(), lift_int.min()))
    hi = float(max(lift_ext.max(), lift_int.max()))
    diagonal = np.array([lo, hi])
    fit_x = np.array([float(lift_ext.min()), float(lift_ext.max())])
    fit_y = 10 ** (fit.log_intercept + fit.log_slope * np.log10(fit_x))

    hover_cols, hovertemplate = _hover(secondary_pos_col)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=diagonal,
            y=diagonal,
            mode="lines",
            line=dict(color="grey", dash="dash", width=1),
            name="y = x",
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=fit_x,
            y=fit_y,
            mode="lines",
            line=dict(color="firebrick", width=1.5),
            name=f"log-log fit (slope {fit.log_slope:.2f})",
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=lift_ext,
            y=lift_int,
            mode="markers",
            name="variants",
            marker=dict(
                size=7,
                color=variants[_MAX_PRIOR_PIP_COL].to_numpy(),
                colorscale="Viridis",
                cmin=0.0,
                colorbar=dict(title="max prior PIP"),
                line=dict(width=0.3, color="black"),
            ),
            customdata=variants.select(hover_cols)
            .with_columns(
                pl.col(c).cast(pl.String).fill_null("none")
                for c in (CS_PF_EXT_COL, CS_PF_INT_COL, CS_U_COL)
            )
            .to_numpy(),
            hovertemplate=hovertemplate,
        )
    )
    fig.update_xaxes(type="log", title="lift, external prior")
    fig.update_yaxes(type="log", title="lift, internal prior")
    fig.update_layout(
        template="plotly_white",
        title=(
            f"Prior lift, external vs internal "
            f"(Spearman {fit.spearman:.2f}, log-log slope {fit.log_slope:.2f}, "
            f"{variants.height} variants)"
        ),
        legend=dict(orientation="h", yanchor="bottom", y=1.0, x=0.0),
    )
    return fig


def _hover(secondary_pos_col: str | None) -> tuple[list[str], str]:
    """The customdata columns and the hover template reading them by index."""
    rows = [(_VARIANT_LABEL_COL, "<b>%{customdata[0]}</b>")]
    if secondary_pos_col is not None:
        rows.append((secondary_pos_col, f"{secondary_pos_col}: %{{customdata[1]}}"))
    labelled = [
        (PIP_PF_EXT_COL, "PIP external", ":.3f"),
        (PIP_PF_INT_COL, "PIP internal", ":.3f"),
        (PIP_U_COL, "PIP uniform", ":.3f"),
        (CS_PF_EXT_COL, "CS external", ""),
        (CS_PF_INT_COL, "CS internal", ""),
        (CS_U_COL, "CS uniform", ""),
        (LIFT_EXT_COL, "lift external", ":.3g"),
        (LIFT_INT_COL, "lift internal", ":.3g"),
    ]
    for col, text, fmt in labelled:
        rows.append((col, f"{text}: %{{customdata[{len(rows)}]{fmt}}}"))
    cols = [col for col, _ in rows]
    template = "<br>".join(line for _, line in rows) + "<extra></extra>"
    return cols, template
