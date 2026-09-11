"""
Task that draws the LD-score-regression diagnostic plot for a heritability analysis.

LD score regression fits E[chi^2_j] = intercept + (N * h2 / M) * ld_j but, because there are
millions of variants, the fit is never shown -- standard regression diagnostics are essentially
never performed (see the discussion in gwaslab_snp_heritability_by_ldsc_task). This task supplies
one: it bins the analysis's variants into equal-count LD-score bins and plots each bin's mean
chi-square against the fitted line, so model misfit (bins drifting off the line) and confounding
(an intercept away from one, or a low-LD bin sitting below one) are visible by eye.

It is a sibling of SNPHeritabilityByLDSCTask rather than a reader of its output: it takes that
task purely to reach a single, consistent definition of the analysis (which summary statistics,
which LD reference, the preprocessing pipe, the sample size), and re-runs LD-score regression on
that definition to draw the fitted line.

The fit is re-run with gwaslab -- the same estimator the heritability task uses -- so the line and
the annotated intercept/heritability match what that task reports, rather than an independent
solver that differs in its treatment of extreme chi-square (gwaslab keeps all variants and uses a
two-step intercept estimator with a cutoff at chi-square 30; the repo's GenomicSEM port cuts
variants above max(0.001*N, 80) and fits a single-step intercept). The re-run passes NaN
prevalences so gwaslab reports observed-scale heritability, sidestepping the liability-scale
rescaling it applies to binary traits, which would otherwise make the fitted slope hard to
reconstruct. Consistently, the binned points keep every merged variant -- no high-chi-square cut --
to match the variant set gwaslab's slope regression uses.

The left y-axis is chi-square; the right y-axis is the same points rescaled to chi^2 * M / N, in
which the fitted slope equals the observed-scale heritability and the axis is comparable across
traits of different sample size. The two axes are one affine relabeling of a single set of marks,
not two independent measures.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Protocol

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from attrs import field, frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.build_system.meta.read_spec.read_sumstats import read_sumstats
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.consolidate_ld_scores_task import (
    LD_SCORE_LD_SCORE_COL,
    LD_SCORE_RSID_COL,
    read_ld_scores,
    total_m_5_50,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    PhenotypeInfo,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_snp_heritability_by_ldsc_task import (
    SNPHeritabilityByLDSCTask,
)
from mecfs_bio.build_system.task.gwaslab.ldsc_degenerate_z import (
    drop_variants_with_degenerate_z,
)
from mecfs_bio.build_system.task.gwaslab.ldsc_diagnostic import (
    LdscDiagnosticBins,
    bin_by_ld_score,
    chi2_to_heritability_units,
    fit_line_chi2,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.ldsc_constants import LDSC_Z_COL
from mecfs_bio.util.plotting.save_fig import write_plots_to_dir
from mecfs_bio.util.type_related.unwrap import unwrap

# Column names in gwaslab's parsed univariate LDSC summary (parse_ldsc_summary): the regression
# intercept and the observed-scale heritability. Their values are strings.
_GWASLAB_LDSC_INTERCEPT_COL = "Intercept"
_GWASLAB_LDSC_H2_OBS_COL = "h2_obs"


class SumstatsData(Protocol):
    """The slice of a gwaslab Sumstats the task reads directly: its mutable variant table. Kept
    narrow so the production reader can hand back a real gwaslab Sumstats while tests inject a
    lightweight stand-in."""

    data: pd.DataFrame


class LdscEstimableSumstats(Protocol):
    """The slice gwaslab_observed_fit drives: enough to run observed-scale LD-score regression and
    read back its summary. gwaslab's Sumstats satisfies this at runtime; a test can supply a fake
    that records the call and serves a canned summary."""

    meta: dict
    ldsc_h2: pd.DataFrame | None

    def infer_build(self) -> None: ...

    def estimate_h2_by_ldsc(
        self,
        *,
        ref_ld_chr: str,
        w_ld_chr: str,
        samp_prev: float | None,
        pop_prev: float | None,
    ) -> None: ...


# One data series (the binned points) plus two reference lines, so the three are separated by
# mark type -- markers, a solid line, a dashed line -- as well as by these colorblind-safe hues.
_POINT_COLOR = "#1f77b4"
_FIT_COLOR = "#d62728"
_REFERENCE_COLOR = "#7f7f7f"


@frozen
class LdscDiagnosticPlotConfig:
    """How the diagnostic is binned and framed.

    x_range and chi2_m_over_n_range are left unset for per-trait autoscaling; pin them to a shared
    window to make a set of plots comparable across traits (equal LD-score and heritability-unit
    axes). chi2_m_over_n_range is given in the right-axis units (chi^2 * M / N)."""

    n_bins: int = 20
    show_error_bars: bool = False
    x_range: tuple[float, float] | None = None
    chi2_m_over_n_range: tuple[float, float] | None = None
    title: str | None = None


@frozen
class LdscFit:
    """The observed-scale LD-score-regression fit needed to draw the line and annotate it."""

    intercept: float
    h2_obs: float


@frozen
class MergedDiagnosticData:
    """Per-variant chi-square and LD score for the variants shared by the summary statistics and
    the LD reference. Both arrays are 1-D of the same length, aligned row-for-row."""

    chi2: np.ndarray
    ld: np.ndarray

    def __attrs_post_init__(self) -> None:
        assert self.chi2.ndim == 1, f"chi2 must be 1-D, got shape {self.chi2.shape}"
        assert self.chi2.shape == self.ld.shape, (
            f"chi2 and ld must have the same shape, got {self.chi2.shape} and {self.ld.shape}"
        )


def compute_chi2(sumstats: pd.DataFrame) -> np.ndarray:
    """The per-variant chi-square. Prefer an existing Z column -- the same precedence
    drop_variants_with_degenerate_z uses -- and otherwise form Z from BETA / SE."""
    if LDSC_Z_COL in sumstats.columns:
        z = sumstats[LDSC_Z_COL].to_numpy(dtype=float)
    elif GWASLAB_BETA_COL in sumstats.columns and GWASLAB_SE_COL in sumstats.columns:
        z = sumstats[GWASLAB_BETA_COL].to_numpy(dtype=float) / sumstats[
            GWASLAB_SE_COL
        ].to_numpy(dtype=float)
    else:
        raise AssertionError(
            f"cannot compute chi-square: need a {LDSC_Z_COL} column, or both "
            f"{GWASLAB_BETA_COL} and {GWASLAB_SE_COL}"
        )
    return z * z


def merge_chi2_with_ld_scores(
    sumstats: pd.DataFrame, ld_scores: pd.DataFrame
) -> MergedDiagnosticData:
    """Inner-join the summary statistics to the LD reference on rsID, pairing each variant's
    chi-square with its LD score. Variants missing from either side drop out, which restricts the
    plot to the LD reference's variant set (HapMap3, in the standard reference)."""
    joined = pd.DataFrame(
        {
            GWASLAB_RSID_COL: sumstats[GWASLAB_RSID_COL].to_numpy(),
            "chi2": compute_chi2(sumstats),
        }
    ).merge(
        ld_scores[[LD_SCORE_RSID_COL, LD_SCORE_LD_SCORE_COL]],
        left_on=GWASLAB_RSID_COL,
        right_on=LD_SCORE_RSID_COL,
        how="inner",
    )
    return MergedDiagnosticData(
        chi2=joined["chi2"].to_numpy(dtype=float),
        ld=joined[LD_SCORE_LD_SCORE_COL].to_numpy(dtype=float),
    )


def gwaslab_observed_fit(
    sumstats: LdscEstimableSumstats, ref_ld_chr: str, build: GenomeBuild
) -> LdscFit:
    """Run gwaslab's LD-score regression on the observed scale and return its intercept and
    heritability -- the same estimator the heritability task uses, so the drawn line matches what
    that task reports. None (not NaN) prevalences keep gwaslab on the observed scale: its summary
    labels the result 'Liability' whenever both prevalences are not None, which renames the h2
    column, whereas None reports plain observed-scale h2 -- the scale in which the chi-square-units
    slope is reconstructable.

    sumstats is a gwaslab Sumstats already carrying the analysis's preprocessing and sample size."""
    sumstats.infer_build()
    assert sumstats.meta["gwaslab"]["genome_build"] == build, (
        f"sumstats build {sumstats.meta['gwaslab']['genome_build']} does not match {build}"
    )
    sumstats.estimate_h2_by_ldsc(
        ref_ld_chr=ref_ld_chr,
        w_ld_chr=ref_ld_chr,
        samp_prev=None,
        pop_prev=None,
    )
    ldsc_h2 = unwrap(sumstats.ldsc_h2)
    assert isinstance(ldsc_h2, pd.DataFrame)
    return LdscFit(
        intercept=float(ldsc_h2[_GWASLAB_LDSC_INTERCEPT_COL].iloc[0]),
        h2_obs=float(ldsc_h2[_GWASLAB_LDSC_H2_OBS_COL].iloc[0]),
    )


def resolve_sample_size(
    sumstats: pd.DataFrame, set_n: int | None, phenotype_info: PhenotypeInfo
) -> float:
    """The scalar N the fit and the right-axis rescaling use, mirroring how the heritability task
    sources it: an explicit override, then a per-variant N column (averaged), then the phenotype's
    total sample size."""
    if set_n is not None:
        return float(set_n)
    if GWASLAB_SAMPLE_SIZE_COLUMN in sumstats.columns:
        return float(sumstats[GWASLAB_SAMPLE_SIZE_COLUMN].mean())
    assert phenotype_info.total_sample_size is not None, (
        "no sample size available: set a sample size on the heritability task, provide an "
        f"{GWASLAB_SAMPLE_SIZE_COLUMN} column, or a phenotype total_sample_size"
    )
    return float(phenotype_info.total_sample_size)


def build_diagnostic_figure(
    bins: LdscDiagnosticBins,
    fit: LdscFit,
    n: float,
    m: float,
    config: LdscDiagnosticPlotConfig,
) -> go.Figure:
    """Assemble the interactive figure: binned mean chi-square, the fitted line, a chi^2 = 1
    reference, and a right-hand axis in chi^2 * M / N units whose slope reads as heritability."""
    x_min = config.x_range[0] if config.x_range else float(bins.mean_ld.min())
    x_max = config.x_range[1] if config.x_range else float(bins.mean_ld.max())
    ld_line = np.array([x_min, x_max])
    y_fit = fit_line_chi2(fit.intercept, fit.h2_obs, n, m, ld_line)

    if config.chi2_m_over_n_range is not None:
        secondary_range = list(config.chi2_m_over_n_range)
        primary_range = [v * n / m for v in secondary_range]
    else:
        candidates = np.concatenate([bins.mean_chi2, y_fit, [1.0]])
        low, high = float(candidates.min()), float(candidates.max())
        pad = 0.05 * (high - low) if high > low else 0.05 * high
        primary_range = [low - pad, high + pad]
        secondary_range = [v * m / n for v in primary_range]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=bins.mean_ld,
            y=bins.mean_chi2,
            mode="markers",
            name="Binned mean χ²",
            marker=dict(size=9, color=_POINT_COLOR),
            error_y=dict(
                type="data", array=bins.se_chi2, visible=config.show_error_bars
            ),
            customdata=np.stack(
                [chi2_to_heritability_units(bins.mean_chi2, n, m), bins.count], axis=-1
            ),
            hovertemplate=(
                "LD score = %{x:.2f}<br>"
                "χ² = %{y:.4f}<br>"
                "χ²·M/N = %{customdata[0]:.3e}<br>"
                "variants in bin = %{customdata[1]}<extra></extra>"
            ),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=ld_line,
            y=y_fit,
            mode="lines",
            name="LDSC fit",
            line=dict(color=_FIT_COLOR, width=2),
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[x_min, x_max],
            y=[1.0, 1.0],
            mode="lines",
            name="χ² = 1 (no confounding)",
            line=dict(color=_REFERENCE_COLOR, width=2, dash="dash"),
            hoverinfo="skip",
        )
    )
    fig.update_layout(
        template="plotly_white",
        title=config.title,
        xaxis=dict(title="LD score", range=[x_min, x_max]),
        yaxis=dict(title="mean χ²", range=primary_range),
        yaxis2=dict(
            title="χ² · M / N  (slope = h²)",
            overlaying="y",
            side="right",
            range=secondary_range,
            showgrid=False,
        ),
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        annotations=[
            dict(
                xref="paper",
                yref="paper",
                x=0.99,
                y=0.02,
                xanchor="right",
                yanchor="bottom",
                showarrow=False,
                align="right",
                bgcolor="rgba(255,255,255,0.8)",
                text=(
                    f"h² (observed) = {fit.h2_obs:.4f}<br>"
                    f"intercept = {fit.intercept:.4f}"
                ),
            )
        ],
    )
    return fig


def _ensure_sample_size_column(
    data: pd.DataFrame, set_n: int | None, phenotype_info: PhenotypeInfo
) -> None:
    """Set the N column on the sumstats the way the heritability task does before its gwaslab run,
    so gwaslab sees the same per-variant sample size: an explicit override wins, otherwise fall
    back to the phenotype total only when no N column is present."""
    if set_n is not None:
        data[GWASLAB_SAMPLE_SIZE_COLUMN] = set_n
    elif GWASLAB_SAMPLE_SIZE_COLUMN not in data.columns:
        assert phenotype_info.total_sample_size is not None, (
            "no sample size available: set a sample size on the heritability task, provide an "
            f"{GWASLAB_SAMPLE_SIZE_COLUMN} column, or a phenotype total_sample_size"
        )
        data[GWASLAB_SAMPLE_SIZE_COLUMN] = phenotype_info.total_sample_size


def _read_sumstats(asset: Asset) -> SumstatsData:
    return read_sumstats(asset)


@frozen
class LdscDiagnosticPlotTask(Task):
    """Draw the LD-score-regression diagnostic plot for a SNPHeritabilityByLDSCTask analysis."""

    meta: Meta
    ldsc_task: SNPHeritabilityByLDSCTask
    config: LdscDiagnosticPlotConfig
    sumstats_reader: Callable[[Asset], SumstatsData] = field(default=_read_sumstats)
    fit_estimator: Callable[..., LdscFit] = field(default=gwaslab_observed_fit)

    @property
    def _source_sumstats_id(self) -> AssetId:
        return self.ldsc_task.source_sumstats_task.asset_id

    @property
    def _ld_ref_id(self) -> AssetId:
        return self.ldsc_task.ld_ref_task.asset_id

    @property
    def deps(self) -> list[Task]:
        return [self.ldsc_task.source_sumstats_task, self.ldsc_task.ld_ref_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        sumstats = self.sumstats_reader(fetch(self._source_sumstats_id))
        sumstats.data = self.ldsc_task.pipe.process_pandas(sumstats.data)
        sumstats.data = drop_variants_with_degenerate_z(sumstats.data)
        _ensure_sample_size_column(
            sumstats.data, self.ldsc_task.set_N, self.ldsc_task.phenotype_info
        )
        n = resolve_sample_size(
            sumstats.data, self.ldsc_task.set_N, self.ldsc_task.phenotype_info
        )

        ld_asset = fetch(self._ld_ref_id)
        assert isinstance(ld_asset, DirectoryAsset)
        ld_frame = read_ld_scores(ld_asset.path).collect()
        total_m = total_m_5_50(ld_frame)

        # Bin from the processed sumstats before the fit: gwaslab's estimate mutates the Sumstats,
        # and the bins must reflect the variants as they enter the regression.
        merged = merge_chi2_with_ld_scores(sumstats.data, ld_frame.to_pandas())
        bins = bin_by_ld_score(merged.chi2, merged.ld, self.config.n_bins)

        ref_ld_chr = str(ld_asset.path) + self.ldsc_task.ld_file_filename_pattern
        fit = self.fit_estimator(sumstats, ref_ld_chr, self.ldsc_task.build)

        fig = build_diagnostic_figure(
            bins=bins, fit=fit, n=n, m=total_m, config=self.config
        )
        out_dir = scratch_dir / "ldsc_diagnostic_plot"
        write_plots_to_dir(out_dir, {"ldsc_diagnostic": fig})
        return DirectoryAsset(out_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        ldsc_task: SNPHeritabilityByLDSCTask,
        config: LdscDiagnosticPlotConfig = LdscDiagnosticPlotConfig(),
        sumstats_reader: Callable[[Asset], SumstatsData] = _read_sumstats,
        fit_estimator: Callable[..., LdscFit] = gwaslab_observed_fit,
    ) -> LdscDiagnosticPlotTask:
        source_meta = ldsc_task.meta
        assert isinstance(source_meta, ResultTableMeta)
        meta = GWASPlotDirectoryMeta(
            trait=source_meta.trait,
            project=source_meta.project,
            id=AssetId(asset_id),
        )
        return cls(
            meta=meta,
            ldsc_task=ldsc_task,
            config=config,
            sumstats_reader=sumstats_reader,
            fit_estimator=fit_estimator,
        )
