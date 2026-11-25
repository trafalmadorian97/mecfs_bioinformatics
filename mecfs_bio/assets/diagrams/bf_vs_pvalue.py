"""
Diagram comparing the coloc approximate Bayes factor with the p value of the
same variant.

Plots log10 BF(i), the Wakefield approximate Bayes factor that coloc.abf
computes for one variant, against the two-sided -log10 p value, for a few fixed estimates beta_hat. Along
each curve the sample size N varies, so the sampling variance V runs from large
to small and both z and the Bayes factor grow.

Units: sd(y) = 1, so beta_hat and omega are in phenotypic standard deviations.

Glossary:
    beta_hat: estimated per-allele effect of the variant.
    omega: prior standard deviation of the true effect under the alternative.
    f: allele frequency.
    N: sample size.
    V = 1 / (2 N f (1 - f)): sampling variance of beta_hat under Hardy-Weinberg
        equilibrium.
    z = beta_hat / sqrt(V).
    r = omega^2 / (omega^2 + V): shrinkage factor.
    log BF = (log(1 - r) + r z^2) / 2, natural log.
    p = 2 (1 - Phi(|z|)), computed through the log survival function so that it
        does not underflow for large |z|.
"""

import attrs
import matplotlib
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from scipy.stats import norm

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.diagram_file_meta import DiagramFileMeta
from mecfs_bio.build_system.task.matplotlib_diagram_task import (
    MatplotlibDiagramTask,
)

PRIOR_SD_OMEGA = 0.15
ALLELE_FREQUENCY = 0.3
BETA_HATS = (0.02, 0.05, 0.1, 0.3)
MARKED_SAMPLE_SIZES = (1e3, 1e4, 1e5)
MARKER_SHAPES = ("o", "s", "^")
N_CURVE_POINTS = 4000
LOG10_SAMPLE_SIZE_RANGE = (1.0, 7.0)
GENOME_WIDE_NEG_LOG10_P = -np.log10(5e-8)
X_MAX = 12.0
Y_MIN = -1.0
Y_MAX = 10.0

# Reference categorical palette, light mode, slots 1-4 (validated for adjacent pairs).
SERIES_COLORS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID_COLOR = "#e4e3df"
REFERENCE_COLOR = "#8a8984"


@attrs.frozen(slots=True)
class CurvePoints:
    """
    Points (x, y) = (-log10 p, log10 BF) for one beta_hat, one entry per sample
    size. Both arrays have shape (n_sample_sizes,).
    """

    neg_log10_p: np.ndarray
    log10_bf: np.ndarray

    def __attrs_post_init__(self) -> None:
        assert self.neg_log10_p.ndim == 1, "neg_log10_p must be 1-dimensional"
        assert self.neg_log10_p.shape == self.log10_bf.shape, (
            "neg_log10_p and log10_bf must have the same shape"
        )


def sampling_variance(sample_size: np.ndarray, allele_frequency: float) -> np.ndarray:
    """
    V = 1 / (2 N f (1 - f)), the variance of beta_hat for sd(y) = 1 under
    Hardy-Weinberg equilibrium.
    """
    return 1.0 / (2.0 * sample_size * allele_frequency * (1.0 - allele_frequency))


def curve_points(
    beta_hat: float,
    sample_size: np.ndarray,
    prior_sd_omega: float,
    allele_frequency: float,
) -> CurvePoints:
    """
    Evaluate the approximate Bayes factor and the two-sided p value at each
    sample size for a fixed beta_hat.
    """
    variance = sampling_variance(sample_size, allele_frequency=allele_frequency)
    shrinkage_r = prior_sd_omega**2 / (prior_sd_omega**2 + variance)
    z = beta_hat / np.sqrt(variance)
    log_bf = 0.5 * (np.log1p(-shrinkage_r) + shrinkage_r * z**2)
    log_p = np.log(2.0) + norm.logsf(np.abs(z))
    return CurvePoints(
        neg_log10_p=-log_p / np.log(10.0), log10_bf=log_bf / np.log(10.0)
    )


def _draw_reference_lines(ax: Axes) -> None:
    """
    BF = 1, BF = 1 / p, and genome-wide significance, all recessive.
    """
    ax.axhline(0.0, color=REFERENCE_COLOR, linewidth=0.8)
    ax.plot(
        [0.0, X_MAX],
        [0.0, X_MAX],
        color=REFERENCE_COLOR,
        linewidth=0.8,
        linestyle=(0, (4, 3)),
    )
    ax.text(
        8.75,
        9.05,
        r"$\mathrm{BF} = 1/p$",
        color=TEXT_SECONDARY,
        fontsize=8,
        ha="right",
    )
    ax.axvline(
        GENOME_WIDE_NEG_LOG10_P,
        color=REFERENCE_COLOR,
        linewidth=0.8,
        linestyle=(0, (1, 2)),
    )
    ax.text(
        GENOME_WIDE_NEG_LOG10_P + 0.12,
        Y_MIN + 0.25,
        r"$p = 5\times10^{-8}$",
        color=TEXT_SECONDARY,
        fontsize=8,
    )


def _draw_curve(
    ax: Axes,
    beta_hat: float,
    color: str,
    prior_sd_omega: float,
    allele_frequency: float,
) -> None:
    """
    One beta_hat: the curve over sample size, plus markers at the marked sample
    sizes that fall inside the axes.
    """
    curve = curve_points(
        beta_hat,
        np.logspace(*LOG10_SAMPLE_SIZE_RANGE, N_CURVE_POINTS),
        prior_sd_omega=prior_sd_omega,
        allele_frequency=allele_frequency,
    )
    ax.plot(
        curve.neg_log10_p,
        curve.log10_bf,
        color=color,
        linewidth=1.6,
        label=rf"$\hat\beta = {beta_hat}$",
    )
    marked = curve_points(
        beta_hat,
        np.array(MARKED_SAMPLE_SIZES),
        prior_sd_omega=prior_sd_omega,
        allele_frequency=allele_frequency,
    )
    for x, y, shape in zip(
        marked.neg_log10_p, marked.log10_bf, MARKER_SHAPES, strict=True
    ):
        if x <= X_MAX and Y_MIN <= y <= Y_MAX:
            ax.plot(
                x,
                y,
                marker=shape,
                markersize=5.5,
                color=color,
                markeredgecolor="white",
                markeredgewidth=0.8,
            )


def _draw_marker_legend(ax: Axes) -> None:
    """
    Legend entries for the sample-size markers, in neutral ink so they do not
    read as a series.
    """
    handles = [
        Line2D(
            [], [], linestyle="none", marker=shape, markersize=5.5, color=TEXT_SECONDARY
        )
        for shape in MARKER_SHAPES
    ]
    labels = [rf"$N = 10^{{{int(np.log10(n))}}}$" for n in MARKED_SAMPLE_SIZES]
    legend = ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.0, 0.72),
        frameon=False,
        fontsize=8,
    )
    ax.add_artist(legend)


def _style_axes(ax: Axes) -> None:
    ax.set_xlim(0.0, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_xlabel(r"$-\log_{10} p$ (two-sided)", color=TEXT_PRIMARY)
    ax.set_ylabel(r"$\log_{10}\,\mathrm{BF}(i)$", color=TEXT_PRIMARY)
    ax.grid(color=GRID_COLOR, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(REFERENCE_COLOR)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)


@attrs.frozen(slots=True)
class BayesFactorVsPValueFigure:
    """
    Builds the log10 BF versus -log10 p figure for fixed prior_sd_omega and
    allele_frequency.
    """

    prior_sd_omega: float = PRIOR_SD_OMEGA
    allele_frequency: float = ALLELE_FREQUENCY

    def __attrs_post_init__(self) -> None:
        assert self.prior_sd_omega > 0, "prior_sd_omega must be positive"
        assert 0 < self.allele_frequency < 1, "allele_frequency must lie in (0, 1)"

    def build(self) -> Figure:
        with matplotlib.rc_context(
            {"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9}
        ):
            fig = Figure(figsize=(6.3, 4.0))
            ax = fig.add_subplot()
            _draw_reference_lines(ax)
            for beta_hat, color in zip(BETA_HATS, SERIES_COLORS, strict=True):
                _draw_curve(
                    ax,
                    beta_hat,
                    color,
                    prior_sd_omega=self.prior_sd_omega,
                    allele_frequency=self.allele_frequency,
                )
            _draw_marker_legend(ax)
            ax.legend(loc="upper left", frameon=False, fontsize=8)
            _style_axes(ax)
            fig.tight_layout()
        return fig


BF_VS_PVALUE_DIAGRAM = MatplotlibDiagramTask(
    meta=DiagramFileMeta(AssetId("coloc_bf_vs_pvalue_diagram")),
    figure=BayesFactorVsPValueFigure(),
)
