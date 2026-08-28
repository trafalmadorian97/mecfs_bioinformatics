"""Stacked explainability figure for a polyfun-vs-uniform SUSIE result.

Panels, top to bottom, sharing the genomic-position x-axis:
  1. Manhattan (-log10 p), points colored by LD r^2 with the min-p lead
     variant, with local recombination rate (cM/Mb from the hg19 genetic map)
     on a secondary axis.
  2..(1+x). One panel per important family: sum_c gamma_raw_c * a_ic across
     the locus (raw scaled value; the contrast is the profile, not plotted
     here).
  (2+x). Prior fold m*pi_i (log scale).
  (3+x). PIP, uniform run.
  (4+x). PIP, polyfun run.
  (5+x). Genes.

Writes both explain_plot.png and explain_plot.svg. Inspired by
SusieStackPlotTask but independent of it.

The important families and the focal variant come from the contrast task's
selection.json, so the figure agrees with the contrast task's tables. The
per-variant family-scaled tracks are recomputed here (rather than read back
from the contrast task's output) directly from the ridge weights and the
annotation parquet, using the same allele-aware join and the same
_family_scaled helper the contrast task itself uses, so the values plotted
match the contrast task's family_scaled table exactly.
"""

import json
from pathlib import Path, PurePath

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import narwhals as nw  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
from attrs import frozen  # noqa: E402
from matplotlib import gridspec  # noqa: E402

from mecfs_bio.build_system.asset.base_asset import Asset  # noqa: E402
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset  # noqa: E402
from mecfs_bio.build_system.meta.asset_id import AssetId  # noqa: E402
from mecfs_bio.build_system.meta.meta import Meta  # noqa: E402
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (  # noqa: E402
    scan_dataframe_asset,
)
from mecfs_bio.build_system.meta.result_directory_meta import (  # noqa: E402
    ResultDirectoryMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch  # noqa: E402
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (  # noqa: E402
    ANNOTATION_COL,
    GAMMA_RAW_COL,
)
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (  # noqa: E402
    FAMILY_COL as WEIGHTS_FAMILY_COL,
)
from mecfs_bio.build_system.task.base_task import Task  # noqa: E402
from mecfs_bio.build_system.task.genetic_map.parse_genetic_map_task import (  # noqa: E402
    GMAP_POS_COL,
    GMAP_RATE_COL,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import (  # noqa: E402
    DataProcessingPipe,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe  # noqa: E402
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (  # noqa: E402
    _ANNOT_KEY,
    FAMILY_COL,
    FAMILY_SCALED_COL,
    SELECTION_JSON_FILENAME,
    _family_scaled,
    _load_annotations,
    _load_run_variants,
    _load_weights,
)
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (  # noqa: E402
    FILTERED_GWAS_FILENAME,
    FILTERED_LD_FILENAME,
    PIP_COLUMN,
    PRIOR_FILENAME,
    PRIOR_WEIGHT_COLUMN,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (  # noqa: E402
    GENE_INFO_END_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
)
from mecfs_bio.build_system.wf.base_wf import WF  # noqa: E402
from mecfs_bio.constants.gwaslab_constants import (  # noqa: E402
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)

PLOT_PNG_FILENAME = "explain_plot.png"
PLOT_SVG_FILENAME = "explain_plot.svg"
# _ANNOT_KEY (imported from the contrast task) is the allele-aware join key
# (CHR, POS, unordered-allele-key); _load_annotations and _load_run_variants both
# attach that key, so the family-scaled values plotted here match the tables.


@frozen
class PolyfunExplainPlotTask(Task):
    """Render the 8-panel polyfun-vs-uniform explainability figure.

    Depends directly on the two SUSIE run dirs, the contrast task (read only
    for selection.json, so the plotted focal variant and important families
    match the contrast task's tables), the ridge weights, the annotation
    parquet, and a gene-info source. The per-variant family-scaled tracks are
    recomputed from the ridge weights and annotation parquet rather than read
    back from the contrast task, using the same helpers and join key the
    contrast task uses, so the values agree.
    """

    meta: Meta
    susie_uniform_task: Task
    susie_polyfun_task: Task
    contrast_task: Task
    annotation_parquet_task: Task
    gene_info_task: Task
    ridge_weights_task: Task
    genetic_map_task: Task
    gene_info_pipe: DataProcessingPipe = IdentityPipe()
    n_family_panels: int = 3

    @property
    def deps(self) -> list["Task"]:
        return [
            self.susie_uniform_task,
            self.susie_polyfun_task,
            self.contrast_task,
            self.ridge_weights_task,
            self.annotation_parquet_task,
            self.gene_info_task,
            self.genetic_map_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        uni_dir = _dir(fetch, self.susie_uniform_task)
        pf_dir = _dir(fetch, self.susie_polyfun_task)
        contrast_dir = _dir(fetch, self.contrast_task)
        selection = json.loads((contrast_dir / SELECTION_JSON_FILENAME).read_text())
        families = selection["important_families"][: self.n_family_panels]

        pf_variants = _load_run_variants(pf_dir).sort(GWASLAB_POS_COL)
        uni_variants = _load_run_variants(uni_dir)
        pf_full = pl.read_parquet(pf_dir / FILTERED_GWAS_FILENAME)
        ld = np.load(pf_dir / FILTERED_LD_FILENAME)
        prior_w = pl.read_parquet(pf_dir / PRIOR_FILENAME)[
            PRIOR_WEIGHT_COLUMN
        ].to_numpy()

        chrom = int(pf_variants[GWASLAB_CHROM_COL][0])
        bp_min = int(pf_variants[GWASLAB_POS_COL].to_numpy().min())
        bp_max = int(pf_variants[GWASLAB_POS_COL].to_numpy().max())

        weights = _load_weights(fetch, self.ridge_weights_task)
        annot_cols = weights[ANNOTATION_COL].to_list()
        annot = _load_annotations(
            fetch, self.annotation_parquet_task, chrom, bp_min, bp_max, annot_cols
        )
        # Allele-aware join on (CHR, POS, allele-key) -- mirrors the contrast
        # task's own join exactly, so the family-scaled values plotted here match
        # its per_family_contrast/family_scaled tables.
        pf_annot = pf_variants.join(annot, on=_ANNOT_KEY, how="inner")
        gamma = dict(zip(weights[ANNOTATION_COL], weights[GAMMA_RAW_COL]))
        family = dict(zip(weights[ANNOTATION_COL], weights[WEIGHTS_FAMILY_COL]))
        scaled = _family_scaled(pf_annot, annot_cols, gamma, family)

        genes = (
            self.gene_info_pipe.process(
                scan_dataframe_asset(
                    fetch(self.gene_info_task.asset_id), self.gene_info_task.meta
                )
            )
            .collect()
            .to_polars()
        )

        _render(
            scratch_dir=scratch_dir,
            pf_variants=pf_variants,
            uni_variants=uni_variants,
            pf_full=pf_full,
            ld=ld,
            prior_w=prior_w,
            families=families,
            scaled=scaled,
            genes=genes,
            fetch=fetch,
            genetic_map_task=self.genetic_map_task,
            chrom=chrom,
            bp_min=bp_min,
            bp_max=bp_max,
        )
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        susie_uniform_task: Task,
        susie_polyfun_task: Task,
        contrast_task: Task,
        annotation_parquet_task: Task,
        gene_info_task: Task,
        ridge_weights_task: Task,
        genetic_map_task: Task,
        gene_info_pipe: DataProcessingPipe = IdentityPipe(),
        n_family_panels: int = 3,
    ) -> "PolyfunExplainPlotTask":
        source_meta = susie_polyfun_task.meta
        if not isinstance(source_meta, ResultDirectoryMeta):
            raise ValueError(f"Unknown meta for polyfun susie task: {source_meta}")
        meta = ResultDirectoryMeta(
            id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=PurePath("analysis"),
        )
        return cls(
            meta=meta,
            susie_uniform_task=susie_uniform_task,
            susie_polyfun_task=susie_polyfun_task,
            contrast_task=contrast_task,
            annotation_parquet_task=annotation_parquet_task,
            gene_info_task=gene_info_task,
            ridge_weights_task=ridge_weights_task,
            genetic_map_task=genetic_map_task,
            gene_info_pipe=gene_info_pipe,
            n_family_panels=n_family_panels,
        )


def _dir(fetch: Fetch, task: Task) -> Path:
    asset = fetch(task.asset_id)
    assert isinstance(asset, DirectoryAsset)
    return asset.path


def _norm_sf(z: np.ndarray) -> np.ndarray:
    from scipy.stats import norm

    return norm.sf(z)


def _load_recomb(
    fetch: Fetch, task: Task, chrom: int, bp_min: int, bp_max: int
) -> pl.DataFrame:
    """Locus-windowed recombination rate (cM/Mb) from the hg19 genetic map."""
    return (
        scan_dataframe_asset(fetch(task.asset_id), task.meta)
        .filter(
            (nw.col(GWASLAB_CHROM_COL) == chrom)
            & (nw.col(GMAP_POS_COL) >= bp_min)
            & (nw.col(GMAP_POS_COL) <= bp_max)
        )
        .select(GWASLAB_CHROM_COL, GMAP_POS_COL, GMAP_RATE_COL)
        .collect()
        .to_polars()
        .sort(GMAP_POS_COL)
    )


def _plot_recomb(ax, recomb: pl.DataFrame) -> None:
    bp = recomb[GMAP_POS_COL].to_numpy().astype(float)
    rate = recomb[GMAP_RATE_COL].to_numpy().astype(float)
    if len(bp) >= 1:
        ax.plot(bp, rate, color="tab:red", alpha=0.4, linewidth=1)
    ax.set_ylabel("cM/Mb", color="tab:red")


def _plot_genes(ax, genes: pl.DataFrame, bp_min: int, bp_max: int) -> None:
    row = 0
    for g in genes.iter_rows(named=True):
        start = max(int(g[GENE_INFO_START_COL]), bp_min)
        end = min(int(g[GENE_INFO_END_COL]), bp_max)
        ax.plot([start, end], [row, row], linewidth=4)
        ax.text((start + end) / 2, row + 0.1, g[GENE_INFO_NAME_COL], ha="center")
        row += 1
    ax.set_yticks([])
    ax.set_ylabel("genes")


def _render(
    scratch_dir: Path,
    pf_variants: pl.DataFrame,
    uni_variants: pl.DataFrame,
    pf_full: pl.DataFrame,
    ld: np.ndarray,
    prior_w: np.ndarray,
    families: list[str],
    scaled: pl.DataFrame,
    genes: pl.DataFrame,
    fetch: Fetch,
    genetic_map_task: Task,
    chrom: int,
    bp_min: int,
    bp_max: int,
) -> None:
    n_panels = 1 + len(families) + 3 + 1  # manhattan + families + lift + 2 pip + genes
    fig = plt.figure(figsize=(10, 2.0 * n_panels))
    gs = gridspec.GridSpec(nrows=n_panels, ncols=1, hspace=0.4)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(n_panels)]
    x = pf_full[GWASLAB_POS_COL].to_numpy()

    # Panel 1: Manhattan colored by LD with lead (min-p) variant + recomb rate.
    z = (pf_full[GWASLAB_BETA_COL] / pf_full[GWASLAB_SE_COL]).to_numpy()
    neglogp = -np.log10(2.0 * _norm_sf(np.abs(z)))
    lead = int(np.argmax(np.abs(z)))
    r2 = ld[lead, :] ** 2
    ax0 = axes[0]
    sc = ax0.scatter(x, neglogp, c=r2, cmap="viridis", vmin=0, vmax=1, s=12)
    ax0.set_ylabel("-log10 p")
    fig.colorbar(sc, ax=ax0, label="r2 w/ lead")

    recomb = _load_recomb(fetch, genetic_map_task, chrom, bp_min, bp_max)
    ax0b = ax0.twinx()
    _plot_recomb(ax0b, recomb)

    # Family panels: sum_c gamma_raw_c * a_ic across the locus.
    for k, fam in enumerate(families):
        axk = axes[1 + k]
        fam_series = scaled.filter(pl.col(FAMILY_COL) == fam).sort(GWASLAB_POS_COL)
        axk.plot(
            fam_series[GWASLAB_POS_COL].to_numpy(),
            fam_series[FAMILY_SCALED_COL].to_numpy(),
        )
        axk.set_ylabel(fam)

    # Prior fold (log).
    lift = pf_variants.height * prior_w / prior_w.sum()
    ax_lift = axes[1 + len(families)]
    ax_lift.scatter(pf_full[GWASLAB_POS_COL].to_numpy(), lift, s=12)
    ax_lift.set_yscale("log")
    ax_lift.set_ylabel("prior fold")

    # PIP uniform / polyfun.
    ax_u = axes[2 + len(families)]
    ax_u.scatter(
        uni_variants[GWASLAB_POS_COL].to_numpy(),
        uni_variants[PIP_COLUMN].to_numpy(),
        s=12,
    )
    ax_u.set_ylabel("PIP uniform")
    ax_pf = axes[3 + len(families)]
    ax_pf.scatter(
        pf_variants[GWASLAB_POS_COL].to_numpy(),
        pf_variants[PIP_COLUMN].to_numpy(),
        s=12,
    )
    ax_pf.set_ylabel("PIP polyfun")

    # Genes.
    _plot_genes(axes[-1], genes, bp_min, bp_max)
    axes[-1].set_xlabel(f"chr{chrom} position (bp)")

    fig.savefig(scratch_dir / PLOT_PNG_FILENAME, dpi=150, bbox_inches="tight")
    fig.savefig(scratch_dir / PLOT_SVG_FILENAME, bbox_inches="tight")
    plt.close(fig)
