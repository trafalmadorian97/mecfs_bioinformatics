"""Fit a ridge surrogate of the polyfun prior on the baseline-LF annotations.

Regresses snpvar_bin (the polyfun binned per-SNP heritability prior actually used
in fine mapping, not the original S-LDSC tau_c) on the 187 annotations, genome
wide. The fit is done from per-chromosome cross-product sufficient statistics, so
the full design matrix is never held in memory; alpha is chosen by
leave-one-chromosome-out. Outputs raw-scale coefficients gamma_raw (used by the
explainability contrast) and standardized coefficients gamma_standardized (for
global importance ranking), plus each annotation's family.
"""

import json
from pathlib import Path

import numpy as np
import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.annotation_weights.build_baseline_lf_annotation_parquet_task import (
    ANNOT_KEY_COLUMNS,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.polyfun_annotation_families import family_for_annotation

logger = structlog.get_logger()

WEIGHTS_PARQUET_FILENAME = "weights.parquet"
DIAGNOSTICS_JSON_FILENAME = "diagnostics.json"
ANNOTATION_COL = "annotation"
GAMMA_RAW_COL = "gamma_raw"
GAMMA_STANDARDIZED_COL = "gamma_standardized"
FAMILY_COL = "family"
SNPVAR_COL = "snpvar_bin"
_CHR_COL = "CHR"
_SNP_COL = "SNP"

_DEFAULT_ALPHAS: tuple[float, ...] = (0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0)


@frozen
class _ChromStats:
    """Cross-product sufficient statistics for one chromosome (standardized-free,
    raw annotation scale). p = number of annotations."""

    n: int
    sx: np.ndarray  # (p,) sum of annotations
    sxx: np.ndarray  # (p, p) sum of a a^T
    sxy: np.ndarray  # (p,) sum of a * y
    sy: float
    syy: float


@frozen
class RidgeAnnotationWeightsTask(Task):
    meta: Meta
    annotation_parquet_task: Task
    snpvar_meta_task: Task
    alphas: tuple[float, ...] = _DEFAULT_ALPHAS

    @property
    def deps(self) -> list["Task"]:
        return [self.annotation_parquet_task, self.snpvar_meta_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        annot_asset = fetch(self.annotation_parquet_task.asset_id)
        assert isinstance(annot_asset, FileAsset)
        meta_asset = fetch(self.snpvar_meta_task.asset_id)
        meta = (
            scan_dataframe_asset(meta_asset, self.snpvar_meta_task.meta)
            .select(_SNP_COL, SNPVAR_COL)
            .collect()
            .to_polars()
            .unique(subset=_SNP_COL)
        )

        annot_columns = _annotation_columns(annot_asset.path)
        per_chrom = _accumulate_per_chromosome(annot_asset.path, annot_columns, meta)
        alpha, mean_r2, r2_per_chrom = _select_alpha_loco(per_chrom, self.alphas)
        gamma_raw, gamma_std, intercept = _fit_all(per_chrom, alpha)

        weights = pl.DataFrame(
            {
                ANNOTATION_COL: annot_columns,
                GAMMA_RAW_COL: gamma_raw,
                GAMMA_STANDARDIZED_COL: gamma_std,
                FAMILY_COL: [family_for_annotation(c) for c in annot_columns],
            }
        )
        weights.write_parquet(scratch_dir / WEIGHTS_PARQUET_FILENAME)
        diagnostics = {
            "alpha": alpha,
            "intercept": intercept,
            "heldout_r2_per_chrom": r2_per_chrom,
            "mean_heldout_r2": mean_r2,
            "n_variants": int(sum(s.n for s in per_chrom.values())),
        }
        (scratch_dir / DIAGNOSTICS_JSON_FILENAME).write_text(
            json.dumps(diagnostics, indent=2, sort_keys=True)
        )
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        annotation_parquet_task: Task,
        snpvar_meta_task: Task,
        alphas: tuple[float, ...] = _DEFAULT_ALPHAS,
    ) -> "RidgeAnnotationWeightsTask":
        # Derive the output directory meta from the primary dependency (the
        # annotation parquet), reusing group/sub_group/sub_folder - the
        # CompressedCSVToParquetTask.create pattern applied to a DirMeta.
        source_meta = annotation_parquet_task.meta
        if not isinstance(source_meta, ReferenceFileMeta):
            raise ValueError(f"Unknown meta for annotation parquet task: {source_meta}")
        meta = ReferenceDataDirectoryMeta(
            group=source_meta.group,
            sub_group=source_meta.sub_group,
            sub_folder=source_meta.sub_folder,
            id=AssetId(asset_id),
        )
        return cls(
            meta=meta,
            annotation_parquet_task=annotation_parquet_task,
            snpvar_meta_task=snpvar_meta_task,
            alphas=alphas,
        )


def _annotation_columns(annot_path: Path) -> list[str]:
    schema = pl.scan_parquet(annot_path).collect_schema()
    return [c for c in schema.names() if c not in ANNOT_KEY_COLUMNS]


def _accumulate_per_chromosome(
    annot_path: Path, annot_columns: list[str], meta: pl.DataFrame
) -> dict[int, _ChromStats]:
    chroms = (
        pl.scan_parquet(annot_path)
        .select(_CHR_COL)
        .unique()
        .collect()
        .to_series()
        .to_list()
    )
    per_chrom: dict[int, _ChromStats] = {}
    for chrom in sorted(chroms):
        frame = (
            pl.scan_parquet(annot_path)
            .filter(pl.col(_CHR_COL) == chrom)
            .collect()
            .join(meta, on=_SNP_COL, how="inner")
        )
        x = frame.select(annot_columns).to_numpy().astype(np.float64)
        y = frame.select(SNPVAR_COL).to_numpy().ravel().astype(np.float64)
        per_chrom[chrom] = _ChromStats(
            n=x.shape[0],
            sx=x.sum(0),
            sxx=x.T @ x,
            sxy=x.T @ y,
            sy=float(y.sum()),
            syy=float(y @ y),
        )
    return per_chrom


def _combine(stats: list[_ChromStats]) -> _ChromStats:
    """Sum per-chromosome sufficient statistics into one aggregate."""
    return _ChromStats(
        n=sum(s.n for s in stats),
        sx=sum((s.sx for s in stats), start=np.zeros_like(stats[0].sx)),
        sxx=sum((s.sxx for s in stats), start=np.zeros_like(stats[0].sxx)),
        sxy=sum((s.sxy for s in stats), start=np.zeros_like(stats[0].sxy)),
        sy=sum(s.sy for s in stats),
        syy=sum(s.syy for s in stats),
    )


@frozen
class _StandardizedSystem:
    """The centered+standardized ridge system for a set of annotations."""

    g_std: np.ndarray  # (p, p) standardized Gram
    b_std: np.ndarray  # (p,) standardized cross-term
    mean: np.ndarray  # (p,) per-annotation mean
    sd: np.ndarray  # (p,) per-annotation std (zeros replaced by 1)
    mean_y: float


def _standardized_system(stats: _ChromStats) -> _StandardizedSystem:
    """Build the centered+standardized ridge system from raw cross-products."""
    n = stats.n
    mean = stats.sx / n
    var = np.diag(stats.sxx) / n - mean**2
    sd = np.sqrt(np.maximum(var, 0.0))
    sd[sd == 0] = 1.0
    centered_gram = stats.sxx - n * np.outer(mean, mean)
    g_std = centered_gram / np.outer(sd, sd)
    b_std = (stats.sxy - mean * stats.sy) / sd
    return _StandardizedSystem(
        g_std=g_std, b_std=b_std, mean=mean, sd=sd, mean_y=stats.sy / n
    )


def _solve(system: _StandardizedSystem, alpha: float) -> np.ndarray:
    """Solve (G_std + alpha I) gamma_std = b_std for the standardized coefficients."""
    p = system.g_std.shape[0]
    return np.linalg.solve(system.g_std + alpha * np.eye(p), system.b_std)


def _heldout_r2(
    held: _ChromStats,
    gamma_std: np.ndarray,
    train_mean: np.ndarray,
    train_sd: np.ndarray,
    train_mean_y: float,
) -> float:
    """R^2 of the standardized model on a held-out chromosome, standardizing the
    held-out annotations with the TRAIN mean/sd.

    Prediction for raw annotations x is pred = train_mean_y + gamma_std . z, where
    z = (x - train_mean) / train_sd. Everything is expanded from the held-out
    chromosome's raw cross-product sufficient statistics; no per-SNP data needed.
    """
    n = held.n
    c = train_mean_y
    tm = train_mean
    ts = train_sd
    # SS over (y - c): sum (y - c)^2
    ss_res_y = held.syy - 2.0 * c * held.sy + n * c * c
    # sum z_i (y_i - c)  and  sum z_i z_i^T  in terms of raw stats
    z_r = (held.sxy - c * held.sx - tm * held.sy + n * c * tm) / ts
    zz = (
        held.sxx - np.outer(tm, held.sx) - np.outer(held.sx, tm) + n * np.outer(tm, tm)
    ) / np.outer(ts, ts)
    ss_res = ss_res_y - 2.0 * float(gamma_std @ z_r) + float(gamma_std @ zz @ gamma_std)
    mean_y_held = held.sy / n
    ss_tot = held.syy - n * mean_y_held * mean_y_held
    return 1.0 - ss_res / ss_tot


def _select_alpha_loco(
    per_chrom: dict[int, _ChromStats], alphas: tuple[float, ...]
) -> tuple[float, float, dict[str, float]]:
    """Pick alpha by leave-one-chromosome-out; return (alpha, mean_r2, r2_per_chrom)."""
    chroms = sorted(per_chrom)
    best_alpha = alphas[0]
    best_mean = -np.inf
    best_r2_per_chrom: dict[str, float] = {}
    for alpha in alphas:
        r2s: dict[str, float] = {}
        for held in chroms:
            train = _combine([per_chrom[c] for c in chroms if c != held])
            system = _standardized_system(train)
            gamma_std = _solve(system, alpha)
            r2s[str(held)] = _heldout_r2(
                per_chrom[held],
                gamma_std=gamma_std,
                train_mean=system.mean,
                train_sd=system.sd,
                train_mean_y=system.mean_y,
            )
        mean_r2 = float(np.mean(list(r2s.values())))
        if mean_r2 > best_mean:
            best_mean, best_alpha, best_r2_per_chrom = mean_r2, alpha, r2s
    return best_alpha, best_mean, best_r2_per_chrom


def _fit_all(
    per_chrom: dict[int, _ChromStats], alpha: float
) -> tuple[list[float], list[float], float]:
    """Fit on all chromosomes; return (gamma_raw, gamma_standardized, intercept)."""
    system = _standardized_system(_combine(list(per_chrom.values())))
    gamma_std = _solve(system, alpha)
    gamma_raw = gamma_std / system.sd
    intercept = float(system.mean_y - gamma_raw @ system.mean)
    return gamma_raw.tolist(), gamma_std.tolist(), intercept
