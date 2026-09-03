"""Fit a ridge surrogate of the polyfun prior on the baseline-LF annotations.

Regresses snpvar_bin (the polyfun binned per-SNP heritability prior actually used
in fine mapping, not the original S-LDSC tau_c) on the 187 annotations, genome
wide. The fit is done from per-chromosome cross-product sufficient statistics, so
the full design matrix is never held in memory; alpha is chosen by
leave-one-chromosome-out. Outputs raw-scale coefficients gamma_raw (used by the
explainability contrast) and standardized coefficients gamma_standardized (for
global importance ranking), plus each annotation's family.


NOTE: the annotations are standardized before fitting the ridge regression model.  The resulting annotation ridge regression
coefficients are then un-standardized.  This approach was chosen so that the coefficients of all annotations are shrunk
by the ridge penalty on the same standardized scale.
"""

import json
from pathlib import Path

import narwhals
import numpy as np
import polars as pl
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
from mecfs_bio.build_system.task.dataframe_output import (
    ParquetOutFormat,
    write_df_according_to_format,
)
from mecfs_bio.build_system.task.ppp_database.allele_key import unordered_allele_key
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.polyfun_annotation_families import family_for_annotation

WEIGHTS_PARQUET_FILENAME = "weights.parquet"
DIAGNOSTICS_JSON_FILENAME = "diagnostics.json"
ANNOTATION_COL = "annotation"
GAMMA_RAW_COL = "gamma_raw"
GAMMA_STANDARDIZED_COL = "gamma_standardized"
FAMILY_COL = "family"
SNPVAR_COL = "snpvar_bin"
_CHR_COL = "CHR"
_BP_COL = "BP"
_A1_COL = "A1"
_A2_COL = "A2"
_ALLELE_KEY_COL = "allele_key"
# Both the annotation matrix and snpvar_meta carry alleles (A1/A2), so the
# annotation<->snpvar join is allele-precise: it matches on (CHR, BP,
# unordered-allele-key). This pairs each allele of a multiallelic site with its
# own snpvar_bin, rather than the old SNP-keyed join that (with a SNP dedup on
# each side) arbitrarily dropped one allele.
_JOIN_KEYS = [_CHR_COL, _BP_COL, _ALLELE_KEY_COL]

_DEFAULT_ALPHAS: tuple[float, ...] = (0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0)


@frozen
class _ChromStats:
    """Cross-product sufficient statistics for one chromosome (standardized-free,
    raw annotation scale). p = number of annotations.
    NOTE
    Compute
    X \\in \\mathbb{R}^{n \times p}

    sx= \\sum rows of annotation matrix X
    sxx \\sum x_i x_i^T
    sxy: X^T y
    sy: \\sum (y)\\in\\mathbb{R}
    syy: ||y||^2

    """

    n: int
    sx: np.ndarray  # (p,) sum of annotations
    sxx: np.ndarray  # (p, p) sum of a a^T
    sxy: np.ndarray  # (p,) sum of a * y
    sy: float
    syy: float


@frozen
class _StandardizedSystem:
    """The centered+standardized ridge system for a set of annotations."""

    g_std: np.ndarray  # (p, p) standardized Gram
    b_std: np.ndarray  # (p,) standardized cross-term
    mean: np.ndarray  # (p,) per-annotation mean
    sd: np.ndarray  # (p,) per-annotation std (zeros replaced by 1)
    mean_y: float


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
            .select(_CHR_COL, _BP_COL, _A1_COL, _A2_COL, SNPVAR_COL)
            .collect()
            .to_polars()
            .with_columns(unordered_allele_key(_A1_COL, _A2_COL).alias(_ALLELE_KEY_COL))
            .select(_CHR_COL, _BP_COL, _ALLELE_KEY_COL, SNPVAR_COL)
            .unique(subset=[_CHR_COL, _BP_COL, _ALLELE_KEY_COL])
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
                # Any parquet column not classifiable by family_for_annotation
                # hard-fails here by design: this is the de-facto guard that the
                # built parquet's columns match the known baseline-LF annotation set.
                FAMILY_COL: [family_for_annotation(c) for c in annot_columns],
            }
        )
        write_df_according_to_format(
            df=narwhals.from_native(weights).lazy(),
            out_path=scratch_dir / WEIGHTS_PARQUET_FILENAME,
            out_format=ParquetOutFormat(),
        )
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
        annot_chrom = (
            pl.scan_parquet(annot_path)
            .filter(pl.col(_CHR_COL) == chrom)
            .collect()
            .with_columns(unordered_allele_key(_A1_COL, _A2_COL).alias(_ALLELE_KEY_COL))
        )
        frame = annot_chrom.join(meta, on=_JOIN_KEYS, how="inner")
        # meta is unique on the join key and the annotation matrix is unique on
        # (CHR, BP, allele_key), so the inner join must not multiply rows.
        assert frame.height <= annot_chrom.height, (
            f"annotation<->snpvar join multiplied rows on chr{chrom}: "
            f"{annot_chrom.height} -> {frame.height}"
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


def _standardized_system(stats: _ChromStats) -> _StandardizedSystem:
    """Build the centered+standardized ridge system from raw cross-products.

    NOTE:
        to fit ridge regression we need
        X^TX: Gram matrix
            =sum x_i x_i^T where x_i is the ith row of x
        X^ty:
            = sum x_i^T y where x_i is the ith row of X


    Goal of this function is to convert from ridge data for unstandardized
     system to ridge data for standardized system

     Derivation of Centered Gram matrix:

     sum_i (x_i - mean)(x_i - mean)^T
          = sum_i x_ix_i^T  - x_i mean^T - mean x_i^T + mean mean^T
          =  (sum_i x_ix_i^T) - 2*n*mean mean^T + n * mean mean^T
          = (sum_i x_ix_i^T) - n * mean mean^T


    """
    n = stats.n
    mean = stats.sx / n
    var = np.diag(stats.sxx) / n - mean**2  # Since var(z)= Ez^2 - (E(z))^2
    sd = np.sqrt(np.maximum(var, 0.0))
    sd[sd == 0] = 1.0  # p-length vector
    centered_gram = stats.sxx - n * np.outer(mean, mean)  # sum_i (x_i-mean)(x_i-mean)^T
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


    NOTE:
        Definition of R^2
        R^2 = 1- (sum of squared residuals)/(sum of squared total)
        - sum of squared residuals = sum_i (y_i-f_i)^2
        - sum of squared total = sum_i (y_i-mean(y))^2


    Derivation:

        Let the held-out chromosome have n SNPs. For SNP i, x_i is its raw
        p-vector of annotations and y_i its snpvar_bin. The model was fit on
        the training chromosomes, so it carries the train intercept
        c = train_mean_y and the train standardizer (train_mean tm, train_sd
        ts, both p-vectors). The held-out annotations are standardized with the
        TRAIN moments, not their own:

            z_i = (x_i - tm) / ts       (elementwise)
            f_i = c + gamma_std . z_i

        Sum of squared residuals is then a quadratic form in gamma_std:

            SS_res = sum_i (y_i - f_i)^2
                   = sum_i [ (y_i - c) - gamma_std . z_i ]^2
                   = sum_i (y_i - c)^2                         (call it ss_res_y)
                     - 2 gamma_std . [ sum_i z_i (y_i - c) ]   (call it z_r)
                     + gamma_std^T [ sum_i z_i z_i^T ] gamma_std  (call it zz)

        Each of the three pieces expands into the raw cross-product statistics
        carried on held (n, sx, sxx, sxy, sy, syy), so no per-SNP array is ever
        rebuilt.

        ss_res_y (scalar):
            sum_i (y_i - c)^2 = sum_i y_i^2 - 2 c sum_i y_i + n c^2
                              = syy - 2 c sy + n c^2

        z_r (p-vector), via component j with z_ij = (x_ij - tm_j) / ts_j:
            sum_i z_ij (y_i - c)
              = (1 / ts_j) sum_i (x_ij - tm_j)(y_i - c)
              = (1 / ts_j) [ sxy_j - c sx_j - tm_j sy + n c tm_j ]
            => z_r = (sxy - c sx - tm sy + n c tm) / ts

        zz (pxp), via component (j, k):
            sum_i z_ij z_ik
              = (1 / (ts_j ts_k)) sum_i (x_ij - tm_j)(x_ik - tm_k)
              = (1 / (ts_j ts_k)) [ sxx_jk - tm_j sx_k - tm_k sx_j
                                    + n tm_j tm_k ]
            => zz = (sxx - outer(tm, sx) - outer(sx, tm) + n outer(tm, tm))
                    / outer(ts, ts)

        so SS_res = ss_res_y - 2 gamma_std . z_r + gamma_std^T zz gamma_std.

        The denominator uses the held-out chromosome's OWN mean ybar = sy / n
        (the R^2 null model is "predict the held-out mean"), unlike SS_res,
        which is taken about the train intercept c:

            SS_tot = sum_i (y_i - ybar)^2 = syy - n ybar^2

        and R^2 = 1 - SS_res / SS_tot.
    """
    n = held.n
    c = train_mean_y
    tm = train_mean
    ts = train_sd
    # SS over (y - c): sum (y - c)^2
    ss_res_y = held.syy - 2.0 * c * held.sy + n * c * c  # sum (y-c)^2
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
