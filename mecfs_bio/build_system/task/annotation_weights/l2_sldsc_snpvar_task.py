"""Trait-specific per-SNP heritability (snpvar) prior via L2-regularized S-LDSC.

PolyFun's Approach 2: estimate per-annotation heritability coefficients tau from
the GWAS being fine-mapped, then score every reference variant as
snpvar_i = sum_c a_ic tau_c over the baseline-LF annotations a. The output is a
drop-in replacement for the precomputed (Approach 1) prior at the SUSIE PriorInfo
seam; SUSIE normalizes the prior within each locus, so only relative snpvar
matters.

Model. Stratified LD score regression says

    E[chi2_i] = 1 + N sum_c tau_c ell(i, c)

where ell(i, c) is variant i's LD-score for annotation c. We regress chi2_i on the
187 annotation LD-scores with weighted ridge (free intercept), then
tau = beta_raw / N. The LDSC regression weights are

    omega_i = 1 / (het_i * oc_i),  het_i = (1 + N h2bar l_i / M)^2,  oc_i = max(l_i, 1)

with l_i the variant's total LD-score and M the reference variant count. The
distributed LD-score members have no all-ones base column, but the 20 MAFbin_*
annotations partition every variant, so l_i is the sum of the 20 MAFbin_*
LD-scores. h2bar is estimated from the univariate LDSC slope of chi2 on l
(slope = N h2bar / M).

No leakage. Chromosomes are split by parity. Everything that scores the even
chromosomes (h2bar, the ridge penalty chosen by leave-one-chromosome-out
cross-validation, and tau) is fit on the odd chromosomes only, and vice versa.

Memory. Nothing genome-wide and dense is materialized: each pass streams one
chromosome's LD-score member (or annotation-matrix slice) at a time, and the
regression is carried by per-chromosome sufficient statistics (see
chromosome_blocked_ridge).

Outputs: snpvar.parquet (CHR, POS, NEA, EA, snpvar), one weights table per half
holding its tau in the RidgeAnnotationWeightsTask schema, and diagnostics.json.
Since snpvar is exactly linear in the annotations, the weights table of the tau
that scored a chromosome explains that chromosome's snpvar exactly.

All per-variant joins use the exact (CHR, POS, NEA, EA) key. The PolyFun inputs are
REF-oriented (A1 == REF == NEA, A2 == ALT == EA), so the sumstats must be
genome-reference harmonized (NEA == hg19 REF).
"""

import json
import re
from collections.abc import Mapping, Sequence
from pathlib import Path, PurePath

import narwhals as nw
import numpy as np
import polars as pl
import pyarrow.parquet as pq
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.annotation_weights.chromosome_blocked_ridge import (
    AlphaSelection,
    ChromRidgeBlock,
    accumulate_block,
    combine,
    fit,
    select_alpha_loco,
)
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    ANNOTATION_COL,
    FAMILY_COL,
    GAMMA_RAW_COL,
    GAMMA_STANDARDIZED_COL,
)
from mecfs_bio.build_system.task.annotation_weights.stream_extract_annotation_parquets_task import (
    LDSCORE_PARQUET_MEMBER_RE,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
    GWASLAB_Z_COL,
)
from mecfs_bio.constants.polyfun_annotation_families import family_for_annotation
from mecfs_bio.constants.polyfun_constants import (
    POLYFUN_A1_COL,
    POLYFUN_A2_COL,
    POLYFUN_BP_COL,
    POLYFUN_CHR_COL,
    POLYFUN_SNP_COL,
    POLYFUN_TO_GWASLAB_KEY_RENAME,
)

logger = structlog.get_logger()

SNPVAR_PARQUET_FILENAME = "snpvar.parquet"
DIAGNOSTICS_JSON_FILENAME = "diagnostics.json"
SNPVAR_COL = "snpvar"
# Per-half tau in the RidgeAnnotationWeightsTask weights schema (annotation,
# gamma_raw, gamma_standardized, family), so the polyfun explainability contrast
# can consume them in place of the ridge surrogate of the precomputed prior.
# tau fit on odd chromosomes scores the even chromosomes, and vice versa.
TAU_ODD_WEIGHTS_FILENAME = "weights_tau_odd.parquet"
TAU_EVEN_WEIGHTS_FILENAME = "weights_tau_even.parquet"
MAFBIN_LDSCORE_RE = re.compile(r"^MAFbin_(lowfreq|frequent)_\d+$")
_N_MAFBINS = 20
# chi2_i = Z_i^2, the regression target.
_CHI2_COL = "chi2"
# l_i = sum over the 20 MAFbin_* LD-score columns = total LD-score of variant i.
_TOTAL_LDSCORE_COL = "l"
# PolyFun's MAX_CHI2: variants above this are dropped from the regression.
_CHI2_CAP = 80.0
_H2BAR_FLOOR = 1e-8
_POLYFUN_KEY_COLS = [
    POLYFUN_CHR_COL,
    POLYFUN_BP_COL,
    POLYFUN_SNP_COL,
    POLYFUN_A1_COL,
    POLYFUN_A2_COL,
]
_JOIN_KEYS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_COL,
]
# Half-decade steps from 1e-1 to 1e12. The per-SNP chi2 regression is noise
# dominated, so cross-validation can favor very heavy shrinkage; if it selects the
# largest value, widen this grid and refit. Scoring extra penalties is cheap: each
# solves one p x p system per cross-validation fold.
_DEFAULT_ALPHAS: tuple[float, ...] = tuple(10.0 ** (k / 2) for k in range(-2, 25))


@frozen
class _ParityFit:
    """Everything fit on one parity's chromosomes; it scores the other parity."""

    h2bar: float
    selection: AlphaSelection
    tau: np.ndarray
    # tau on the standardized-LD-score scale (beta_std / N), comparable across
    # annotations.
    tau_standardized: np.ndarray
    n_regression_variants: int


@frozen
class L2RegularizedSldscSnpvarTask(Task):
    meta: Meta
    sumstats_task: Task
    annotation_ldscore_members_task: Task
    annotation_matrix_task: Task
    effective_sample_size: float
    alphas: tuple[float, ...] = _DEFAULT_ALPHAS

    def __attrs_post_init__(self) -> None:
        assert self.effective_sample_size > 0, "effective_sample_size must be positive"
        assert len(self.alphas) > 0, "need at least one candidate ridge penalty"

    @property
    def deps(self) -> list[Task]:
        return [
            self.sumstats_task,
            self.annotation_ldscore_members_task,
            self.annotation_matrix_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> DirectoryAsset:
        ldscore_asset = fetch(self.annotation_ldscore_members_task.asset_id)
        assert isinstance(ldscore_asset, DirectoryAsset)
        annot_asset = fetch(self.annotation_matrix_task.asset_id)
        assert isinstance(annot_asset, FileAsset)
        sumstats = _load_sumstats(
            scan_dataframe_asset(
                fetch(self.sumstats_task.asset_id), self.sumstats_task.meta
            )
        )

        ldscore_paths = _ldscore_member_paths(ldscore_asset.path)
        annot_cols = _annotation_columns(annot_asset.path)
        _assert_ldscore_columns_match(ldscore_paths, annot_cols)
        mafbins = [c for c in annot_cols if _is_mafbin(c)]
        assert len(mafbins) == _N_MAFBINS, (
            f"expected {_N_MAFBINS} MAFbin columns, got {len(mafbins)}"
        )
        m_ref = _reference_variant_count(ldscore_paths)
        n_bar = self.effective_sample_size

        # Pass A: per-chromosome univariate (l, chi2) blocks for per-parity h2bar.
        univariate: dict[int, ChromRidgeBlock] = {}
        for chrom, path in sorted(ldscore_paths.items()):
            frame = _ldscore_regression_frame(
                ldscore_path=path,
                sumstats=sumstats.filter(pl.col(GWASLAB_CHROM_COL) == chrom),
                ldscore_cols=mafbins,
            )
            univariate[chrom] = accumulate_block(
                x=frame.select(_TOTAL_LDSCORE_COL).to_numpy(),
                y=frame[_CHI2_COL].to_numpy(),
            )
        h2bar = {
            parity: _h2bar(
                combine(list(_of_parity(univariate, parity).values())),
                m_ref=m_ref,
                n_bar=n_bar,
            )
            for parity in (0, 1)
        }

        # Pass B: per-chromosome weighted blocks over the annotation LD-scores,
        # each chromosome weighted with its own parity's h2bar.
        weighted: dict[int, ChromRidgeBlock] = {}
        for chrom, path in sorted(ldscore_paths.items()):
            frame = _ldscore_regression_frame(
                ldscore_path=path,
                sumstats=sumstats.filter(pl.col(GWASLAB_CHROM_COL) == chrom),
                ldscore_cols=annot_cols,
            )
            total_ld = frame[_TOTAL_LDSCORE_COL].to_numpy()
            weighted[chrom] = accumulate_block(
                x=frame.select(annot_cols).to_numpy(),
                y=frame[_CHI2_COL].to_numpy(),
                w=_ldsc_weights(
                    total_ld=total_ld,
                    h2bar=h2bar[chrom % 2],
                    m_ref=m_ref,
                    n_bar=n_bar,
                ),
            )
            logger.info(
                "accumulated weighted S-LDSC block",
                chromosome=chrom,
                n_regression_variants=frame.height,
            )

        fits = {
            parity: _fit_parity(
                _of_parity(weighted, parity),
                h2bar=h2bar[parity],
                # Pass A and pass B filter identically, so the unweighted pass-A
                # blocks count the regression variants.
                n_regression_variants=sum(
                    _block_count(b) for b in _of_parity(univariate, parity).values()
                ),
                alphas=self.alphas,
                n_bar=n_bar,
            )
            for parity in (0, 1)
        }

        # Pass C: score every annotation-matrix variant with the opposite parity's tau.
        n_scored = _write_snpvar(
            annot_path=annot_asset.path,
            annot_cols=annot_cols,
            tau_by_parity={parity: fits[1 - parity].tau for parity in (0, 1)},
            out_path=scratch_dir / SNPVAR_PARQUET_FILENAME,
        )
        for parity, filename in (
            (1, TAU_ODD_WEIGHTS_FILENAME),
            (0, TAU_EVEN_WEIGHTS_FILENAME),
        ):
            _tau_weights_table(fits[parity], annot_cols).write_parquet(
                scratch_dir / filename
            )
        (scratch_dir / DIAGNOSTICS_JSON_FILENAME).write_text(
            json.dumps(
                _diagnostics(
                    fits=fits,
                    annot_cols=annot_cols,
                    n_bar=n_bar,
                    m_ref=m_ref,
                    n_scored=n_scored,
                ),
                indent=2,
                sort_keys=True,
            )
        )
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        *,
        sumstats_task: Task,
        effective_sample_size: float,
        annotation_ldscore_members_task: Task,
        annotation_matrix_task: Task,
        alphas: tuple[float, ...] = _DEFAULT_ALPHAS,
    ) -> "L2RegularizedSldscSnpvarTask":
        source_meta = sumstats_task.meta
        assert isinstance(source_meta, (FilteredGWASDataMeta, ResultTableMeta)), (
            f"sumstats task must carry trait/project metadata, got {type(source_meta)}"
        )
        meta = ResultDirectoryMeta(
            id=asset_id,
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=PurePath("analysis") / "polyfun_l2_sldsc_prior",
        )
        return cls(
            meta=meta,
            sumstats_task=sumstats_task,
            annotation_ldscore_members_task=annotation_ldscore_members_task,
            annotation_matrix_task=annotation_matrix_task,
            effective_sample_size=effective_sample_size,
            alphas=alphas,
        )


def _is_mafbin(col: str) -> bool:
    return MAFBIN_LDSCORE_RE.match(col) is not None


def _of_parity(
    blocks: Mapping[int, ChromRidgeBlock], parity: int
) -> dict[int, ChromRidgeBlock]:
    return {c: b for c, b in blocks.items() if c % 2 == parity}


def _block_count(block: ChromRidgeBlock) -> int:
    """Number of variants in an unweighted block (its total weight)."""
    return round(block.sw)


def _cast_key(frame: pl.DataFrame) -> pl.DataFrame:
    return frame.with_columns(
        pl.col(GWASLAB_CHROM_COL).cast(pl.Int64),
        pl.col(GWASLAB_POS_COL).cast(pl.Int64),
        pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).cast(pl.String),
        pl.col(GWASLAB_EFFECT_ALLELE_COL).cast(pl.String),
    )


def _load_sumstats(frame: nw.LazyFrame) -> pl.DataFrame:
    """The four-part key and Z, keeping only variants with a finite Z and a key that
    occurs exactly once (a duplicated key cannot be matched unambiguously).

    Z is resolved in gwaslab's order: an existing Z column, otherwise BETA / SE.
    A non-finite Z (e.g. SE == 0) is dropped.
    """
    names = frame.collect_schema().names()
    missing_key = [c for c in _JOIN_KEYS if c not in names]
    assert not missing_key, f"sumstats lack join-key columns {missing_key}"
    if GWASLAB_Z_COL in names:
        z = nw.col(GWASLAB_Z_COL)
    else:
        assert GWASLAB_BETA_COL in names and GWASLAB_SE_COL in names, (
            f"sumstats need {GWASLAB_Z_COL}, or {GWASLAB_BETA_COL} and "
            f"{GWASLAB_SE_COL}; found {names}"
        )
        z = nw.col(GWASLAB_BETA_COL) / nw.col(GWASLAB_SE_COL)
    loaded = _cast_key(
        frame.select(*_JOIN_KEYS, z.alias(GWASLAB_Z_COL)).collect().to_polars()
    )
    finite = loaded.filter(pl.col(GWASLAB_Z_COL).is_finite())
    unique = finite.unique(subset=_JOIN_KEYS, keep="none")
    logger.info(
        "loaded sumstats for S-LDSC",
        n_loaded=loaded.height,
        n_dropped_non_finite_z=loaded.height - finite.height,
        n_dropped_duplicate_key=finite.height - unique.height,
    )
    return unique


def _ldscore_member_paths(directory: Path) -> dict[int, Path]:
    paths: dict[int, Path] = {}
    for path in directory.iterdir():
        match = LDSCORE_PARQUET_MEMBER_RE.search(path.name)
        if match is not None:
            paths[int(match.group(1))] = path
    assert paths, f"no LD-score members found in {directory}"
    return paths


def _annotation_columns(annot_path: Path) -> list[str]:
    names = pl.scan_parquet(annot_path).collect_schema().names()
    missing = [c for c in _POLYFUN_KEY_COLS if c not in names]
    assert not missing, f"annotation matrix lacks key columns {missing}"
    return [c for c in names if c not in _POLYFUN_KEY_COLS]


def _assert_ldscore_columns_match(
    ldscore_paths: Mapping[int, Path], annot_cols: Sequence[str]
) -> None:
    """tau is fit on LD-score columns and applied to annotation columns by name, so
    every member must carry the four-part key and exactly the annotation columns."""
    for chrom, path in ldscore_paths.items():
        names = pl.scan_parquet(path).collect_schema().names()
        missing = [c for c in _POLYFUN_KEY_COLS if c not in names]
        assert not missing, f"chr{chrom} LD-score member lacks key columns {missing}"
        ld_cols = [c for c in names if c not in _POLYFUN_KEY_COLS]
        assert set(ld_cols) == set(annot_cols), (
            f"chr{chrom} LD-score columns differ from annotation columns: "
            f"only in LD-scores {sorted(set(ld_cols) - set(annot_cols))}, "
            f"only in annotations {sorted(set(annot_cols) - set(ld_cols))}"
        )


def _reference_variant_count(ldscore_paths: Mapping[int, Path]) -> int:
    return sum(
        pl.scan_parquet(path).select(pl.len()).collect().item()
        for path in ldscore_paths.values()
    )


def _read_ldscore_member(path: Path, columns: Sequence[str]) -> pl.DataFrame:
    """Read the given columns of one LD-score member with buffered reads.

    The members (~29GB) are candidates for an asset-store remap onto a network-like
    mount (e.g. a WSL2 DrvFs drive), where polars' memory-mapped parquet reader
    faults the file in one page at a time and runs orders of magnitude slower than
    large sequential reads. pyarrow with memory_map=False reads whole column chunks.
    """
    table = pq.read_table(path, columns=list(columns), memory_map=False)
    frame = pl.from_arrow(table)
    assert isinstance(frame, pl.DataFrame)
    return frame


def _ldscore_regression_frame(
    ldscore_path: Path, sumstats: pl.DataFrame, ldscore_cols: Sequence[str]
) -> pl.DataFrame:
    """One chromosome's LD-scores joined to its sumstats on the four-part key.

    Reads the key plus ldscore_cols (which must include the MAFbin columns) from
    the member, renaming its REF-oriented A1/A2 to NEA/EA. Adds chi2 (capped at
    _CHI2_CAP) and the MAFbin-sum total LD-score l.
    """
    ld = _cast_key(
        _read_ldscore_member(
            ldscore_path,
            columns=[
                POLYFUN_CHR_COL,
                POLYFUN_BP_COL,
                POLYFUN_A1_COL,
                POLYFUN_A2_COL,
                *ldscore_cols,
            ],
        ).rename(POLYFUN_TO_GWASLAB_KEY_RENAME)
    )
    mafbins = [c for c in ldscore_cols if _is_mafbin(c)]
    assert len(mafbins) == _N_MAFBINS, (
        f"expected {_N_MAFBINS} MAFbin columns, got {len(mafbins)}"
    )
    ld = ld.with_columns(pl.sum_horizontal(mafbins).alias(_TOTAL_LDSCORE_COL))
    frame = ld.join(sumstats, on=_JOIN_KEYS, how="inner")
    # sumstats are unique on the key, so the join cannot multiply LD-score rows.
    assert frame.height <= ld.height, "LD-score/sumstats join multiplied rows"
    return frame.with_columns((pl.col(GWASLAB_Z_COL) ** 2).alias(_CHI2_COL)).filter(
        pl.col(_CHI2_COL) < _CHI2_CAP
    )


def _h2bar(univariate: ChromRidgeBlock, m_ref: int, n_bar: float) -> float:
    """Mean per-variant heritability scale from the LDSC slope of chi2 on l:
    E[chi2] = 1 + (N h2bar / M) l, so h2bar = slope M / N (floored above zero)."""
    slope = float(fit(univariate, alpha=0.0).beta_raw[0])
    return max(slope * m_ref / n_bar, _H2BAR_FLOOR)


def _ldsc_weights(
    total_ld: np.ndarray, h2bar: float, m_ref: int, n_bar: float
) -> np.ndarray:
    """omega_i = 1 / (het_i * oc_i), het_i = (1 + N h2bar l_i / M)^2, oc_i = max(l_i, 1)."""
    het = (1.0 + n_bar * h2bar * total_ld / m_ref) ** 2
    oc = np.maximum(total_ld, 1.0)
    return 1.0 / (het * oc)


def _fit_parity(
    blocks: Mapping[int, ChromRidgeBlock],
    h2bar: float,
    n_regression_variants: int,
    alphas: Sequence[float],
    n_bar: float,
) -> _ParityFit:
    """Choose alpha by leave-one-chromosome-out CV within this parity's blocks only,
    then refit on all of them; tau = beta_raw / N."""
    selection = select_alpha_loco(blocks, alphas)
    ridge = fit(combine(list(blocks.values())), alpha=selection.alpha)
    return _ParityFit(
        h2bar=h2bar,
        selection=selection,
        tau=ridge.beta_raw / n_bar,
        tau_standardized=ridge.beta_std / n_bar,
        n_regression_variants=n_regression_variants,
    )


def _tau_weights_table(
    parity_fit: _ParityFit, annot_cols: Sequence[str]
) -> pl.DataFrame:
    """One half's tau as a ridge-weights table. family_for_annotation hard-fails on
    an annotation outside the known baseline-LF set."""
    return pl.DataFrame(
        {
            ANNOTATION_COL: list(annot_cols),
            GAMMA_RAW_COL: parity_fit.tau,
            GAMMA_STANDARDIZED_COL: parity_fit.tau_standardized,
            FAMILY_COL: [family_for_annotation(c) for c in annot_cols],
        }
    )


def _write_snpvar(
    annot_path: Path,
    annot_cols: Sequence[str],
    tau_by_parity: Mapping[int, np.ndarray],
    out_path: Path,
) -> int:
    """Stream the annotation matrix one chromosome at a time, scoring
    snpvar_i = a_i . tau_by_parity[chrom % 2], and write (CHR, POS, NEA, EA, snpvar)."""
    chroms = (
        pl.scan_parquet(annot_path)
        .select(POLYFUN_CHR_COL)
        .unique()
        .collect()
        .to_series()
        .to_list()
    )
    n_scored = 0
    writer: pq.ParquetWriter | None = None
    try:
        for chrom in sorted(chroms):
            annot = (
                pl.scan_parquet(annot_path)
                .filter(pl.col(POLYFUN_CHR_COL) == chrom)
                .collect()
            )
            snpvar = annot.select(annot_cols).to_numpy() @ tau_by_parity[chrom % 2]
            scored = _cast_key(
                annot.select(
                    POLYFUN_CHR_COL, POLYFUN_BP_COL, POLYFUN_A1_COL, POLYFUN_A2_COL
                ).rename(POLYFUN_TO_GWASLAB_KEY_RENAME)
            ).with_columns(pl.Series(SNPVAR_COL, snpvar))
            table = scored.to_arrow()
            if writer is None:
                writer = pq.ParquetWriter(out_path, table.schema)
            writer.write_table(table)
            n_scored += scored.height
    finally:
        if writer is not None:
            writer.close()
    assert n_scored > 0, f"annotation matrix {annot_path} had no variants to score"
    return n_scored


def _diagnostics(
    fits: Mapping[int, _ParityFit],
    annot_cols: Sequence[str],
    n_bar: float,
    m_ref: int,
    n_scored: int,
) -> dict[str, object]:
    out: dict[str, object] = {
        "n_bar": n_bar,
        "m_ref": m_ref,
        "n_scored_variants": n_scored,
    }
    for parity, name in ((1, "odd"), (0, "even")):
        f = fits[parity]
        out[f"h2bar_{name}"] = f.h2bar
        out[f"alpha_{name}"] = f.selection.alpha
        out[f"mean_heldout_r2_{name}"] = f.selection.mean_r2
        out[f"heldout_r2_per_chrom_{name}"] = {
            str(c): r2 for c, r2 in sorted(f.selection.r2_per_chrom.items())
        }
        out[f"n_regression_variants_{name}"] = f.n_regression_variants
        out[f"tau_{name}"] = dict(zip(annot_cols, f.tau.tolist()))
    return out
