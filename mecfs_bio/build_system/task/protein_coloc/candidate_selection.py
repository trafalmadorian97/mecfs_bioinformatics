"""
Candidate (trait locus, protein) pairs for colocalization, and their selection within a budget.

A protein is a candidate at a trait locus when one of its hit segments (from the protein hit
table) overlaps the locus window on the same chromosome. Excluded loci get no candidates.

Each candidate is classified by where the protein's gene lies (Sun et al. 2023 ST3, GRCh38):
- cis: the gene interval lies within cis_window_bp of the locus window;
- trans: the gene is elsewhere;
- unknown: the protein has no gene coordinates (for example a multi-gene complex).

Each pair gets a score, and pairs are taken in order of decreasing score:
    SCORE = SCORE_TRAIT + SCORE_PROTEIN - SCORE_PENALTY
    SCORE_TRAIT = min(trait lead -log10 p, clip)
    SCORE_PROTEIN = min(protein lead -log10 p, clip)
    SCORE_PENALTY = w * log2(1 + N_OTHER_PROTEINS) when the pair is not cis, else 0
N_OTHER_PROTEINS is the number of other proteins with a segment overlapping the same window. The
penalty demotes trans pairs at pleiotropic hotspots, where many proteins share a pQTL; cis pairs
are never penalized. Unknown is penalized like trans. Ties are broken by PROTEIN_ID, then
LOCUS_ID. A pair is selected while fewer than max_pairs are selected and its locus has fewer
than max_per_locus selected pairs.
"""

from collections.abc import Mapping

import polars as pl
from attrs import frozen

from mecfs_bio.constants.ppp_database_constants import Oid
from mecfs_bio.constants.protein_coloc_constants import (
    CIS_STATUS_CIS,
    CIS_STATUS_COL,
    CIS_STATUS_TRANS,
    CIS_STATUS_UNKNOWN,
    HIT_CHR_COL,
    HIT_GENE_COL,
    HIT_LEAD_MLOG10P_COL,
    HIT_SEG_END_COL,
    HIT_SEG_START_COL,
    LOCUS_CHR_COL,
    LOCUS_END_COL,
    LOCUS_ID_COL,
    LOCUS_LEAD_MLOG10P_COL,
    LOCUS_START_COL,
    N_OTHER_PROTEINS_COL,
    PROTEIN_ID_COL,
    PROTEIN_MLOG10P_COL,
    SCORE_COL,
    SCORE_PENALTY_COL,
    SCORE_PROTEIN_COL,
    SCORE_TRAIT_COL,
    STATUS_CANDIDATE,
    STATUS_COL,
    STATUS_NOT_SELECTED_BUDGET,
    STATUS_SELECTED,
)
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

_GENE_CHR = "_gene_chr"
_GENE_START = "_gene_start"
_GENE_END = "_gene_end"

CANDIDATE_SCHEMA = pl.Schema(
    {
        LOCUS_ID_COL: pl.String,
        PROTEIN_ID_COL: pl.String,
        HIT_GENE_COL: pl.String,
        PROTEIN_MLOG10P_COL: pl.Float64,
        CIS_STATUS_COL: pl.String,
        N_OTHER_PROTEINS_COL: pl.Int64,
    }
)


@frozen(slots=True)
class SelectionConfig:
    max_pairs: int = 30
    max_per_locus: int = 5
    mlog10p_clip: float = 50.0
    trans_penalty_weight: float = 2.0
    cis_window_bp: int = 1_000_000

    def __attrs_post_init__(self) -> None:
        assert self.max_pairs >= 1 and self.max_per_locus >= 1
        assert self.mlog10p_clip > 0 and self.trans_penalty_weight >= 0
        assert self.cis_window_bp >= 0


def _gene_coords_frame(gene_coords: Mapping[Oid, GenomicInterval]) -> pl.DataFrame:
    return pl.DataFrame(
        {
            PROTEIN_ID_COL: [str(oid) for oid in gene_coords],
            _GENE_CHR: [interval.chrom for interval in gene_coords.values()],
            _GENE_START: [interval.start for interval in gene_coords.values()],
            _GENE_END: [interval.end for interval in gene_coords.values()],
        },
        schema={
            PROTEIN_ID_COL: pl.String,
            _GENE_CHR: pl.Int64,
            _GENE_START: pl.Int64,
            _GENE_END: pl.Int64,
        },
    )


def _cis_status_expr(cis_window_bp: int) -> pl.Expr:
    near = (
        (pl.col(_GENE_CHR) == pl.col(LOCUS_CHR_COL))
        & (pl.col(_GENE_START) <= pl.col(LOCUS_END_COL) + cis_window_bp)
        & (pl.col(_GENE_END) >= pl.col(LOCUS_START_COL) - cis_window_bp)
    )
    return (
        pl.when(pl.col(_GENE_CHR).is_null())
        .then(pl.lit(CIS_STATUS_UNKNOWN))
        .when(near)
        .then(pl.lit(CIS_STATUS_CIS))
        .otherwise(pl.lit(CIS_STATUS_TRANS))
        .alias(CIS_STATUS_COL)
    )


def build_candidates(
    loci: pl.DataFrame,
    hit_segments: pl.DataFrame,
    gene_coords: Mapping[Oid, GenomicInterval],
    cis_window_bp: int,
) -> pl.DataFrame:
    """One row per (candidate locus, protein with a segment overlapping its window)."""
    windows = loci.filter(pl.col(STATUS_COL) == STATUS_CANDIDATE).select(
        LOCUS_ID_COL,
        pl.col(LOCUS_CHR_COL).cast(pl.Int64),
        LOCUS_START_COL,
        LOCUS_END_COL,
    )
    segments = hit_segments.select(
        PROTEIN_ID_COL,
        HIT_GENE_COL,
        pl.col(HIT_CHR_COL).cast(pl.Int64).alias(LOCUS_CHR_COL),
        HIT_SEG_START_COL,
        HIT_SEG_END_COL,
        HIT_LEAD_MLOG10P_COL,
    )
    overlapping = windows.join(segments, on=LOCUS_CHR_COL, how="inner").filter(
        pl.col(HIT_SEG_START_COL) <= pl.col(LOCUS_END_COL),
        pl.col(HIT_SEG_END_COL) >= pl.col(LOCUS_START_COL),
    )
    pairs = overlapping.group_by(
        LOCUS_ID_COL, LOCUS_CHR_COL, LOCUS_START_COL, LOCUS_END_COL, PROTEIN_ID_COL
    ).agg(
        pl.col(HIT_GENE_COL).first(),
        pl.col(HIT_LEAD_MLOG10P_COL).max().alias(PROTEIN_MLOG10P_COL),
    )
    return (
        pairs.with_columns(
            (pl.len().over(LOCUS_ID_COL) - 1).cast(pl.Int64).alias(N_OTHER_PROTEINS_COL)
        )
        .join(_gene_coords_frame(gene_coords), on=PROTEIN_ID_COL, how="left")
        .with_columns(_cis_status_expr(cis_window_bp))
        .select(CANDIDATE_SCHEMA.names())
        .cast(CANDIDATE_SCHEMA)
        .sort(LOCUS_ID_COL, PROTEIN_ID_COL)
    )


def _score_columns(config: SelectionConfig) -> list[pl.Expr]:
    penalty = (
        pl.when(pl.col(CIS_STATUS_COL) == CIS_STATUS_CIS)
        .then(pl.lit(0.0))
        .otherwise(
            config.trans_penalty_weight * (1 + pl.col(N_OTHER_PROTEINS_COL)).log(2)
        )
    )
    return [
        pl.min_horizontal(pl.col(LOCUS_LEAD_MLOG10P_COL), config.mlog10p_clip).alias(
            SCORE_TRAIT_COL
        ),
        pl.min_horizontal(pl.col(PROTEIN_MLOG10P_COL), config.mlog10p_clip).alias(
            SCORE_PROTEIN_COL
        ),
        penalty.alias(SCORE_PENALTY_COL),
    ]


def _budget_statuses(ranked_locus_ids: list[str], config: SelectionConfig) -> list[str]:
    """Walk the ranking, selecting while the pair budget and the locus cap allow."""
    per_locus: dict[str, int] = {}
    n_selected = 0
    statuses: list[str] = []
    for locus_id in ranked_locus_ids:
        taken = per_locus.get(locus_id, 0)
        if n_selected < config.max_pairs and taken < config.max_per_locus:
            per_locus[locus_id] = taken + 1
            n_selected += 1
            statuses.append(STATUS_SELECTED)
        else:
            statuses.append(STATUS_NOT_SELECTED_BUDGET)
    return statuses


def select_candidates(
    candidates: pl.DataFrame, loci: pl.DataFrame, config: SelectionConfig
) -> pl.DataFrame:
    """Score the candidates and mark each selected or not_selected_budget, in ranked order."""
    ranked = (
        candidates.join(
            loci.select(LOCUS_ID_COL, LOCUS_LEAD_MLOG10P_COL),
            on=LOCUS_ID_COL,
            how="left",
        )
        .with_columns(_score_columns(config))
        .with_columns(
            (
                pl.col(SCORE_TRAIT_COL)
                + pl.col(SCORE_PROTEIN_COL)
                - pl.col(SCORE_PENALTY_COL)
            ).alias(SCORE_COL)
        )
        .drop(LOCUS_LEAD_MLOG10P_COL)
        .sort(
            [SCORE_COL, PROTEIN_ID_COL, LOCUS_ID_COL],
            descending=[True, False, False],
        )
    )
    return ranked.with_columns(
        pl.Series(
            STATUS_COL,
            _budget_statuses(ranked[LOCUS_ID_COL].to_list(), config),
            dtype=pl.String,
        )
    )
