"""
Genome-wide significant loci of a trait, as windows around lead variants.

Greedy distance clumping: take the most significant remaining variant below the p threshold as
a lead, open the window lead +- half_window_bp, remove every significant variant inside it, and
repeat. Windows on one chromosome that then overlap are merged; a merged locus keeps the most
significant of its leads. This needs no LD.

A locus whose window overlaps any excluded interval (by default the extended MHC, whose long
range LD breaks the single-causal-variant assumption) gets status excluded_mhc and is reported
but not tested.

p is the two-sided normal p-value of BETA / SE, kept in log space.
"""

from collections.abc import Sequence

import numpy as np
import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    two_sided_mlog10p,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_MLOG10P_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    LOCUS_CHR_COL,
    LOCUS_END_COL,
    LOCUS_ID_COL,
    LOCUS_LEAD_EA_COL,
    LOCUS_LEAD_MLOG10P_COL,
    LOCUS_LEAD_NEA_COL,
    LOCUS_LEAD_POS_COL,
    LOCUS_START_COL,
    STATUS_CANDIDATE,
    STATUS_COL,
    STATUS_EXCLUDED_MHC,
)
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

LOCI_SCHEMA = pl.Schema(
    {
        LOCUS_ID_COL: pl.String,
        LOCUS_CHR_COL: pl.Int64,
        LOCUS_START_COL: pl.Int64,
        LOCUS_END_COL: pl.Int64,
        LOCUS_LEAD_POS_COL: pl.Int64,
        LOCUS_LEAD_EA_COL: pl.String,
        LOCUS_LEAD_NEA_COL: pl.String,
        LOCUS_LEAD_MLOG10P_COL: pl.Float64,
        STATUS_COL: pl.String,
    }
)


@frozen(slots=True)
class _Lead:
    chrom: int
    pos: int
    ea: str
    nea: str
    mlog10p: float
    start: int
    end: int


def _significant(trait: pl.DataFrame, p_threshold: float) -> pl.DataFrame:
    with np.errstate(invalid="ignore", divide="ignore"):
        mlog10p = two_sided_mlog10p(
            (trait[GWASLAB_BETA_COL] / trait[GWASLAB_SE_COL]).to_numpy()
        )
    return (
        trait.select(
            GWASLAB_CHROM_COL,
            GWASLAB_POS_COL,
            GWASLAB_EFFECT_ALLELE_COL,
            GWASLAB_NON_EFFECT_ALLELE_COL,
        )
        .with_columns(pl.Series(GWASLAB_MLOG10P_COL, mlog10p))
        .filter(
            pl.col(GWASLAB_MLOG10P_COL).is_finite(),
            pl.col(GWASLAB_MLOG10P_COL) > -np.log10(p_threshold),
        )
        .sort(
            [GWASLAB_MLOG10P_COL, GWASLAB_CHROM_COL, GWASLAB_POS_COL],
            descending=[True, False, False],
        )
    )


def _greedy_leads(significant: pl.DataFrame, half_window_bp: int) -> list[_Lead]:
    """Leads in decreasing significance; a variant within an earlier lead's window is skipped."""
    leads: list[_Lead] = []
    for chrom, pos, ea, nea, mlog10p in significant.iter_rows():
        if any(
            lead.chrom == chrom and abs(lead.pos - pos) <= half_window_bp
            for lead in leads
        ):
            continue
        leads.append(
            _Lead(
                chrom=chrom,
                pos=pos,
                ea=ea,
                nea=nea,
                mlog10p=mlog10p,
                start=max(0, pos - half_window_bp),
                end=pos + half_window_bp,
            )
        )
    return leads


def _merge_overlapping(leads: list[_Lead]) -> list[_Lead]:
    merged: list[_Lead] = []
    for lead in sorted(leads, key=lambda x: (x.chrom, x.start)):
        previous = merged[-1] if merged else None
        if (
            previous is None
            or previous.chrom != lead.chrom
            or lead.start > previous.end
        ):
            merged.append(lead)
            continue
        best = previous if previous.mlog10p >= lead.mlog10p else lead
        merged[-1] = _Lead(
            chrom=best.chrom,
            pos=best.pos,
            ea=best.ea,
            nea=best.nea,
            mlog10p=best.mlog10p,
            start=previous.start,
            end=max(previous.end, lead.end),
        )
    return merged


def _overlaps_any(lead: _Lead, exclude: Sequence[GenomicInterval]) -> bool:
    return any(
        lead.chrom == interval.chrom
        and lead.start <= interval.end
        and lead.end >= interval.start
        for interval in exclude
    )


def find_trait_loci(
    trait: pl.DataFrame,
    p_threshold: float,
    half_window_bp: int,
    exclude: Sequence[GenomicInterval],
) -> pl.DataFrame:
    """One row per merged locus window, in LOCI_SCHEMA, sorted by chromosome and start."""
    loci = _merge_overlapping(
        _greedy_leads(_significant(trait, p_threshold), half_window_bp)
    )
    return pl.DataFrame(
        [
            {
                LOCUS_ID_COL: f"chr{lead.chrom}:{lead.start}-{lead.end}",
                LOCUS_CHR_COL: lead.chrom,
                LOCUS_START_COL: lead.start,
                LOCUS_END_COL: lead.end,
                LOCUS_LEAD_POS_COL: lead.pos,
                LOCUS_LEAD_EA_COL: lead.ea,
                LOCUS_LEAD_NEA_COL: lead.nea,
                LOCUS_LEAD_MLOG10P_COL: lead.mlog10p,
                STATUS_COL: (
                    STATUS_EXCLUDED_MHC
                    if _overlaps_any(lead, exclude)
                    else STATUS_CANDIDATE
                ),
            }
            for lead in loci
        ],
        schema=LOCI_SCHEMA,
    )
