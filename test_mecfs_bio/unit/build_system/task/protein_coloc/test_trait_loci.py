import polars as pl

from mecfs_bio.build_system.task.protein_coloc.trait_loci import find_trait_loci
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    LOCUS_END_COL,
    LOCUS_ID_COL,
    LOCUS_LEAD_EA_COL,
    LOCUS_LEAD_POS_COL,
    LOCUS_START_COL,
    STATUS_CANDIDATE,
    STATUS_COL,
    STATUS_EXCLUDED_MHC,
)
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

_SE = 0.01


def _trait(rows: list[tuple[int, int, float]]) -> pl.DataFrame:
    """Rows of (CHR, POS, z)."""
    return pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [r[0] for r in rows],
            GWASLAB_POS_COL: [r[1] for r in rows],
            GWASLAB_EFFECT_ALLELE_COL: ["A"] * len(rows),
            GWASLAB_NON_EFFECT_ALLELE_COL: ["G"] * len(rows),
            GWASLAB_BETA_COL: [r[2] * _SE for r in rows],
            GWASLAB_SE_COL: [_SE] * len(rows),
        }
    )


def test_nearby_hits_merge_into_one_locus() -> None:
    # 1.6 Mb is outside the first lead's +-500 kb, so it leads its own window, which then
    # overlaps the first and merges into it; the null variant never leads.
    trait = _trait([(1, 1_000_000, 12.0), (1, 1_600_000, 8.0), (1, 3_000_000, 1.0)])
    loci = find_trait_loci(trait, p_threshold=5e-8, half_window_bp=500_000, exclude=[])
    assert loci.height == 1
    row = loci.row(0, named=True)
    assert row[LOCUS_LEAD_POS_COL] == 1_000_000
    assert row[LOCUS_START_COL] == 500_000
    assert row[LOCUS_END_COL] == 2_100_000
    assert row[LOCUS_ID_COL] == "chr1:500000-2100000"
    assert row[LOCUS_LEAD_EA_COL] == "A"
    assert row[STATUS_COL] == STATUS_CANDIDATE


def test_distant_hits_and_other_chromosomes_are_separate_loci() -> None:
    trait = _trait([(1, 1_000_000, 12.0), (1, 5_000_000, 8.0), (2, 1_000_000, 9.0)])
    loci = find_trait_loci(trait, p_threshold=5e-8, half_window_bp=500_000, exclude=[])
    assert loci.height == 3


def test_variant_inside_a_lead_window_does_not_lead() -> None:
    trait = _trait([(1, 1_000_000, 12.0), (1, 1_400_000, 8.0)])
    loci = find_trait_loci(trait, p_threshold=5e-8, half_window_bp=500_000, exclude=[])
    assert loci[LOCUS_END_COL].to_list() == [1_500_000]


def test_mhc_window_marked_excluded() -> None:
    trait = _trait([(6, 30_000_000, 20.0), (1, 1_000_000, 9.0)])
    loci = find_trait_loci(
        trait,
        p_threshold=5e-8,
        half_window_bp=500_000,
        exclude=[GenomicInterval(chrom=6, start=25_000_000, end=33_000_000)],
    )
    status = dict(zip(loci[LOCUS_LEAD_POS_COL], loci[STATUS_COL], strict=True))
    assert status == {30_000_000: STATUS_EXCLUDED_MHC, 1_000_000: STATUS_CANDIDATE}


def test_no_significant_variant_gives_empty_loci() -> None:
    loci = find_trait_loci(
        _trait([(1, 1_000, 1.0)]), p_threshold=5e-8, half_window_bp=500_000, exclude=[]
    )
    assert loci.height == 0
    assert LOCUS_ID_COL in loci.columns
