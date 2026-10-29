from collections.abc import Sequence

import polars as pl

from mecfs_bio.build_system.task.protein_coloc.candidate_selection import (
    SelectionConfig,
    build_candidates,
    select_candidates,
)
from mecfs_bio.constants.ppp_database_constants import Oid
from mecfs_bio.constants.protein_coloc_constants import (
    CIS_STATUS_CIS,
    CIS_STATUS_COL,
    CIS_STATUS_TRANS,
    CIS_STATUS_UNKNOWN,
    HIT_CHR_COL,
    HIT_GENE_COL,
    HIT_LEAD_MLOG10P_COL,
    HIT_LEAD_POS_COL,
    HIT_N_SIG_VARIANTS_COL,
    HIT_SEG_END_COL,
    HIT_SEG_START_COL,
    LOCUS_CHR_COL,
    LOCUS_END_COL,
    LOCUS_ID_COL,
    LOCUS_LEAD_EA_COL,
    LOCUS_LEAD_MLOG10P_COL,
    LOCUS_LEAD_NEA_COL,
    LOCUS_LEAD_POS_COL,
    LOCUS_START_COL,
    N_OTHER_PROTEINS_COL,
    PROTEIN_ID_COL,
    PROTEIN_MLOG10P_COL,
    SCORE_PROTEIN_COL,
    STATUS_CANDIDATE,
    STATUS_COL,
    STATUS_EXCLUDED_MHC,
    STATUS_NOT_SELECTED_BUDGET,
    STATUS_SELECTED,
)
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval


def _loci(rows: list[tuple[str, int, int, int, float, str]]) -> pl.DataFrame:
    """Rows of (LOCUS_ID, CHR, START, END, LEAD_MLOG10P, STATUS)."""
    return pl.DataFrame(
        {
            LOCUS_ID_COL: [r[0] for r in rows],
            LOCUS_CHR_COL: [r[1] for r in rows],
            LOCUS_START_COL: [r[2] for r in rows],
            LOCUS_END_COL: [r[3] for r in rows],
            LOCUS_LEAD_POS_COL: [(r[2] + r[3]) // 2 for r in rows],
            LOCUS_LEAD_EA_COL: ["A"] * len(rows),
            LOCUS_LEAD_NEA_COL: ["G"] * len(rows),
            LOCUS_LEAD_MLOG10P_COL: [r[4] for r in rows],
            STATUS_COL: [r[5] for r in rows],
        }
    )


def _segments(rows: list[tuple[str, int, int, int, int, float]]) -> pl.DataFrame:
    """Rows of (PROTEIN_ID, CHR, SEG_START, SEG_END, LEAD_POS, LEAD_MLOG10P)."""
    return pl.DataFrame(
        {
            PROTEIN_ID_COL: [r[0] for r in rows],
            HIT_GENE_COL: [f"gene_{r[0]}" for r in rows],
            HIT_CHR_COL: [r[1] for r in rows],
            HIT_SEG_START_COL: [r[2] for r in rows],
            HIT_SEG_END_COL: [r[3] for r in rows],
            HIT_LEAD_POS_COL: [r[4] for r in rows],
            HIT_LEAD_MLOG10P_COL: [r[5] for r in rows],
            HIT_N_SIG_VARIANTS_COL: [1] * len(rows),
        }
    )


def _candidates(
    rows: Sequence[tuple[str, str, float, str, int]],
) -> pl.DataFrame:
    """Rows of (LOCUS_ID, PROTEIN_ID, PROTEIN_MLOG10P, CIS_STATUS, N_OTHER_PROTEINS)."""
    return pl.DataFrame(
        {
            LOCUS_ID_COL: [r[0] for r in rows],
            PROTEIN_ID_COL: [r[1] for r in rows],
            HIT_GENE_COL: [f"gene_{r[1]}" for r in rows],
            PROTEIN_MLOG10P_COL: [r[2] for r in rows],
            CIS_STATUS_COL: [r[3] for r in rows],
            N_OTHER_PROTEINS_COL: [r[4] for r in rows],
        }
    )


_ONE_LOCUS = _loci([("L1", 1, 0, 450_000, 10.0, STATUS_CANDIDATE)])


def _status_by_protein(selected: pl.DataFrame) -> dict[str, str]:
    return dict(zip(selected[PROTEIN_ID_COL], selected[STATUS_COL], strict=True))


def test_segment_overlapping_window_edge_is_a_candidate() -> None:
    # The segment's lead is outside the window, but the segment itself overlaps it.
    candidates = build_candidates(
        _ONE_LOCUS,
        _segments(
            [
                ("OID1", 1, 400_000, 700_000, 650_000, 20.0),
                ("OID2", 1, 460_000, 700_000, 650_000, 20.0),
                ("OID3", 2, 0, 700_000, 100_000, 20.0),
            ]
        ),
        gene_coords={},
        cis_window_bp=1_000_000,
    )
    assert candidates[PROTEIN_ID_COL].to_list() == ["OID1"]
    assert candidates[PROTEIN_MLOG10P_COL].to_list() == [20.0]


def test_cis_status_and_other_protein_count() -> None:
    candidates = build_candidates(
        _ONE_LOCUS,
        _segments(
            [
                ("OID1", 1, 100_000, 200_000, 150_000, 20.0),
                ("OID1", 1, 300_000, 400_000, 350_000, 30.0),
                ("OID2", 1, 100_000, 200_000, 150_000, 20.0),
                ("OID3", 1, 100_000, 200_000, 150_000, 20.0),
            ]
        ),
        gene_coords={
            Oid("OID1"): GenomicInterval(chrom=1, start=1_200_000, end=1_300_000),
            Oid("OID2"): GenomicInterval(chrom=1, start=5_000_000, end=5_100_000),
        },
        cis_window_bp=1_000_000,
    ).sort(PROTEIN_ID_COL)
    assert candidates[CIS_STATUS_COL].to_list() == [
        CIS_STATUS_CIS,
        CIS_STATUS_TRANS,
        CIS_STATUS_UNKNOWN,
    ]
    assert candidates[N_OTHER_PROTEINS_COL].to_list() == [2, 2, 2]
    # A protein with two overlapping segments scores its stronger lead.
    assert candidates[PROTEIN_MLOG10P_COL].to_list() == [30.0, 20.0, 20.0]


def test_excluded_loci_get_no_candidates() -> None:
    loci = _loci([("MHC", 6, 0, 1_000_000, 50.0, STATUS_EXCLUDED_MHC)])
    candidates = build_candidates(
        loci,
        _segments([("OID1", 6, 100, 200, 150, 20.0)]),
        gene_coords={},
        cis_window_bp=1_000_000,
    )
    assert candidates.height == 0


def test_cis_pair_not_penalized_at_hotspot() -> None:
    selected = select_candidates(
        _candidates(
            [
                ("L1", "OID_TRANS", 20.0, CIS_STATUS_TRANS, 50),
                ("L1", "OID_CIS", 20.0, CIS_STATUS_CIS, 50),
            ]
        ),
        _ONE_LOCUS,
        SelectionConfig(max_pairs=1),
    )
    assert _status_by_protein(selected) == {
        "OID_CIS": STATUS_SELECTED,
        "OID_TRANS": STATUS_NOT_SELECTED_BUDGET,
    }


def test_unknown_gene_coords_treated_like_trans() -> None:
    selected = select_candidates(
        _candidates(
            [
                ("L1", "OID_TRANS", 20.0, CIS_STATUS_TRANS, 50),
                ("L1", "OID_UNKNOWN", 20.0, CIS_STATUS_UNKNOWN, 50),
                ("L1", "OID_CIS", 20.0, CIS_STATUS_CIS, 50),
            ]
        ),
        _ONE_LOCUS,
        SelectionConfig(),
    )
    score = dict(zip(selected[PROTEIN_ID_COL], selected["SCORE"], strict=True))
    assert score["OID_UNKNOWN"] == score["OID_TRANS"] < score["OID_CIS"]


def test_per_locus_cap_respected() -> None:
    candidates = _candidates(
        [("L1", f"OID{i}", 20.0 + i, CIS_STATUS_CIS, 7) for i in range(8)]
    )
    selected = select_candidates(
        candidates, _ONE_LOCUS, SelectionConfig(max_per_locus=5)
    )
    assert selected.filter(pl.col(STATUS_COL) == STATUS_SELECTED).height == 5
    # The five strongest proteins win.
    assert set(
        selected.filter(pl.col(STATUS_COL) == STATUS_SELECTED)[PROTEIN_ID_COL]
    ) == {f"OID{i}" for i in range(3, 8)}


def test_pair_budget_respected_and_rest_marked() -> None:
    loci = _loci(
        [
            (f"L{i}", 1, i * 10_000_000, i * 10_000_000 + 1, 10.0, STATUS_CANDIDATE)
            for i in range(10)
        ]
    )
    candidates = _candidates(
        [
            (f"L{i}", f"OID{i}_{j}", 20.0, CIS_STATUS_CIS, 3)
            for i in range(10)
            for j in range(4)
        ]
    )
    selected = select_candidates(candidates, loci, SelectionConfig(max_pairs=30))
    counts = dict(selected.group_by(STATUS_COL).len().rows())
    assert counts == {STATUS_SELECTED: 30, STATUS_NOT_SELECTED_BUDGET: 10}


def test_clipping_limits_protein_term() -> None:
    selected = select_candidates(
        _candidates(
            [
                ("L1", "OID_A", 300.0, CIS_STATUS_CIS, 1),
                ("L1", "OID_B", 60.0, CIS_STATUS_CIS, 1),
            ]
        ),
        _ONE_LOCUS,
        SelectionConfig(mlog10p_clip=50.0),
    )
    assert selected[SCORE_PROTEIN_COL].to_list() == [50.0, 50.0]


def test_selection_is_deterministic_under_ties() -> None:
    rows = [("L1", f"OID{i}", 20.0, CIS_STATUS_CIS, 9) for i in (3, 1, 2)]
    first = select_candidates(
        _candidates(rows), _ONE_LOCUS, SelectionConfig(max_pairs=2)
    )
    second = select_candidates(
        _candidates(list(reversed(rows))), _ONE_LOCUS, SelectionConfig(max_pairs=2)
    )
    assert _status_by_protein(first) == _status_by_protein(second)
    assert _status_by_protein(first) == {
        "OID1": STATUS_SELECTED,
        "OID2": STATUS_SELECTED,
        "OID3": STATUS_NOT_SELECTED_BUDGET,
    }
