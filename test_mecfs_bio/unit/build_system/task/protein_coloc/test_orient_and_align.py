from pathlib import Path

import polars as pl
import pytest

from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    DROP_NOT_ON_REFERENCE,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.task.protein_coloc.orient_and_align import (
    OrientedWindow,
    align_windows,
    orient_window,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    BETA_PROTEIN_COL,
    BETA_TRAIT_COL,
    WINDOW_SCHEMA,
)
from test_mecfs_bio.unit.build_system.task.genome_reference_harmonization.genome_reference_fixtures import (
    TEST_OPTIONS,
    write_fasta,
)

# Synthetic chr1 (genome_reference_fixtures): position 2 is C, position 1 is A, and
# positions 52-54 read TCA.
_MIRROR_POS = 52


def _window(rows: list[tuple[int, str, str, float, float]]) -> pl.DataFrame:
    """Rows of (POS, EA, NEA, BETA, EAF) on chr1 in WINDOW_SCHEMA."""
    return pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [1] * len(rows),
            GWASLAB_POS_COL: [r[0] for r in rows],
            GWASLAB_EFFECT_ALLELE_COL: [r[1] for r in rows],
            GWASLAB_NON_EFFECT_ALLELE_COL: [r[2] for r in rows],
            GWASLAB_BETA_COL: [r[3] for r in rows],
            "SE": [0.01] * len(rows),
            GWASLAB_EFFECT_ALLELE_FREQ_COL: [r[4] for r in rows],
            "N": [1000.0] * len(rows),
        },
        schema=WINDOW_SCHEMA,
    )


def _panel(path: Path, records: list[tuple[int, str, str, float]]) -> Path:
    pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [1] * len(records),
            GWASLAB_POS_COL: [r[0] for r in records],
            PANEL_REF_COL: [r[1] for r in records],
            PANEL_ALT_COL: [r[2] for r in records],
            PANEL_AF_COL: [r[3] for r in records],
        },
        schema={
            GWASLAB_CHROM_COL: pl.Int32,
            GWASLAB_POS_COL: pl.Int32,
            PANEL_REF_COL: pl.String,
            PANEL_ALT_COL: pl.String,
            PANEL_AF_COL: pl.Float32,
        },
    ).write_parquet(path)
    return path


def _orient(
    tmp_path: Path,
    rows: pl.DataFrame,
    panel_records: list[tuple[int, str, str, float]] | None = None,
) -> OrientedWindow:
    fasta_dir = tmp_path / "fasta"
    write_fasta(fasta_dir)
    panel_path = _panel(tmp_path / "panel.parquet", panel_records or [])
    return orient_window(
        rows,
        chrom=1,
        fasta=IndexedFasta.open(fasta_dir),
        panel_path=panel_path,
        panel_af_col=PANEL_AF_COL,
        options=TEST_OPTIONS,
    )


def test_swapped_snv_is_flipped(tmp_path: Path) -> None:
    oriented = _orient(tmp_path, _window([(2, "C", "T", 0.2, 0.3)]))
    row = oriented.frame.row(0, named=True)
    assert row[GWASLAB_NON_EFFECT_ALLELE_COL] == "C"
    assert row[GWASLAB_EFFECT_ALLELE_COL] == "T"
    assert row[GWASLAB_BETA_COL] == pytest.approx(-0.2)
    assert row[GWASLAB_EFFECT_ALLELE_FREQ_COL] == pytest.approx(0.7)
    assert oriented.frame.schema == WINDOW_SCHEMA


def test_unorientable_snv_is_dropped_and_counted(tmp_path: Path) -> None:
    oriented = _orient(
        tmp_path, _window([(1, "G", "C", 0.2, 0.3), (2, "T", "C", 0.1, 0.3)])
    )
    assert oriented.frame[GWASLAB_POS_COL].to_list() == [2]
    assert oriented.n_input == 2
    assert oriented.drop_counts == {DROP_NOT_ON_REFERENCE: 1}


def test_flipped_deletion_does_not_join_trait_insertion(tmp_path: Path) -> None:
    # Both T>TCA (insertion) and TCA>T (deletion) are on-reference at _MIRROR_POS.
    trait = _window([(_MIRROR_POS, "TCA", "T", 0.3, 0.05)])
    # The protein carries the deletion (panel REF=TCA, ALT=T, AF 0.1), stored flipped.
    protein = _window([(_MIRROR_POS, "TCA", "T", 0.4, 0.9)])
    oriented = _orient(
        tmp_path, protein, panel_records=[(_MIRROR_POS, "TCA", "T", 0.1)]
    )
    assert oriented.frame.select(
        GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
    ).rows() == [("T", "TCA")]
    aligned = align_windows(trait=trait, protein=oriented.frame)
    assert aligned.frame.filter(pl.col(GWASLAB_POS_COL) == _MIRROR_POS).height == 0


def test_duplicate_keys_are_dropped_on_both_sides() -> None:
    trait = _window(
        [(2, "T", "C", 0.1, 0.3), (2, "T", "C", 0.2, 0.3), (5, "G", "T", 0.1, 0.3)]
    )
    protein = _window(
        [(5, "G", "T", 0.5, 0.3), (5, "G", "T", 0.6, 0.3), (2, "T", "C", 0.1, 0.3)]
    )
    aligned = align_windows(trait=trait, protein=protein)
    assert aligned.n_duplicate_keys_dropped == 2
    assert aligned.frame.height == 0


def test_alignment_is_exact_on_ordered_alleles() -> None:
    trait = _window([(2, "T", "C", 0.1, 0.3), (5, "G", "T", 0.2, 0.3)])
    # Same position and allele pair at 5 but reversed order: a different oriented variant.
    protein = _window([(2, "T", "C", 0.5, 0.3), (5, "T", "G", 0.6, 0.3)])
    aligned = align_windows(trait=trait, protein=protein)
    assert aligned.frame[GWASLAB_POS_COL].to_list() == [2]
    assert aligned.frame[BETA_TRAIT_COL].to_list() == [0.1]
    assert aligned.frame[BETA_PROTEIN_COL].to_list() == [0.5]
    assert aligned.n_trait == 2 and aligned.n_protein == 2


def test_null_eaf_rows_survive_orientation(tmp_path: Path) -> None:
    # Position 1 is A: A/T is a palindrome, unresolvable without a frequency. Position 2 is C:
    # C/T (stored flipped) is oriented by the reference allele alone.
    rows = _window([(1, "T", "A", 0.2, 0.0), (2, "C", "T", 0.1, 0.0)]).with_columns(
        pl.lit(None, dtype=pl.Float64).alias(GWASLAB_EFFECT_ALLELE_FREQ_COL)
    )
    oriented = _orient(tmp_path, rows)
    assert oriented.frame.select(
        GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
    ).rows() == [(2, "T", "C")]
    assert oriented.frame[GWASLAB_BETA_COL].item() == pytest.approx(-0.1)
    assert oriented.frame[GWASLAB_EFFECT_ALLELE_FREQ_COL].is_null().all()
    assert sum(oriented.drop_counts.values()) == 1


def test_eaf_discordant_variants_are_dropped_and_counted() -> None:
    # Position 2 differs by 0.3 (dropped), 5 by 0.15 (kept); 7 has no protein EAF (kept).
    trait = _window(
        [(2, "T", "C", 0.1, 0.3), (5, "G", "T", 0.2, 0.3), (7, "A", "G", 0.3, 0.3)]
    )
    protein = _window(
        [(2, "T", "C", 0.5, 0.6), (5, "G", "T", 0.6, 0.45), (7, "A", "G", 0.7, 0.3)]
    ).with_columns(
        pl.when(pl.col(GWASLAB_POS_COL) == 7)
        .then(None)
        .otherwise(pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL))
        .alias(GWASLAB_EFFECT_ALLELE_FREQ_COL)
    )
    aligned = align_windows(trait=trait, protein=protein, max_eaf_difference=0.2)
    assert aligned.frame[GWASLAB_POS_COL].to_list() == [5, 7]
    assert aligned.n_eaf_discordant_dropped == 1

    unchecked = align_windows(trait=trait, protein=protein, max_eaf_difference=None)
    assert unchecked.frame.height == 3
    assert unchecked.n_eaf_discordant_dropped == 0


def test_multi_record_trait_positions_are_dropped_only_when_asked() -> None:
    # Position 2 has two trait records (a split multi-allelic site); 5 has one.
    trait = _window(
        [(2, "T", "C", 0.1, 0.3), (2, "G", "C", 0.2, 0.1), (5, "G", "T", 0.3, 0.3)]
    )
    protein = _window([(2, "T", "C", 0.5, 0.3), (5, "G", "T", 0.6, 0.3)])
    strict = align_windows(
        trait=trait, protein=protein, drop_multi_record_trait_positions=True
    )
    assert strict.frame[GWASLAB_POS_COL].to_list() == [5]
    assert strict.n_trait_multi_record_dropped == 2

    default = align_windows(trait=trait, protein=protein)
    assert default.frame[GWASLAB_POS_COL].to_list() == [2, 5]
    assert default.n_trait_multi_record_dropped == 0
