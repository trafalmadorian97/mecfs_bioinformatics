"""Task-level tests of GenomeReferenceHarmonizationTask on a synthetic genome and panel."""

from pathlib import Path

import attrs
import narwhals
import polars as pl
import pytest
from attrs import frozen

from mecfs_bio.build_system.task.genome_reference_harmonization.flip import (
    FLIP_NEGATE,
    ExtraColumnRule,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_HAZARD_RATIO_95L_COL,
    GWASLAB_HAZARD_RATIO_95U_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_ODDS_RATIO_95L_COL,
    GWASLAB_ODDS_RATIO_95U_COL,
    GWASLAB_ODDS_RATIO_COL,
    GWASLAB_POS_COL,
)
from mecfs_bio.constants.regenie_constants import REGENIE_A1FREQ_CASES_COL
from test_mecfs_bio.unit.build_system.task.genome_reference_harmonization.genome_reference_fixtures import (
    CONSISTENT_INDEL,
    CONSISTENT_SNV,
    INCONSISTENT_SNV,
    TEST_OPTIONS,
    PanelRecord,
    Variant,
    positions,
    row_at,
    run_harmonization,
    sumstats_frame,
)

EA = GWASLAB_EFFECT_ALLELE_COL
NEA = GWASLAB_NON_EFFECT_ALLELE_COL
EAF = GWASLAB_EFFECT_ALLELE_FREQ_COL
BETA = GWASLAB_BETA_COL
OR_95L = GWASLAB_ODDS_RATIO_95L_COL
OR_95U = GWASLAB_ODDS_RATIO_95U_COL
HR_95L = GWASLAB_HAZARD_RATIO_95L_COL
HR_95U = GWASLAB_HAZARD_RATIO_95U_COL
_UNREGISTERED_COLUMN = "MYSTERY_STATISTIC"
_EXTRA_COLUMN = "BETA_ALTERNATIVE_MODEL"


def test_snvs_are_oriented_against_the_reference(tmp_path: Path) -> None:
    variants = [
        CONSISTENT_SNV,  # 1: NEA is the reference base -> kept
        INCONSISTENT_SNV,  # 2: EA is the reference base -> swapped
        Variant(
            pos=3, ea="A", nea="C", beta=0.3
        ),  # complement of NEA is ref -> complemented
        Variant(
            pos=11, ea="G", nea="A", beta=0.4
        ),  # complement of EA is ref -> complemented and swapped
        Variant(pos=12, ea="C", nea="G"),  # nothing matches -> dropped
    ]
    result = run_harmonization(tmp_path / "run", sumstats_frame(variants))
    assert result.select(GWASLAB_POS_COL, EA, NEA, BETA, EAF).rows() == [
        (1, "G", "A", pytest.approx(0.1), pytest.approx(0.3)),
        (2, "T", "C", pytest.approx(-0.2), pytest.approx(0.7)),
        (3, "T", "G", pytest.approx(0.3), pytest.approx(0.3)),
        (11, "T", "C", pytest.approx(-0.4), pytest.approx(0.7)),
    ]


def test_indels_with_one_matching_allele_are_oriented_by_it(tmp_path: Path) -> None:
    variants = [
        CONSISTENT_SNV,
        CONSISTENT_INDEL,  # 21: only NEA G matches -> kept
        Variant(pos=22, ea="A", nea="AG", beta=0.5),  # only EA matches -> swapped
        Variant(pos=23, ea="C", nea="CA"),  # neither matches -> dropped
    ]
    result = run_harmonization(tmp_path / "run", sumstats_frame(variants))
    assert positions(result) == [1, 21, 22]
    swapped = row_at(result, 22)
    assert (swapped[EA], swapped[NEA]) == ("AG", "A")
    assert swapped[BETA] == pytest.approx(-0.5)


def test_trusted_table_keeps_ambiguous_indel_in_source_orientation(
    tmp_path: Path,
) -> None:
    ambiguous = Variant(pos=5, ea="T", nea="TG")  # both alleles match the T-G-run
    result = run_harmonization(
        tmp_path / "run", sumstats_frame([CONSISTENT_SNV, CONSISTENT_INDEL, ambiguous])
    )
    kept = row_at(result, 5)
    assert (kept[EA], kept[NEA]) == ("T", "TG")


@pytest.mark.parametrize("trusted", [True, False])
def test_palindromic_mnp_is_kept_only_when_trusted(
    tmp_path: Path, trusted: bool
) -> None:
    palindromic_mnp = Variant(pos=41, ea="GT", nea="AC")
    variants = [CONSISTENT_SNV, CONSISTENT_INDEL, palindromic_mnp]
    if not trusted:
        variants.append(INCONSISTENT_SNV)
    result = run_harmonization(tmp_path / "run", sumstats_frame(variants))
    assert (41 in positions(result)) == trusted


def test_invalid_alleles_are_dropped_and_lowercase_is_accepted(tmp_path: Path) -> None:
    variants = [
        CONSISTENT_SNV,
        Variant(pos=3, ea="N", nea="C"),
        Variant(pos=4, ea="T", nea="T"),
        Variant(pos=51, ea="a", nea="g"),
    ]
    result = run_harmonization(tmp_path / "run", sumstats_frame(variants))
    assert positions(result) == [1, 51]
    assert (row_at(result, 51)[EA], row_at(result, 51)[NEA]) == ("A", "G")


def test_flip_covers_allele_frequency_columns_and_confidence_bounds(
    tmp_path: Path,
) -> None:
    frame = sumstats_frame([CONSISTENT_SNV, INCONSISTENT_SNV]).with_columns(
        pl.Series(REGENIE_A1FREQ_CASES_COL, [0.3, 0.25]),
        pl.Series(GWASLAB_ODDS_RATIO_COL, [2.0, 4.0]),
        pl.Series(OR_95L, [1.5, 2.0]),
        pl.Series(OR_95U, [2.5, 8.0]),
        pl.Series(HR_95L, [1.2, 1.25]),
        pl.Series(HR_95U, [1.8, 5.0]),
    )
    result = run_harmonization(tmp_path / "run", frame)
    unchanged = row_at(result, 1)
    assert (unchanged[OR_95L], unchanged[OR_95U]) == (
        pytest.approx(1.5),
        pytest.approx(2.5),
    )
    swapped = row_at(result, 2)
    assert swapped[REGENIE_A1FREQ_CASES_COL] == pytest.approx(0.75)
    assert swapped[GWASLAB_ODDS_RATIO_COL] == pytest.approx(0.25)
    assert (swapped[OR_95L], swapped[OR_95U]) == (
        pytest.approx(0.125),
        pytest.approx(0.5),
    )
    assert (swapped[HR_95L], swapped[HR_95U]) == (
        pytest.approx(0.2),
        pytest.approx(0.8),
    )


def test_unregistered_column_fails(tmp_path: Path) -> None:
    frame = sumstats_frame([CONSISTENT_SNV]).with_columns(
        pl.lit(1.0).alias(_UNREGISTERED_COLUMN)
    )
    with pytest.raises(AssertionError):
        run_harmonization(tmp_path / "run", frame)


def test_extra_column_rule_is_applied(tmp_path: Path) -> None:
    frame = sumstats_frame([CONSISTENT_SNV, INCONSISTENT_SNV]).with_columns(
        pl.Series(_EXTRA_COLUMN, [0.5, 0.5])
    )
    options = attrs.evolve(
        TEST_OPTIONS,
        extra_column_rules=(ExtraColumnRule(column=_EXTRA_COLUMN, rule=FLIP_NEGATE),),
    )
    result = run_harmonization(tmp_path / "run", frame, options=options)
    assert row_at(result, 1)[_EXTRA_COLUMN] == pytest.approx(0.5)
    assert row_at(result, 2)[_EXTRA_COLUMN] == pytest.approx(-0.5)


def test_output_does_not_depend_on_input_row_order(tmp_path: Path) -> None:
    chr1_keep, chr1_swap = CONSISTENT_SNV, INCONSISTENT_SNV
    chr2_keep = Variant(chrom=2, pos=5, ea="G", nea="A")
    chr2_swap = Variant(chrom=2, pos=9, ea="C", nea="T")
    ordered = run_harmonization(
        tmp_path / "ordered",
        sumstats_frame([chr1_keep, chr1_swap, chr2_keep, chr2_swap]),
    )
    interleaved = run_harmonization(
        tmp_path / "interleaved",
        sumstats_frame([chr2_swap, chr1_keep, chr2_keep, chr1_swap]),
    )
    assert interleaved.equals(ordered)
    assert ordered[GWASLAB_CHROM_COL].to_list() == [1, 1, 2, 2]


def test_excluded_chromosome_rows_are_dropped(tmp_path: Path) -> None:
    mitochondrial = Variant(chrom=25, pos=1, ea="G", nea="A")
    result = run_harmonization(
        tmp_path / "run", sumstats_frame([CONSISTENT_SNV, mitochondrial])
    )
    assert result[GWASLAB_CHROM_COL].to_list() == [1]


def test_chromosome_missing_from_fasta_fails(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        run_harmonization(
            tmp_path / "run",
            sumstats_frame([CONSISTENT_SNV, Variant(chrom=3, pos=1, ea="G", nea="A")]),
        )


def test_null_position_fails(tmp_path: Path) -> None:
    frame = sumstats_frame([CONSISTENT_SNV]).with_columns(
        pl.lit(None, dtype=pl.Int64).alias(GWASLAB_POS_COL)
    )
    with pytest.raises(AssertionError):
        run_harmonization(tmp_path / "run", frame)


def test_duplicate_variant_after_orientation_fails(tmp_path: Path) -> None:
    variants = [CONSISTENT_SNV, INCONSISTENT_SNV, Variant(pos=2, ea="T", nea="C")]
    with pytest.raises(AssertionError):
        run_harmonization(tmp_path / "run", sumstats_frame(variants))


def test_long_alleles_are_classified_with_a_tiny_gather_budget(tmp_path: Path) -> None:
    long_mnp = Variant(
        pos=31, ea="T" * 10, nea="T" * 9 + "A", beta=0.2
    )  # EA matches -> swap
    options = attrs.evolve(TEST_OPTIONS, max_gather_bytes=1)
    result = run_harmonization(
        tmp_path / "run", sumstats_frame([CONSISTENT_SNV, long_mnp]), options=options
    )
    swapped = row_at(result, 31)
    assert (swapped[EA], swapped[NEA]) == ("T" * 9 + "A", "T" * 10)
    assert swapped[BETA] == pytest.approx(-0.2)


@frozen
class _SwitchToDuckDbPipe(DataProcessingPipe):
    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.collect().lazy(backend="duckdb")


def test_pipe_that_changes_backend_fails(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        run_harmonization(
            tmp_path / "run",
            sumstats_frame([CONSISTENT_SNV]),
            pipe=_SwitchToDuckDbPipe(),
        )


_PALINDROME_PANEL = [
    PanelRecord(pos=13, ref="T", alt="A", af=0.15),
    PanelRecord(pos=14, ref="G", alt="C", af=0.85),
    PanelRecord(pos=18, ref="T", alt="A", af=0.1),
    PanelRecord(pos=19, ref="T", alt="A", af=0.5),
    PanelRecord(pos=27, ref="T", alt="A", af=0.45),
]
_PALINDROMES = [
    Variant(
        pos=13, ea="A", nea="T", eaf=0.1, beta=0.3
    ),  # same side of 0.5 as the panel -> kept
    Variant(
        pos=14, ea="G", nea="C", eaf=0.8, beta=0.3
    ),  # EA is ref: swapped, then strand-flipped
    Variant(
        pos=18, ea="A", nea="T", eaf=0.9, beta=0.3
    ),  # opposite side of 0.5 -> strand-flipped
    Variant(pos=19, ea="A", nea="T", eaf=0.5),  # sumstats MAF above 0.4 -> unresolved
    Variant(pos=24, ea="G", nea="C", eaf=0.1),  # no panel record -> unresolved
    Variant(pos=27, ea="A", nea="T", eaf=0.1),  # panel MAF above 0.4 -> unresolved
    Variant(pos=30, ea="T", nea="A", eaf=None),  # no EAF -> unresolved
]
_UNRESOLVED_PALINDROME_POSITIONS = [19, 24, 27, 30]


def test_untrusted_palindromes_are_resolved_by_panel_frequency(tmp_path: Path) -> None:
    variants = [CONSISTENT_SNV, INCONSISTENT_SNV, *_PALINDROMES]
    result = run_harmonization(
        tmp_path / "run", sumstats_frame(variants), panel=_PALINDROME_PANEL
    )
    assert positions(result) == [1, 2, 13, 14, 18]
    kept = row_at(result, 13)
    assert (kept[EA], kept[NEA], kept[BETA]) == ("A", "T", pytest.approx(0.3))
    swapped_then_flipped = row_at(result, 14)
    assert (swapped_then_flipped[EA], swapped_then_flipped[NEA]) == ("C", "G")
    assert swapped_then_flipped[BETA] == pytest.approx(0.3)
    assert swapped_then_flipped[EAF] == pytest.approx(0.8)
    flipped = row_at(result, 18)
    assert (flipped[EA], flipped[NEA]) == ("A", "T")
    assert (flipped[BETA], flipped[EAF]) == (pytest.approx(-0.3), pytest.approx(0.1))


def test_unresolved_palindromes_can_be_kept(tmp_path: Path) -> None:
    variants = [CONSISTENT_SNV, INCONSISTENT_SNV, *_PALINDROMES]
    options = attrs.evolve(TEST_OPTIONS, keep_unresolved_palindromes=True)
    result = run_harmonization(
        tmp_path / "run",
        sumstats_frame(variants),
        panel=_PALINDROME_PANEL,
        options=options,
    )
    assert set(_UNRESOLVED_PALINDROME_POSITIONS) <= set(positions(result))


def test_trusted_palindromes_keep_source_strand(tmp_path: Path) -> None:
    opposite_side = Variant(pos=18, ea="A", nea="T", eaf=0.9, beta=0.3)
    result = run_harmonization(
        tmp_path / "run",
        sumstats_frame([CONSISTENT_SNV, CONSISTENT_INDEL, opposite_side]),
        panel=_PALINDROME_PANEL,
    )
    assert row_at(result, 18)[BETA] == pytest.approx(0.3)


_AMBIGUOUS_DISTANCE = 0.1
_AMBIGUOUS_MARGIN = 0.2
_STRINGENT_OPTIONS = attrs.evolve(
    TEST_OPTIONS,
    indel_max_af_distance=_AMBIGUOUS_DISTANCE,
    indel_min_af_margin=_AMBIGUOUS_MARGIN,
)
# Every variant below is T/TT inside the chr1 T homopolymer at 31-40, so both alleles
# match the genome. "keep" reads the row as REF=TT, ALT=T; "flip" as REF=T, ALT=TT.
_AMBIGUOUS_INDELS = [
    Variant(pos=31, ea="T", nea="TT", eaf=None),  # no EAF -> dropped
    Variant(pos=32, ea="T", nea="TT", eaf=0.3),  # no panel record -> dropped
    Variant(pos=33, ea="T", nea="TT", eaf=0.3),  # keep record fits -> kept
    Variant(pos=34, ea="T", nea="TT", eaf=0.3, beta=0.2),  # flip record fits -> swapped
    Variant(pos=35, ea="T", nea="TT", eaf=0.3),  # only record misfits -> dropped
    Variant(pos=36, ea="T", nea="TT", eaf=0.3),  # both records, keep decisive -> kept
    Variant(
        pos=37, ea="T", nea="TT", eaf=0.5
    ),  # both records, margin too small -> dropped
]
_AMBIGUOUS_PANEL = [
    PanelRecord(pos=33, ref="TT", alt="T", af=0.32),
    PanelRecord(pos=34, ref="T", alt="TT", af=0.68),
    PanelRecord(pos=35, ref="TT", alt="T", af=0.8),
    PanelRecord(pos=36, ref="TT", alt="T", af=0.31),
    # Polymorphic flip record (AF 0.4 -> predicted EAF 0.6): the keep reading still wins
    # decisively. A monomorphic (AF 0) record here would be dropped by the AF filter, so a
    # polymorphic misfit keeps the "both records, keep decisive" case genuine.
    PanelRecord(pos=36, ref="T", alt="TT", af=0.4),
    PanelRecord(pos=37, ref="TT", alt="T", af=0.5),
    PanelRecord(pos=37, ref="T", alt="TT", af=0.45),
]


def test_untrusted_ambiguous_indels_follow_the_stringent_rules(tmp_path: Path) -> None:
    variants = [CONSISTENT_SNV, INCONSISTENT_SNV, *_AMBIGUOUS_INDELS]
    result = run_harmonization(
        tmp_path / "run",
        sumstats_frame(variants),
        panel=_AMBIGUOUS_PANEL,
        options=_STRINGENT_OPTIONS,
    )
    assert positions(result) == [1, 2, 33, 34, 36]
    swapped = row_at(result, 34)
    assert (swapped[EA], swapped[NEA]) == ("TT", "T")
    assert (swapped[EAF], swapped[BETA]) == (pytest.approx(0.7), pytest.approx(-0.2))


def test_ambiguous_indel_resolution_is_symmetric_in_orientation(tmp_path: Path) -> None:
    as_deletion_label = Variant(pos=36, ea="T", nea="TT", eaf=0.3, beta=0.2)
    as_insertion_label = Variant(pos=36, ea="TT", nea="T", eaf=0.7, beta=-0.2)
    results = [
        row_at(
            run_harmonization(
                tmp_path / name,
                sumstats_frame([CONSISTENT_SNV, INCONSISTENT_SNV, variant]),
                panel=_AMBIGUOUS_PANEL,
                options=_STRINGENT_OPTIONS,
            ),
            36,
        )
        for name, variant in [
            ("deletion", as_deletion_label),
            ("insertion", as_insertion_label),
        ]
    ]
    for row in results:
        assert (row[EA], row[NEA]) == ("T", "TT")
        assert (row[EAF], row[BETA]) == (pytest.approx(0.3), pytest.approx(0.2))


def test_monomorphic_panel_records_are_ignored(tmp_path: Path) -> None:
    # An AF-0 panel record describes a variant absent from EUR, so a match to it is spurious
    # and it is treated as absent. Without that rule the keep reading (EAF 0.05 vs AF 0.0)
    # fits within tolerance and the indel is kept; with it the indel is dropped (not in panel).
    monomorphic = Variant(pos=38, ea="T", nea="TT", eaf=0.05)
    result = run_harmonization(
        tmp_path / "run",
        sumstats_frame([CONSISTENT_SNV, INCONSISTENT_SNV, monomorphic]),
        panel=[PanelRecord(pos=38, ref="TT", alt="T", af=0.0)],
        options=_STRINGENT_OPTIONS,
    )
    assert 38 not in positions(result)


_AMBIGUOUS_WITHOUT_PANEL = Variant(pos=5, ea="T", nea="TG")


@pytest.mark.parametrize(
    ("extra_variants", "options", "expect_trusted"),
    [
        ([], TEST_OPTIONS, True),
        ([INCONSISTENT_SNV], TEST_OPTIONS, False),
        ([], attrs.evolve(TEST_OPTIONS, min_checkable_snvs=2), False),
        ([], attrs.evolve(TEST_OPTIONS, min_checkable_indels=2), False),
    ],
)
def test_trust_decision_controls_ambiguous_indels(
    tmp_path: Path,
    extra_variants: list[Variant],
    options: GenomeReferenceHarmonizationOptions,
    expect_trusted: bool,
) -> None:
    variants = [
        CONSISTENT_SNV,
        CONSISTENT_INDEL,
        _AMBIGUOUS_WITHOUT_PANEL,
        *extra_variants,
    ]
    result = run_harmonization(
        tmp_path / "run", sumstats_frame(variants), options=options
    )
    # Trusted tables keep the ambiguous indel; untrusted ones drop it (no panel record).
    assert (5 in positions(result)) == expect_trusted


def test_suspicious_ambiguous_indels_make_a_table_untrusted(tmp_path: Path) -> None:
    # pos 5 is T at the head of a G run, so T/TG is class BOTH. The panel has only the
    # opposite reading (REF=T, ALT=TG) at AF 0.9, and EAF 0.1 = 1 - 0.9, so the row is
    # suspicious. With min_checkable_ambiguous_indels=1 that makes the table untrusted,
    # and the ambiguous indel is swapped rather than kept in source orientation.
    ambiguous = Variant(pos=5, ea="T", nea="TG", eaf=0.1)
    panel = [PanelRecord(pos=5, ref="T", alt="TG", af=0.9)]
    result = run_harmonization(
        tmp_path / "run",
        sumstats_frame([CONSISTENT_SNV, CONSISTENT_INDEL, ambiguous]),
        panel=panel,
    )
    swapped = row_at(result, 5)
    assert (swapped[EA], swapped[NEA]) == ("TG", "T")


def test_table_without_eaf_is_untrusted(tmp_path: Path) -> None:
    # With no EAF column, ambiguous-indel orientation cannot be assessed, so the table is
    # refused trust and its ambiguous indel is dropped (NO_EAF) rather than kept.
    ambiguous = Variant(pos=5, ea="T", nea="TG")
    frame = sumstats_frame([CONSISTENT_SNV, CONSISTENT_INDEL, ambiguous]).drop(
        GWASLAB_EFFECT_ALLELE_FREQ_COL
    )
    result = run_harmonization(tmp_path / "run", frame)
    assert 5 not in positions(result)
    assert set(positions(result)) >= {1, 21}
