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
