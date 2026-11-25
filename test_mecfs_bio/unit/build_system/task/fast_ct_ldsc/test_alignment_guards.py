"""Alignment-guard thresholds against synthetic diagnostics."""

import pytest

from mecfs_bio.build_system.task.fast_ct_ldsc.alignment_guards import (
    AlignmentDiagnostics,
    AlignmentGuards,
    check_alignment_guards,
)
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild


def _diag(
    *,
    n_found: int = 950,
    n_matched: int = 945,
    n_position_compared: int | None = 945,
    n_position_concordant: int | None = 940,
    declared_build: GenomeBuild | None = "19",
) -> AlignmentDiagnostics:
    return AlignmentDiagnostics(
        n_context=1000,
        n_found=n_found,
        n_matched=n_matched,
        n_position_compared=n_position_compared,
        n_position_concordant=n_position_concordant,
        declared_build=declared_build,
    )


def test_guards_pass_on_clean_diagnostics() -> None:
    check_alignment_guards(_diag(), AlignmentGuards(), "T")


@pytest.mark.parametrize(
    "diag",
    [
        _diag(n_found=300, n_matched=300, n_position_compared=300),
        _diag(n_found=950, n_matched=700, n_position_compared=700),
        _diag(n_position_concordant=100),
    ],
    ids=["low_coverage", "high_mismatch", "low_concordance_build19"],
)
def test_guards_trip(diag: AlignmentDiagnostics) -> None:
    with pytest.raises(AssertionError):
        check_alignment_guards(diag, AlignmentGuards(), "T")


@pytest.mark.parametrize("declared_build", [None, "38"])
def test_position_guard_skipped_without_declared_build_19(
    declared_build: GenomeBuild | None,
) -> None:
    check_alignment_guards(
        _diag(n_position_concordant=100, declared_build=declared_build),
        AlignmentGuards(),
        "T",
    )


def test_permissive_override_lets_tripping_diagnostics_through() -> None:
    check_alignment_guards(
        _diag(n_found=300, n_matched=200, n_position_compared=200),
        AlignmentGuards(
            min_coverage=0.0,
            max_allele_mismatch_fraction=1.0,
            min_position_concordance=0.0,
        ),
        "T",
    )


def test_diagnostics_json_roundtrip() -> None:
    for d in (_diag(), _diag(n_position_compared=None, n_position_concordant=None)):
        assert AlignmentDiagnostics.from_json_dict(d.to_json_dict()) == d


def test_guard_thresholds_validated() -> None:
    with pytest.raises(AssertionError):
        AlignmentGuards(min_coverage=1.5)
