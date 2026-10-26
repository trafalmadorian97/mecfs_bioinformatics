"""
Fail-fast guards against upstream harmonization bugs in a trait aligned for fast CT-LDSC.

The trait is joined to the HapMap3 regression SNPs by rsID, so a trait whose rsIDs were assigned
wrongly upstream (for example from positions in the wrong genome build) would silently yield
z-scores for the wrong variants. Three diagnostics expose that, all over the context SNP set:

- coverage = n_matched / n_context: the fraction of context SNPs the trait supplies with alleles
  that match the index (directly or swapped).
- allele_mismatch_fraction = (n_found - n_matched) / n_found: among context rsIDs the trait
  carries, the fraction whose alleles do not match the index. Correct rsIDs mismatch well under
  1%; rsIDs pointing at the wrong variant mismatch at a large rate.
- position_concordance = n_position_concordant / n_position_compared: among matched SNPs with a
  trait position, the fraction whose (CHR, POS) equals the index's hg19 position. Only meaningful
  when the trait declares build 19, so the guard is applied only then.
"""

from __future__ import annotations

from attrs import frozen

from mecfs_bio.build_system.task.ppp_ldsc.trait_alignment import TraitAligned
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild


@frozen(slots=True)
class AlignmentGuards:
    """Thresholds for the alignment diagnostics; each lies in [0, 1].

    The defaults leave wide margin around the values observed on the 13 traits of the initial
    genetic correlation study (experiments/claude/fast_ct_ldsc_vs_gwaslab): coverage 0.954 to
    0.999, allele_mismatch_fraction at most 0.0007, and position_concordance at least 0.99999.
    """

    min_coverage: float = 0.8
    max_allele_mismatch_fraction: float = 0.01
    min_position_concordance: float = 0.99

    def __attrs_post_init__(self) -> None:
        for name, value in (
            ("min_coverage", self.min_coverage),
            ("max_allele_mismatch_fraction", self.max_allele_mismatch_fraction),
            ("min_position_concordance", self.min_position_concordance),
        ):
            assert 0.0 <= value <= 1.0, f"{name} must lie in [0, 1], got {value}"


@frozen(slots=True)
class AlignmentDiagnostics:
    """Counts from aligning one trait to the context SNP set, plus the trait's declared build."""

    n_context: int
    n_found: int
    n_matched: int
    n_position_compared: int | None
    n_position_concordant: int | None
    declared_build: GenomeBuild | None

    def __attrs_post_init__(self) -> None:
        assert 0 <= self.n_matched <= self.n_found <= self.n_context
        assert (self.n_position_compared is None) == (
            self.n_position_concordant is None
        ), "position counts must be both set or both None"

    @property
    def coverage(self) -> float:
        return self.n_matched / self.n_context if self.n_context else 0.0

    @property
    def allele_mismatch_fraction(self) -> float:
        if self.n_found == 0:
            return 0.0
        return (self.n_found - self.n_matched) / self.n_found

    @property
    def position_concordance(self) -> float | None:
        if (
            self.n_position_compared is None
            or self.n_position_concordant is None
            or self.n_position_compared == 0
        ):
            return None
        return self.n_position_concordant / self.n_position_compared

    def to_json_dict(self) -> dict[str, int | str | None]:
        return {
            "n_context": self.n_context,
            "n_found": self.n_found,
            "n_matched": self.n_matched,
            "n_position_compared": self.n_position_compared,
            "n_position_concordant": self.n_position_concordant,
            "declared_build": self.declared_build,
        }

    @classmethod
    def from_json_dict(cls, d: dict[str, int | str | None]) -> AlignmentDiagnostics:
        return cls(
            n_context=_as_int(d["n_context"]),
            n_found=_as_int(d["n_found"]),
            n_matched=_as_int(d["n_matched"]),
            n_position_compared=_as_optional_int(d["n_position_compared"]),
            n_position_concordant=_as_optional_int(d["n_position_concordant"]),
            declared_build=_as_optional_build(d["declared_build"]),
        )


def _as_optional_build(value: int | str | None) -> GenomeBuild | None:
    if value == "19":
        return "19"
    if value == "38":
        return "38"
    assert value is None, f"expected a genome build or None, got {value!r}"
    return None


def _as_int(value: int | str | None) -> int:
    assert isinstance(value, int), f"expected an int, got {value!r}"
    return value


def _as_optional_int(value: int | str | None) -> int | None:
    return None if value is None else _as_int(value)


def alignment_diagnostics(
    aligned: TraitAligned, declared_build: GenomeBuild | None
) -> AlignmentDiagnostics:
    return AlignmentDiagnostics(
        n_context=int(aligned.z.shape[0]),
        n_found=aligned.n_found,
        n_matched=aligned.n_matched,
        n_position_compared=aligned.n_position_compared,
        n_position_concordant=aligned.n_position_concordant,
        declared_build=declared_build,
    )


def check_alignment_guards(
    diagnostics: AlignmentDiagnostics, guards: AlignmentGuards, alias: str
) -> None:
    """Assert, in order, the coverage, allele-mismatch and (build 19 only) position guards."""
    assert diagnostics.coverage >= guards.min_coverage, (
        f"{alias}: coverage {diagnostics.coverage:.4f} of the HapMap3 regression SNPs is below "
        f"min_coverage={guards.min_coverage}; check the trait's rsID assignment"
    )
    assert (
        diagnostics.allele_mismatch_fraction <= guards.max_allele_mismatch_fraction
    ), (
        f"{alias}: allele_mismatch_fraction {diagnostics.allele_mismatch_fraction:.4f} exceeds "
        f"max_allele_mismatch_fraction={guards.max_allele_mismatch_fraction}; the trait's rsIDs "
        "likely point at the wrong variants (e.g. assigned from the wrong genome build)"
    )
    concordance = diagnostics.position_concordance
    if diagnostics.declared_build == "19" and concordance is not None:
        assert concordance >= guards.min_position_concordance, (
            f"{alias}: position_concordance {concordance:.4f} with the hg19 index is below "
            f"min_position_concordance={guards.min_position_concordance}, although the trait "
            "declares build 19; the declared build is likely wrong"
        )
