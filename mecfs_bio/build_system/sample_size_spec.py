"""
Specification of how the sample size (N) of a GWAS should be supplied to
downstream analyses.

Most GWAS have a single scalar sample size, but meta-analyses and synthetic GWAS
produced by GenomicSEM have a sample size that varies per variant. The two cases
are modelled by :class:`ScalarSampleSize` and :class:`PerVariantSampleSize`
respectively, so that call sites and tasks can branch on a named type rather than
passing a bare ``int`` and a separate flag.
"""

from attrs import frozen

from mecfs_bio.constants.gwaslab_constants import GWASLAB_SAMPLE_SIZE_COLUMN


@frozen
class ScalarSampleSize:
    """A single sample size shared by every variant in the GWAS."""

    n: int


@frozen
class PerVariantSampleSize:
    """The sample size varies per variant and is read from ``column``.

    The column must already be present on the GWAS data (e.g. the per-variant
    effective sample size of a meta-analysis) and survive sumstats creation.
    """

    column: str = GWASLAB_SAMPLE_SIZE_COLUMN


SampleSizeSpec = ScalarSampleSize | PerVariantSampleSize


def coerce_sample_size(value: "int | SampleSizeSpec") -> SampleSizeSpec:
    """Wrap a bare ``int`` as a :class:`ScalarSampleSize` for backwards
    compatibility with call sites that pass a scalar sample size directly."""
    if isinstance(value, int):
        return ScalarSampleSize(value)
    return value


def per_variant_column(value: SampleSizeSpec) -> str | None:
    """The data column holding the per-variant sample size, or ``None`` for a
    scalar sample size (so callers know not to write an N column)."""
    if isinstance(value, PerVariantSampleSize):
        return value.column
    return None
