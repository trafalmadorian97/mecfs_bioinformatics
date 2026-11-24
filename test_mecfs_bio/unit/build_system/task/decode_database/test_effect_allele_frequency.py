import polars as pl
import pytest

from mecfs_bio.build_system.task.decode_database.decode_sumstats_format import (
    chrom_code,
    chrom_code_expr,
)
from mecfs_bio.build_system.task.decode_database.effect_allele_frequency import (
    resolve_effect_allele_frequency,
)
from mecfs_bio.constants.decode_somascan_constants import (
    VI_A1_COL,
    VI_CHR_COL,
    VI_FREQ_COL,
    VI_POS_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_POS_COL,
)

_ROWS_SCHEMA = {
    GWASLAB_CHROM_COL: pl.Int8,
    GWASLAB_POS_COL: pl.Int64,
    GWASLAB_EFFECT_ALLELE_COL: pl.String,
    "ImpMAF": pl.Float64,
}
_VI_SCHEMA = {
    VI_CHR_COL: pl.Int8,
    VI_POS_COL: pl.Int64,
    VI_A1_COL: pl.String,
    VI_FREQ_COL: pl.Float64,
}

# One position per case: (pos, EA, ImpMAF, [(A1, f), ...], expected EAF).
_CASES = [
    (100, "A", 0.3, [("A", 0.7)], 0.7),  # single agreeing candidate
    (200, "C", 0.05, [("C", 0.2), ("C", 0.05)], 0.05),  # several, one agrees
    (300, "G", 0.4, [("G", 0.4), ("G", 0.4)], 0.4),  # several agree, equal f
    (400, "T", 0.3, [("T", 0.3), ("T", 0.7)], None),  # f and 1-f both agree
    (500, "A", 0.1, [("C", 0.1)], None),  # no candidate with A1 == EA
    (600, "A", 0.4, [("A", 0.2)], None),  # candidate disagrees
]


def _rows() -> pl.DataFrame:
    return pl.DataFrame(
        [(1, pos, ea, maf) for pos, ea, maf, _, _ in _CASES],
        schema=_ROWS_SCHEMA,
        orient="row",
    )


def _variant_info() -> pl.DataFrame:
    return pl.DataFrame(
        [(1, pos, a1, f) for pos, _, _, cands, _ in _CASES for a1, f in cands],
        schema=_VI_SCHEMA,
        orient="row",
    )


def test_eaf_rule_cases() -> None:
    eaf = resolve_effect_allele_frequency(_rows(), _variant_info())
    assert eaf.name == GWASLAB_EFFECT_ALLELE_FREQ_COL
    assert eaf.dtype == pl.Float64
    assert eaf.to_list() == pytest.approx([c[-1] for c in _CASES], nan_ok=False)


def test_eaf_preserves_row_order() -> None:
    rows = _rows().with_row_index("k").sample(fraction=1.0, shuffle=True, seed=3)
    eaf = resolve_effect_allele_frequency(rows.drop("k"), _variant_info())
    expected = [_CASES[k][-1] for k in rows["k"]]
    assert eaf.to_list() == pytest.approx(expected)


def test_chrom_code_text_order_and_x() -> None:
    text_order = ["chr1", "chr10", "chr2", "chr22", "chr9", "chrX"]
    assert [chrom_code(c) for c in text_order] == [1, 10, 2, 22, 9, 23]
    codes = pl.DataFrame({"Chrom": text_order}).select(chrom_code_expr("Chrom"))
    assert codes.to_series().to_list() == [1, 10, 2, 22, 9, 23]
    with pytest.raises(AssertionError):
        chrom_code("chrY")
