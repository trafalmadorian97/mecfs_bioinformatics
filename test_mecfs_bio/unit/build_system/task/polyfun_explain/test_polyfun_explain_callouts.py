import polars as pl

from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    FAMILY_COL,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    FAMILY_CONTRAST_COL,
    _callout_families,
    _format_callout_label,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)


def _key(pos: int) -> dict:
    return {
        GWASLAB_CHROM_COL: 1,
        GWASLAB_POS_COL: pos,
        GWASLAB_EFFECT_ALLELE_COL: "T",
        GWASLAB_NON_EFFECT_ALLELE_COL: "A",
    }


def test_callout_families_buckets_and_orders_by_z():
    focal = _key(123)
    per_family = pl.DataFrame(
        [
            {**focal, FAMILY_COL: "conserved", FAMILY_CONTRAST_COL: 3.0},  # z=3 -> ++
            {**focal, FAMILY_COL: "coding", FAMILY_CONTRAST_COL: 1.5},  # z=1.5 -> +
            {
                **focal,
                FAMILY_COL: "repressed",
                FAMILY_CONTRAST_COL: 0.5,
            },  # z=0.5 -> drop
            {
                **focal,
                FAMILY_COL: "histone_marks",
                FAMILY_CONTRAST_COL: -4.0,
            },  # neg -> drop
        ]
    )
    family_sd = {
        "conserved": 1.0,
        "coding": 1.0,
        "repressed": 1.0,
        "histone_marks": 1.0,
    }
    result = _callout_families(per_family, focal, family_sd, max_families=3)
    assert result == [("conserved", "++"), ("coding", "+")]


def test_callout_families_skips_degenerate_sd():
    focal = _key(123)
    per_family = pl.DataFrame(
        [{**focal, FAMILY_COL: "conserved", FAMILY_CONTRAST_COL: 3.0}]
    )
    assert _callout_families(per_family, focal, {"conserved": 0.0}, 3) == []


def test_format_label_with_and_without_families():
    focal = _key(174128548)
    assert (
        _format_callout_label(focal, [("conserved", "++"), ("coding", "+")])
        == "174128548:A:T (conserved ++, coding +)"
    )
    assert _format_callout_label(focal, []) == "174128548:A:T"
