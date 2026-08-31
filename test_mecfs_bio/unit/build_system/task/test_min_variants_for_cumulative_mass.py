import polars as pl

from mecfs_bio.build_system.task.pipes.min_variants_for_cumulative_mass import (
    MinVariantsForCumulativeMass,
)


def test_keeps_minimal_prefix_per_group_crossing_threshold():
    df = pl.DataFrame(
        {
            "cs": ["L1", "L1", "L1", "L2", "L2"],
            "pip": [0.4, 0.3, 0.2, 0.6, 0.4],
            "v": ["a", "b", "c", "d", "e"],
        }
    )
    out = MinVariantsForCumulativeMass(
        group_col="cs", value_col="pip", threshold=0.5
    ).process_eager_polars(df)
    # L1: 0.4 then 0.7 crosses 0.5 -> keep a, b (drop tail c). L2: 0.6 alone -> d.
    assert set(out["v"].to_list()) == {"a", "b", "d"}


def test_group_never_reaching_threshold_keeps_all_rows():
    df = pl.DataFrame({"cs": ["L1", "L1"], "pip": [0.2, 0.15], "v": ["a", "b"]})
    out = MinVariantsForCumulativeMass(
        group_col="cs", value_col="pip", threshold=0.5
    ).process_eager_polars(df)
    assert set(out["v"].to_list()) == {"a", "b"}
