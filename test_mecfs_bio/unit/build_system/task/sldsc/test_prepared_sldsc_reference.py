from pathlib import Path

import numpy as np
import polars as pl
import pytest

from mecfs_bio.build_system.task.sldsc.prepared_sldsc_reference import (
    BP_COL,
    CHR_COL,
    SNP_COL,
    W_LD_COL,
    LdScoreTable,
    PreparedAnnotationSet,
    PreparedBaseline,
    SLDSCAnnotationModel,
    read_prepared_annotation_set,
    read_prepared_baseline,
    write_prepared_annotation_set,
    write_prepared_baseline,
)


def _table(n_cols: int) -> LdScoreTable:
    snps = pl.DataFrame(
        {SNP_COL: ["rs1", "rs2", "rs3"], CHR_COL: [1, 1, 2], BP_COL: [10, 20, 5]}
    )
    return LdScoreTable(
        snps=snps,
        column_ids=tuple(f"ld_{i}" for i in range(n_cols)),
        sources=tuple(f"src:{i}" for i in range(n_cols)),
        matrix=np.arange(3 * n_cols, dtype=np.float32).reshape(3, n_cols),
        m_5_50=np.arange(1, n_cols + 1, dtype=np.float64),
    )


def test_annotation_set_round_trips(tmp_path: Path) -> None:
    original = PreparedAnnotationSet(
        ld=_table(3),
        models=(
            SLDSCAnnotationModel(
                label="a", target_column="ld_0", control_columns=("ld_2",)
            ),
            SLDSCAnnotationModel(
                label="b", target_column="ld_1", control_columns=("ld_2",)
            ),
        ),
    )
    write_prepared_annotation_set(tmp_path, original)
    loaded = read_prepared_annotation_set(tmp_path)
    assert loaded.models == original.models
    assert loaded.ld.column_ids == original.ld.column_ids
    np.testing.assert_array_equal(loaded.ld.matrix, original.ld.matrix)
    np.testing.assert_array_equal(loaded.ld.m_5_50, original.ld.m_5_50)
    assert loaded.ld.snps.equals(original.ld.snps)


def test_baseline_round_trips(tmp_path: Path) -> None:
    w_ld = pl.DataFrame({SNP_COL: ["rs1", "rs2", "rs3"], W_LD_COL: [1.0, 2.0, 3.0]})
    original = PreparedBaseline(ld=_table(2), w_ld=w_ld)
    write_prepared_baseline(tmp_path, original)
    loaded = read_prepared_baseline(tmp_path)
    assert loaded.w_ld.equals(original.w_ld)
    np.testing.assert_array_equal(loaded.ld.matrix, original.ld.matrix)


def test_baseline_with_misaligned_w_ld_is_rejected() -> None:
    w_ld = pl.DataFrame({SNP_COL: ["rs2", "rs1", "rs3"], W_LD_COL: [1.0, 2.0, 3.0]})
    with pytest.raises(AssertionError):
        PreparedBaseline(ld=_table(2), w_ld=w_ld)


def test_model_referencing_missing_column_is_rejected() -> None:
    with pytest.raises(AssertionError):
        PreparedAnnotationSet(
            ld=_table(2),
            models=(
                SLDSCAnnotationModel(
                    label="a", target_column="ld_0", control_columns=("ld_9",)
                ),
            ),
        )


def test_target_among_its_controls_is_rejected() -> None:
    with pytest.raises(AssertionError):
        SLDSCAnnotationModel(label="a", target_column="ld_0", control_columns=("ld_0",))
