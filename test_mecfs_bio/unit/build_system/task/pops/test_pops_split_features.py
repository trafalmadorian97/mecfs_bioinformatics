"""
Unit tests for the streaming POPs feature-file column splitter.
"""

from pathlib import Path

import pandas as pd
import pytest

from mecfs_bio.build_system.task.pops.pops_split_features_task import (
    CHUNK_FILENAME_TEMPLATE,
    split_feature_file,
)


def _write_table(path: Path, num_genes: int, num_features: int) -> pd.DataFrame:
    """Write a small synthetic tab-separated feature table and return it."""
    data = {"ENSGID": [f"ENSG{i:08d}" for i in range(num_genes)]}
    for j in range(num_features):
        data[f"feat.{j}"] = [float(i * 100 + j) for i in range(num_genes)]
    df = pd.DataFrame(data)
    df.to_csv(path, sep="\t", index=False)
    return df


def _read_chunks(out_dir: Path, num_chunks: int) -> list[pd.DataFrame]:
    return [
        pd.read_csv(
            out_dir / CHUNK_FILENAME_TEMPLATE.format(i), sep="\t", index_col="ENSGID"
        )
        for i in range(num_chunks)
    ]


def test_split_preserves_columns_and_values(tmp_path: Path):
    """Concatenating the chunks column-wise reproduces the original table exactly."""
    source = tmp_path / "features.txt"
    original = _write_table(source, num_genes=6, num_features=7)
    out_dir = tmp_path / "chunks"
    out_dir.mkdir()

    num_chunks = split_feature_file(source, out_dir, columns_per_file=3)

    # 7 feature columns in groups of 3 -> 3 chunks.
    assert num_chunks == 3
    chunks = _read_chunks(out_dir, num_chunks)
    # Every chunk carries the ENSGID index and no more than columns_per_file features.
    assert all(chunk.shape[1] <= 3 for chunk in chunks)
    reconstructed = pd.concat(chunks, axis=1)
    expected = original.set_index("ENSGID")
    pd.testing.assert_frame_equal(reconstructed, expected)


def test_split_exact_multiple_of_chunk_size(tmp_path: Path):
    """When features divide evenly, no trailing partial chunk is produced."""
    source = tmp_path / "features.txt"
    _write_table(source, num_genes=4, num_features=6)
    out_dir = tmp_path / "chunks"
    out_dir.mkdir()

    num_chunks = split_feature_file(source, out_dir, columns_per_file=2)

    assert num_chunks == 3
    assert not (out_dir / CHUNK_FILENAME_TEMPLATE.format(3)).exists()


def test_split_rejects_wrong_first_column(tmp_path: Path):
    """A source whose first column is not ENSGID fails fast."""
    source = tmp_path / "bad.txt"
    source.write_text("GENE\tfeat.0\nENSG00000000001\t1.0\n", encoding="utf-8")
    out_dir = tmp_path / "chunks"
    out_dir.mkdir()

    with pytest.raises(AssertionError, match="ENSGID"):
        split_feature_file(source, out_dir, columns_per_file=2)


def test_split_rejects_ragged_row(tmp_path: Path):
    """A data row with the wrong number of fields fails fast."""
    source = tmp_path / "ragged.txt"
    source.write_text(
        "ENSGID\tfeat.0\tfeat.1\nENSG00000000001\t1.0\n", encoding="utf-8"
    )
    out_dir = tmp_path / "chunks"
    out_dir.mkdir()

    with pytest.raises(AssertionError, match="fields"):
        split_feature_file(source, out_dir, columns_per_file=2)
