import gzip
import math
from pathlib import Path

import polars as pl
import pyarrow.parquet as pq
import pytest

from mecfs_bio.build_system.task.csf_database.build_slim_aptamer_parquet_task import (
    CsfAptamerFile,
    align_aptamer_to_index,
    read_aptamer_sumstats,
    write_slim_aptamer_parquet,
)
from mecfs_bio.constants.csf_database_constants import (
    Analyte,
    GcstAccession,
    SeqId,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

_SSF_HEADER = (
    "chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\t"
    "standard_error\teffect_allele_frequency\tneg_log_10_p_value\tvariant_id\t"
    "rs_id\tn"
)


def test_csf_aptamer_file_sumstats_url():
    aptamer = CsfAptamerFile(
        analyte=Analyte("X13681.173"),
        seq_id=SeqId("13681-173"),
        accession=GcstAccession("GCST90421540"),
        entrez_gene_symbol="CSNK2A2",
    )
    assert aptamer.sumstats_url.endswith(
        "GCST90421001-GCST90422000/GCST90421540/GCST90421540.tsv.gz"
    )


def _index(chrom, pos, ea, nea):
    return pl.DataFrame(
        {
            GWASLAB_CHROM_COL: chrom,
            GWASLAB_POS_COL: pos,
            GWASLAB_NON_EFFECT_ALLELE_COL: nea,
            GWASLAB_EFFECT_ALLELE_COL: ea,
        }
    )


def test_align_aptamer_to_index():
    # Index deliberately NOT position-sorted, to prove output follows index row order.
    index = _index([2, 1, 1], [300, 100, 200], ["A", "A", "C"], ["T", "G", "T"])
    aptamer = pl.DataFrame(
        {
            "chromosome": [1, 1, 5],
            "base_pair_location": [100, 200, 999],
            "effect_allele": ["A", "T", "G"],  # row 2 swapped vs index (index EA=C)
            "other_allele": ["G", "C", "A"],
            "beta": [0.5, 0.3, 0.9],
            "standard_error": [0.1, 0.2, 0.4],
            "n": [3400, 3200, 3100],
        }
    )

    out = align_aptamer_to_index(index, aptamer)
    assert out.columns == [
        GWASLAB_BETA_COL,
        GWASLAB_SE_COL,
        GWASLAB_SAMPLE_SIZE_COLUMN,
    ]
    beta = out[GWASLAB_BETA_COL].to_list()
    se = out[GWASLAB_SE_COL].to_list()
    n = out[GWASLAB_SAMPLE_SIZE_COLUMN].to_list()

    # Row 0: index chr2:300 absent from aptamer -> NaN everywhere.
    assert math.isnan(beta[0]) and math.isnan(se[0]) and math.isnan(n[0])
    # Row 1: chr1:100 same orientation -> beta unchanged, per-variant N carried.
    assert beta[1] == pytest.approx(0.5, abs=1e-6)
    assert se[1] == pytest.approx(0.1, abs=1e-6)
    assert n[1] == pytest.approx(3400.0)
    # Row 2: chr1:200 swapped orientation -> beta flips, se and N unaffected.
    assert beta[2] == pytest.approx(-0.3, abs=1e-6)
    assert se[2] == pytest.approx(0.2, abs=1e-6)
    assert n[2] == pytest.approx(3200.0)


def test_write_slim_aptamer_parquet(tmp_path: Path):
    index = _index(
        [1, 1, 1, 2], [100, 200, 250, 300], ["A", "C", "G", "A"], ["G", "T", "C", "T"]
    )
    aptamer = pl.DataFrame(
        {
            "chromosome": [1, 1, 2],
            "base_pair_location": [100, 200, 300],
            "effect_allele": ["A", "T", "A"],  # chr1:200 swapped vs index EA=C
            "other_allele": ["G", "C", "T"],
            "beta": [0.5, 0.3, 0.7],
            "standard_error": [0.1, 0.2, 0.3],
            "n": [3400, 3300, 3200],
        }
    )

    out_path = tmp_path / "slim.parquet"
    write_slim_aptamer_parquet(aptamer, index, out_path)

    out = pl.read_parquet(out_path)
    assert out.columns == [
        GWASLAB_BETA_COL,
        GWASLAB_SE_COL,
        GWASLAB_SAMPLE_SIZE_COLUMN,
    ]
    assert out.height == 4
    beta = out[GWASLAB_BETA_COL].to_list()
    n = out[GWASLAB_SAMPLE_SIZE_COLUMN].to_list()
    assert beta[0] == pytest.approx(0.5, abs=1e-6)
    assert beta[1] == pytest.approx(-0.3, abs=1e-6)  # swapped -> flipped
    assert math.isnan(beta[2])  # chr1:250 absent from aptamer
    assert math.isnan(n[2])
    assert beta[3] == pytest.approx(0.7, abs=1e-6)
    assert n[3] == pytest.approx(3200.0)

    # The output must actually use byte-stream-split on all three columns.
    for col in range(3):
        encodings = pq.ParquetFile(out_path).metadata.row_group(0).column(col).encodings
        assert "BYTE_STREAM_SPLIT" in encodings


def test_read_aptamer_sumstats_ignores_extra_columns(tmp_path: Path):
    # The real files carry EAF / -log10p / variant_id / rs_id too; the reader must
    # select only the alignment columns and tolerate a gzip file with no .gz suffix.
    gz_path = tmp_path / "aptamer_no_extension"
    with gzip.open(gz_path, "wt") as handle:
        handle.write(_SSF_HEADER + "\n")
        handle.write("1\t100\tA\tG\t0.5\t0.1\t0.2\t1.3\t1_100_G_A\trs1\t3400\n")

    frame = read_aptamer_sumstats(gz_path)
    assert set(frame.columns) == {
        "chromosome",
        "base_pair_location",
        "effect_allele",
        "other_allele",
        "beta",
        "standard_error",
        "n",
    }
    assert frame["n"].to_list() == [3400]
