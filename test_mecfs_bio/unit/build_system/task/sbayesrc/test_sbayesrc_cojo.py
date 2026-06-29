"""Unit tests for the COJO .ma writer used by SBayesRC and polypwas."""

from pathlib import Path

import narwhals as nw
import polars as pl
import pytest

from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
    QuantPhenotype,
)
from mecfs_bio.build_system.task.sbayesrc.sbayesrc_cojo import (
    COJO_COLUMN_ORDER,
    write_cojo_ma,
)


def _gwaslab_frame() -> nw.LazyFrame:
    return nw.from_native(
        pl.DataFrame(
            {
                "rsID": ["rs1001", "rs1002", "rs1003"],
                "EA": ["A", "C", "A"],
                "NEA": ["G", "G", "C"],
                "EAF": [0.8493, 0.0306, 0.5128],
                "BETA": [0.0024, 0.0034, 0.0045],
                "SE": [0.0055, 0.0115, 0.0038],
                "P": [0.6653, 0.7659, 0.2319],
                "N": [129850, 129799, 129830],
            }
        )
    ).lazy()


def test_write_cojo_ma_quant_uses_total_sample_size(tmp_path: Path):
    out_path = tmp_path / "trait.ma"
    write_cojo_ma(_gwaslab_frame(), QuantPhenotype(total_sample_size=200000), out_path)

    result = pl.read_csv(out_path, separator="\t")
    assert result.columns == COJO_COLUMN_ORDER
    assert result["SNP"].to_list() == ["rs1001", "rs1002", "rs1003"]
    assert result["A1"].to_list() == ["A", "C", "A"]
    assert result["A2"].to_list() == ["G", "G", "C"]
    # N is filled from the phenotype total sample size, overriding the column.
    assert result["N"].to_list() == [200000, 200000, 200000]


def test_write_cojo_ma_quant_without_total_keeps_existing_n(tmp_path: Path):
    out_path = tmp_path / "trait.ma"
    write_cojo_ma(_gwaslab_frame(), QuantPhenotype(), out_path)

    result = pl.read_csv(out_path, separator="\t")
    assert result.columns == COJO_COLUMN_ORDER
    assert result["N"].to_list() == [129850, 129799, 129830]


def test_write_cojo_ma_quant_without_total_and_no_n_fails_fast(tmp_path: Path):
    frame = nw.from_native(
        pl.DataFrame(
            {
                "rsID": ["rs1"],
                "EA": ["A"],
                "NEA": ["G"],
                "EAF": [0.5],
                "BETA": [0.1],
                "SE": [0.05],
                "P": [0.1],
            }
        )
    ).lazy()
    with pytest.raises(AssertionError):
        write_cojo_ma(frame, QuantPhenotype(), tmp_path / "trait.ma")


def test_write_cojo_ma_binary_uses_effective_sample_size(tmp_path: Path):
    out_path = tmp_path / "trait.ma"
    phenotype = BinaryPhenotypeSampleInfo(
        sample_prevalence=0.5,
        estimated_population_prevalence=0.01,
        total_sample_size=100000,
    )
    write_cojo_ma(_gwaslab_frame(), phenotype, out_path)

    result = pl.read_csv(out_path, separator="\t")
    assert result["N"].unique().to_list() == [phenotype.effective_sample_size]
