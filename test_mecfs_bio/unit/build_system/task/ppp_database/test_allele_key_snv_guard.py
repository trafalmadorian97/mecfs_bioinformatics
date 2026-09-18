import polars as pl
import pytest

from mecfs_bio.build_system.task.ppp_database.allele_key import assert_all_snv


def test_assert_all_snv_passes_for_snvs():
    assert_all_snv(pl.DataFrame({"EA": ["A", "C"], "NEA": ["G", "T"]}), "EA", "NEA")


def test_assert_all_snv_rejects_indels():
    df = pl.DataFrame({"EA": ["A", "TCA"], "NEA": ["G", "T"]})
    with pytest.raises(AssertionError):
        assert_all_snv(df, "EA", "NEA")
