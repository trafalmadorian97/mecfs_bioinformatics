import narwhals
import polars as pl

from mecfs_bio.build_system.task.pipes.drop_palindromes_pipe import DropPalindromesPipe


def test_drops_only_strand_ambiguous_snvs():
    df = pl.DataFrame(
        {
            "EA": ["A", "C", "A", "G"],
            "NEA": ["T", "G", "G", "A"],  # rows 0 (A/T) and 1 (C/G) are palindromic
        }
    )
    out = (
        DropPalindromesPipe(ea_col="EA", nea_col="NEA")
        .process(narwhals.from_native(df).lazy())
        .collect()
        .to_native()
    )
    assert out.to_dicts() == [{"EA": "A", "NEA": "G"}, {"EA": "G", "NEA": "A"}]
