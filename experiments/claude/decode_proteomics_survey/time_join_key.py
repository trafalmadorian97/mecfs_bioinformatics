"""
Join key cost: the index (NAME, CHR, POS, EA, NEA) against the AP4M1 rows at index positions
(about 1.4M candidates, as one_pass keeps), joined on the Name string versus on the 4-column key
(chromosome code, Pos, effectAllele, otherAllele). Both left joins validated one-to-one; five
timed repeats each; checks both give the same matches.

  pixi r python experiments/claude/decode_proteomics_survey/time_join_key.py
"""

import time

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_index import (
    HAPMAP_3_DECODE_SMP_INDEX,
)
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_reference_protein_file import (
    DECODE_SMP_AP4M1_RAW_FILE,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.decode_database.decode_sumstats_format import (
    chrom_code_expr,
)

REPEATS = 5


def _time(label: str, run) -> pl.DataFrame:  # noqa: ANN001
    out = run()
    times = []
    for _ in range(REPEATS):
        t0 = time.perf_counter()
        out = run()
        times.append(time.perf_counter() - t0)
    print(f"{label}: min {min(times):.3f}s median {sorted(times)[REPEATS // 2]:.3f}s")
    return out


def main() -> None:
    assets = DEFAULT_RUNNER.run([DECODE_SMP_AP4M1_RAW_FILE, HAPMAP_3_DECODE_SMP_INDEX])
    raw = assets[DECODE_SMP_AP4M1_RAW_FILE.asset_id]
    index_asset = assets[HAPMAP_3_DECODE_SMP_INDEX.asset_id]
    assert isinstance(raw, FileAsset) and isinstance(index_asset, FileAsset)
    index = pl.read_parquet(index_asset.path, columns=["NAME", "CHR", "POS", "EA", "NEA"])
    rows = pl.read_csv(
        raw.path,
        separator="\t",
        columns=["Chrom", "Pos", "Name", "effectAllele", "otherAllele", "Beta", "SE", "N"],
        null_values=["."],
        quote_char=None,
        schema_overrides={"Pos": pl.Int64},
    )
    candidates = rows.filter(pl.col("Pos").is_in(index["POS"].cast(pl.Int64).implode()))
    candidates = candidates.with_columns(chrom_code_expr("Chrom").cast(pl.Int8).alias("CHR"))
    index = index.with_columns(pl.col("CHR").cast(pl.Int8), pl.col("POS").cast(pl.Int64))
    print(f"index {index.height}, candidates {candidates.height}")

    def by_name() -> pl.DataFrame:
        matched = candidates.filter(pl.col("Name").is_in(index["NAME"].implode()))
        return index.select("NAME").join(
            matched.select("Name", "Beta", "SE", "N"),
            left_on="NAME",
            right_on="Name",
            how="left",
            validate="1:1",
            maintain_order="left",
        )

    def by_key() -> pl.DataFrame:
        right = candidates.select(
            "CHR",
            pl.col("Pos").alias("POS"),
            pl.col("effectAllele").alias("EA"),
            pl.col("otherAllele").alias("NEA"),
            "Beta",
            "SE",
            "N",
        )
        matched = right.join(index.select("CHR", "POS", "EA", "NEA"), on=["CHR", "POS", "EA", "NEA"], how="semi")
        return index.select("CHR", "POS", "EA", "NEA").join(
            matched,
            on=["CHR", "POS", "EA", "NEA"],
            how="left",
            validate="1:1",
            maintain_order="left",
        )

    a = _time("join on Name", by_name)
    b = _time("join on CHR, POS, EA, NEA", by_key)
    assert a.select("Beta", "SE", "N").equals(b.select("Beta", "SE", "N"))
    print("same result")


if __name__ == "__main__":
    main()
