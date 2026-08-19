"""
How much on-disk storage does storing per-variant N cost in a CSF pQTL slim file?

Western et al. 2024 CSF summary statistics carry a per-variant sample size (OBS_CT /
GWAS-SSF `n`) that varies across variants, unlike UKB-PPP where N is a single constant per
protein and is recovered separately by a ranged read. That leaves a design choice for the
CSF database: store N in the slim per-aptamer parquet, or reconstruct/approximate it.

This probe measures the storage side of that choice on one real aptamer, at both scales
that matter:

  - FULL: all variants in the published file (the whole 7.3M-row source).
  - ALIGNED: the HapMap3 subset (~1.0M rows), which is what a slim file aligned to a
    HapMap3-mode variant index would actually store.

and under four encodings, because byte-stream-split is tuned for the high-entropy mantissas
of beta/se and is a poor fit for a column with only a few hundred distinct values:

  A  beta, se                      (float32, byte-stream-split)
  B  beta, se, N                   (all float32, all byte-stream-split)
  C  beta, se BSS + N float32      (N dictionary-encoded instead of split)
  D  beta, se BSS + N int32        (N dictionary-encoded, integer typed)

Sizes are reported both as whole files and as exact per-column compressed byte counts read
back out of the parquet metadata, so the cost of N is measured directly rather than inferred
by differencing.

Run: pixi r python experiments/claude/csf_pqtl_database/n_storage_cost_probe.py
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from pathlib import Path

import gwaslab as gl
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from mecfs_bio.build_system.task.dataframe_output import write_parquet_table  # noqa: E402
from mecfs_bio.constants.gwaslab_constants import (  # noqa: E402
    GWASLAB_HAPMAP3_HG38_SNPLIST_RELPATH,
)

# One arbitrary CSF aptamer, from the GWAS Catalog deposit of PMID 39528825.
ACCESSION = "GCST90421540"
SOURCE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
    f"GCST90421001-GCST90422000/{ACCESSION}/{ACCESSION}.tsv.gz"
)
SOURCE_MD5 = "047befd46b553da2bcecf7c8faa91749"

# GWAS-SSF v1.0 column names.
SSF_CHROM = "chromosome"
SSF_POS = "base_pair_location"
SSF_EA = "effect_allele"
SSF_OA = "other_allele"
SSF_BETA = "beta"
SSF_SE = "standard_error"
SSF_N = "n"

# Max zstd level, per the user's request; pyarrow's ceiling for zstd is 22.
ZSTD_LEVEL = 22

N_APTAMERS = 7008  # published CSF aptamers, for extrapolation


def cached_download(dest: Path) -> Path:
    """Download the source file unless a copy with the published md5 is already present."""
    if dest.exists() and file_md5(dest) == SOURCE_MD5:
        print(f"using cached {dest} ({dest.stat().st_size:,} bytes)")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"downloading {SOURCE_URL}")
    subprocess.run(["curl", "-sSL", SOURCE_URL, "-o", str(dest)], check=True)
    actual = file_md5(dest)
    assert actual == SOURCE_MD5, f"md5 mismatch: got {actual}, expected {SOURCE_MD5}"
    print(f"downloaded {dest.stat().st_size:,} bytes, md5 verified")
    return dest


def file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_sumstats(path: Path) -> pl.DataFrame:
    """Read the alignment columns plus beta/se/N from a GWAS-SSF file."""
    return pl.read_csv(
        path,
        separator="\t",
        columns=[SSF_CHROM, SSF_POS, SSF_EA, SSF_OA, SSF_BETA, SSF_SE, SSF_N],
        schema_overrides={SSF_CHROM: pl.Int32, SSF_POS: pl.Int32, SSF_N: pl.Int32},
    )


def hapmap3_aligned(sumstats: pl.DataFrame) -> pl.DataFrame:
    """The HapMap3 subset of the sumstats, in (chrom, pos) order.

    Matches the gwaslab hg38 HapMap3 snplist allele-aware on the unordered allele set, the
    same rule ConstructPppVariantIndexTask uses, so the row count is what a HapMap3-mode CSF
    index would hold.
    """
    snplist = Path(os.path.dirname(gl.__file__)) / GWASLAB_HAPMAP3_HG38_SNPLIST_RELPATH
    hm3 = pl.read_csv(snplist, separator="\t").select(
        pl.col("#CHROM").cast(pl.Int32).alias(SSF_CHROM),
        pl.col("POS").cast(pl.Int32).alias(SSF_POS),
        allele_key("A1", "A2").alias("key"),
    )
    return (
        sumstats.with_columns(allele_key(SSF_EA, SSF_OA).alias("key"))
        .join(hm3, on=[SSF_CHROM, SSF_POS, "key"], how="semi")
        .sort([SSF_CHROM, SSF_POS])
    )


def allele_key(first: str, second: str) -> pl.Expr:
    """Order-agnostic allele-set key, so a swapped reference orientation still matches."""
    return pl.concat_str(
        [
            pl.min_horizontal(pl.col(first), pl.col(second)),
            pl.max_horizontal(pl.col(first), pl.col(second)),
        ],
        separator="_",
    )


def build_variants(frame: pl.DataFrame) -> dict[str, tuple[pa.Table, list[str]]]:
    """The four (table, byte-stream-split columns) encodings under comparison."""
    beta = pl.col(SSF_BETA).cast(pl.Float32)
    se = pl.col(SSF_SE).cast(pl.Float32)
    beta_se = frame.select(beta, se).to_arrow()
    with_n_float = frame.select(beta, se, pl.col(SSF_N).cast(pl.Float32)).to_arrow()
    with_n_int = frame.select(beta, se, pl.col(SSF_N).cast(pl.Int32)).to_arrow()
    return {
        "A_beta_se": (beta_se, [SSF_BETA, SSF_SE]),
        "B_beta_se_n_f32_bss": (with_n_float, [SSF_BETA, SSF_SE, SSF_N]),
        "C_beta_se_n_f32_dict": (with_n_float, [SSF_BETA, SSF_SE]),
        "D_beta_se_n_i32_dict": (with_n_int, [SSF_BETA, SSF_SE]),
    }


def column_bytes(path: Path) -> dict[str, int]:
    """Compressed bytes per column, summed over row groups, from the parquet metadata."""
    metadata = pq.ParquetFile(path).metadata
    totals: dict[str, int] = {}
    for group in range(metadata.num_row_groups):
        row_group = metadata.row_group(group)
        for column in range(row_group.num_columns):
            chunk = row_group.column(column)
            name = chunk.path_in_schema
            totals[name] = totals.get(name, 0) + chunk.total_compressed_size
    return totals


def run_scale(label: str, frame: pl.DataFrame, out_dir: Path) -> None:
    print(f"\n{'=' * 78}\n{label}: {frame.height:,} rows\n{'=' * 78}")
    baseline: int | None = None
    for name, (table, split_columns) in build_variants(frame).items():
        out_path = out_dir / f"{label}_{name}.parquet.zstd"
        write_parquet_table(
            table=table,
            out_path=out_path,
            compression="zstd",
            compression_level=ZSTD_LEVEL,
            byte_stream_split_columns=split_columns,
        )
        size = out_path.stat().st_size
        if baseline is None:
            baseline = size
        per_column = column_bytes(out_path)
        delta = size - baseline
        print(
            f"\n{name:24} {size:>12,} B  ({size / 1e6:6.2f} MB)"
            f"   delta vs A: {delta:+,} B ({100 * delta / baseline:+.1f}%)"
        )
        for column, byte_count in per_column.items():
            print(f"    {column:20} {byte_count:>12,} B")
        print(
            f"    extrapolated to {N_APTAMERS:,} aptamers: "
            f"{size * N_APTAMERS / 1e9:.1f} GB"
        )


def main() -> None:
    scratch = Path(
        os.environ.get("CSF_PROBE_SCRATCH", Path.home() / ".cache" / "csf_pqtl_probe")
    )
    scratch.mkdir(parents=True, exist_ok=True)
    source = cached_download(scratch / f"{ACCESSION}.tsv.gz")

    sumstats = read_sumstats(source)
    print(f"\nsource: {ACCESSION}, {sumstats.height:,} variants")
    n_series = sumstats[SSF_N]
    print(
        f"per-variant N: min={n_series.min()} max={n_series.max()} "
        f"median={n_series.median()} distinct={n_series.n_unique()}"
    )
    print(f"N is constant across variants: {n_series.n_unique() == 1}")

    run_scale("full", sumstats, scratch)
    run_scale("hapmap3_aligned", hapmap3_aligned(sumstats), scratch)


if __name__ == "__main__":
    main()
