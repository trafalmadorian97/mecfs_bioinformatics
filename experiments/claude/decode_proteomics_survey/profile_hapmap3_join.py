"""
Compare ways to keep HapMap3 rows of one deCODE file's parsed blocks (16 MiB line blocks,
decompressed with ISA-L; decompression and parsing excluded from timings):

  string_pair_join   mockup version: CHR/POS + unordered allele strings, semi-join on 33M rows
  pos_first_then_check semi-join on Pos alone (integer, no string work), then the full
                     CHR/POS/allele check on the survivors
  int_key_then_check semi-join on one Int64 key (CHR * 2^32 + POS), then check the unordered
                     allele pair only on the ~1.2M surviving rows (join to hm3 alleles)

Usage: pixi exec --spec python=3.12 --spec polars --spec python-isal -- \
    python -I profile_hapmap3_join.py <local.txt.gz> <hapmap3.parquet>
"""

import io
import sys
import time
from pathlib import Path

import polars as pl
from isal import isal_zlib

MIB = 1 << 20
HEADER = ["Chrom", "Pos", "Name", "rsids", "effectAllele", "otherAllele", "Beta", "Pval",
          "minus_log10_pval", "SE", "N", "ImpMAF"]
TYPES = [pl.String, pl.Int64, pl.String, pl.String, pl.String, pl.String, pl.Float64,
         pl.Float64, pl.Float64, pl.Float64, pl.Float64, pl.Float64]
SCHEMA = pl.Schema(dict(zip(HEADER, TYPES, strict=True)))
KEEP = ["Chrom", "Pos", "Name", "effectAllele", "otherAllele", "Beta", "SE", "N", "ImpMAF"]


def blocks(path: Path):
    d = isal_zlib.decompressobj(31)
    pending, first = b"", True
    with open(path, "rb") as f:
        while c := f.read(4 * MIB):
            pending += d.decompress(c)
            if len(pending) < 16 * MIB:
                continue
            cut = pending.rindex(b"\n") + 1
            block, pending = pending[:cut], pending[cut:]
            if first:
                block, first = block[block.index(b"\n") + 1:], False
            yield block
    pending += d.flush()
    if pending:
        yield pending


def pair(a, b):
    lo = pl.when(pl.col(a) <= pl.col(b)).then(pl.col(a)).otherwise(pl.col(b))
    hi = pl.when(pl.col(a) <= pl.col(b)).then(pl.col(b)).otherwise(pl.col(a))
    return [lo.alias("_a"), hi.alias("_b")]


def chrom_int(col):
    return pl.col(col).str.strip_prefix("chr").replace({"X": "23"}).cast(pl.Int64)


def string_pair_join(rows, hm3_pairs):
    keyed = rows.with_columns(chrom_int("Chrom").cast(pl.Int32).alias("CHR"),
                              pl.col("Pos").cast(pl.Int32).alias("POS"),
                              *pair("effectAllele", "otherAllele"))
    return keyed.join(hm3_pairs, on=["CHR", "POS", "_a", "_b"], how="semi").drop("CHR", "POS", "_a", "_b")


def int_key_then_check(rows, hm3_keyed):
    key = (chrom_int("Chrom") * (1 << 32) + pl.col("Pos")).alias("_k")
    cand = rows.with_columns(key).join(hm3_keyed.select("_k").unique(), on="_k", how="semi")
    cand = cand.with_columns(*pair("effectAllele", "otherAllele"))
    return cand.join(hm3_keyed, on=["_k", "_a", "_b"], how="semi").drop("_k", "_a", "_b")


def pos_first_then_check(rows, hm3_pairs, hm3_pos):
    cand = rows.join(hm3_pos, left_on="Pos", right_on="POS", how="semi")
    return string_pair_join(cand, hm3_pairs)


def main():
    path, hm3_path = Path(sys.argv[1]), Path(sys.argv[2])
    hm3 = pl.read_parquet(hm3_path)
    hm3_pairs = hm3.select("CHR", "POS", *pair("EA", "NEA"))
    hm3_keyed = hm3.select((pl.col("CHR").cast(pl.Int64) * (1 << 32) + pl.col("POS")).alias("_k"),
                           *pair("EA", "NEA"))
    hm3_pos = hm3.select(pl.col("POS").cast(pl.Int64)).unique()
    t_a = t_b = t_c = 0.0
    out_a, out_b, out_c = [], [], []
    for block in blocks(path):
        rows = pl.read_csv(io.BytesIO(block), separator="\t", has_header=False, schema=SCHEMA,
                           null_values=["."]).select(KEEP)
        t = time.perf_counter(); out_a.append(string_pair_join(rows, hm3_pairs)); t_a += time.perf_counter() - t
        t = time.perf_counter(); out_b.append(int_key_then_check(rows, hm3_keyed)); t_b += time.perf_counter() - t
        t = time.perf_counter(); out_c.append(pos_first_then_check(rows, hm3_pairs, hm3_pos)); t_c += time.perf_counter() - t
    a, b = pl.concat(out_a), pl.concat(out_b)
    print(f"string_pair_join   {t_a:5.1f}s  rows={a.height:,}")
    print(f"int_key_then_check {t_b:5.1f}s  rows={b.height:,}  identical={a.equals(b)}")
    c = pl.concat(out_c)
    print(f"pos_first_then_check {t_c:5.1f}s  rows={c.height:,}  identical={a.equals(c)}")


if __name__ == "__main__":
    main()
