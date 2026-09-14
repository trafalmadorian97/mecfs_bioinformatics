"""Benchmark a vectorized "does column X match the FASTA reference at POS" check.

Implementation under test: memory-map the uncompressed .fa and translate (chrom, pos) to a file
offset with the .fai index (offset + p + p // line_bases * (line_bytes - line_bases)), then gather
bytes with numpy. Alleles are grouped by length so each group is one fancy-index gather + compare.
Nothing is loaded up front; the OS page cache serves the touched pages.

Measures wall time on:
  - a locus-sized table (DecodeME harmonized chr1 locus, ~2k rows), cold-ish and warm;
  - one chromosome of the PolyFun prior (chr22, ~200k rows);
  - the whole PolyFun prior (19.5M rows) and whole DecodeME build-37 table (8.9M rows).

Run:
  pixi r python experiments/claude/polyfun_annot_dedup/benchmark_fasta_ref_check.py \
    2>&1 | tee experiments/claude/polyfun_annot_dedup/benchmark_fasta_ref_check.log
"""

import resource
import time
from pathlib import Path

import numpy as np
import polars as pl
from attrs import frozen

STORE = Path("assets/base_asset_store")
HG19 = Path.home() / ".gwaslab/hg19.fa"


@frozen
class FaiEntry:
    length: int
    offset: int
    line_bases: int
    line_bytes: int


def read_fai(fasta: Path) -> dict[str, FaiEntry]:
    out: dict[str, FaiEntry] = {}
    for line in Path(str(fasta) + ".fai").read_text().splitlines():
        name, length, offset, lb, lbytes = line.split("\t")[:5]
        out[name.removeprefix("chr")] = FaiEntry(int(length), int(offset), int(lb), int(lbytes))
    return out


def ref_match(
    mm: np.ndarray, fai: dict[str, FaiEntry], chrom: np.ndarray, pos: np.ndarray, allele: pl.Series
) -> np.ndarray:
    """Boolean array: allele matches the reference starting at pos (1-based), case-insensitive."""
    result = np.zeros(len(pos), dtype=bool)
    lengths = allele.str.len_bytes().to_numpy()
    # One stable sort by (chrom, allele length); each group is then a contiguous slice.
    order = np.lexsort((lengths, chrom))
    sc, sl = chrom[order], lengths[order]
    change = np.flatnonzero((sc[1:] != sc[:-1]) | (sl[1:] != sl[:-1])) + 1
    bounds = np.concatenate(([0], change, [len(order)]))
    for start, stop in zip(bounds[:-1], bounds[1:], strict=True):
        if start == stop:
            continue
        rows = order[start:stop]
        e = fai[str(sc[start])]
        length = int(sl[start])
        if True:
            p0 = pos[rows].astype(np.int64) - 1
            ok = (p0 >= 0) & (p0 + length <= e.length)
            p = p0[:, None] + np.arange(length)[None, :]
            p = np.clip(p, 0, e.length - 1)
            offs = e.offset + p + (p // e.line_bases) * (e.line_bytes - e.line_bases)
            ref = mm[offs] & 0xDF  # ASCII uppercase for letters
            alt = (
                np.frombuffer("".join(allele.gather(rows).to_list()).encode(), dtype=np.uint8)
                .reshape(len(rows), length)
                & 0xDF
            )
            result[rows] = ok & (ref == alt).all(axis=1)
    return result


def bench(label: str, mm: np.ndarray, fai: dict[str, FaiEntry], df: pl.DataFrame, chr_col: str, pos_col: str, ref_col: str) -> None:
    chrom = df[chr_col].cast(pl.Utf8).to_numpy()
    pos = df[pos_col].to_numpy()
    allele = df[ref_col].cast(pl.Utf8)
    t0 = time.perf_counter()
    m = ref_match(mm, fai, chrom, pos, allele)
    dt = time.perf_counter() - t0
    rss_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    print(
        f"{label}: rows={len(pos):,} time={dt:.3f}s ({1e9 * dt / max(len(pos), 1):.0f} ns/row) "
        f"match={m.mean():.4%} peak_rss_so_far={rss_gb:.2f}GB",
        flush=True,
    )


def main() -> None:
    fai = read_fai(HG19)
    mm = np.memmap(HG19, dtype=np.uint8, mode="r")

    locus = pl.read_parquet(
        STORE / "gwas/ME_CFS/DecodeME/processed/decode_me_polyfun_explainchr1_173500000_174500000_palindromes_keep_gwas_harmonized_with_ref.parquet",
        columns=["CHR", "POS", "NEA"],
    )
    bench("locus table, first call", mm, fai, locus, "CHR", "POS", "NEA")
    bench("locus table, second call", mm, fai, locus, "CHR", "POS", "NEA")

    prior_path = STORE / "reference_data/polyfun/precomputed_prior/raw/polyfun_precomputed_heritability_weight_concat.parquet"
    chr22 = pl.scan_parquet(prior_path).filter(pl.col("CHR") == 22).select("CHR", "BP", "A1").collect()
    bench("PolyFun prior chr22", mm, fai, chr22, "CHR", "BP", "A1")

    t0 = time.perf_counter()
    prior = pl.read_parquet(prior_path, columns=["CHR", "BP", "A1"])
    print(f"(read full prior columns: {time.perf_counter() - t0:.1f}s)")
    bench("PolyFun prior genome-wide", mm, fai, prior, "CHR", "BP", "A1")
    del prior

    decode = pl.read_parquet(
        STORE / "gwas/ME_CFS/DecodeME/processed/decode_me_gwas_1_liftover_to_37_parquet_file.parquet",
        columns=["CHR", "POS", "NEA"],
    )
    bench("DecodeME build-37 genome-wide", mm, fai, decode, "CHR", "POS", "NEA")


if __name__ == "__main__":
    main()
