"""
Check that Western et al. 2024 CSF aptamer GWAS are on a unit-SD scale (sd_y = 1).

For a quantitative trait, Var(beta_hat) ~= sdY^2 / (2 N maf (1 - maf)), so sdY is estimated as in
coloc's sdY.est: regress SE^2 on 1 / (2 N maf (1 - maf)) through the origin and take the square root
of the slope. Done on three 1 Mb windows (remote tabix reads of the GWAS Catalog harmonised files) for a small
random sample of aptamers, MAF >= 0.01, plus genome-wide on the locally cached template aptamer.

Usage: pixi r python experiments/claude/csf_coloc/sdy_check.py <manifest.csv> <template.parquet> [n_aptamers] | tee log
"""

import sys

import numpy as np
import polars as pl
import pysam

FTP = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
_WINDOWS = [("2", 100_000_000, 101_000_000), ("11", 60_000_000, 61_000_000), ("22", 30_000_000, 31_000_000)]


def bucket(acc: str) -> str:
    n = int(acc[4:])
    lo = ((n - 1) // 1000) * 1000 + 1
    w = len(str(n))
    return f"GCST{lo:0{w}d}-GCST{lo + 999:0{w}d}"


def sdy_est(acc: str) -> tuple[float, int]:
    url = f"{FTP}/{bucket(acc)}/{acc}/harmonised/{acc}.h.tsv.gz"
    se, eaf, n = [], [], []
    tabix = pysam.TabixFile(url)
    rows = [r for c, a, b in _WINDOWS for r in tabix.fetch(c, a, b)]
    for row in rows:
        f = row.split("\t")
        try:
            se.append(float(f[5]))
            eaf.append(float(f[6]))
            n.append(float(f[12]))
        except ValueError:
            continue
    return _estimate(np.array(se), np.array(eaf), np.array(n))


def _estimate(se_a: np.ndarray, eaf_a: np.ndarray, n_a: np.ndarray) -> tuple[float, int]:
    maf = np.minimum(eaf_a, 1 - eaf_a)
    keep = maf >= 0.01
    x = 1.0 / (2 * n_a[keep] * maf[keep] * (1 - maf[keep]))
    y = se_a[keep] ** 2
    slope = float((x @ y) / (x @ x))
    return float(np.sqrt(slope)), int(keep.sum())


def main() -> None:
    manifest = pl.read_csv(sys.argv[1])
    t = pl.read_parquet(sys.argv[2], columns=["standard_error", "effect_allele_frequency", "n"])
    est, m = _estimate(t["standard_error"].to_numpy(), t["effect_allele_frequency"].to_numpy(), t["n"].to_numpy().astype(float))
    print(f"template aptamer X13681.173 (CSNK2A2) genome-wide  sdY_est={est:.4f}  variants={m}", flush=True)
    k = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    sample = manifest.sample(n=k, seed=20261106)
    rows = []
    for r in sample.iter_rows(named=True):
        est, m = sdy_est(r["accession"])
        print(f"{r['analyte']:>12} {r['entrez_gene_symbol']:>10} {r['accession']}  sdY_est={est:.4f}  variants={m}", flush=True)
        rows.append(est)
    a = np.array(rows)
    print(f"\nn={len(a)}  mean={a.mean():.4f}  min={a.min():.4f}  max={a.max():.4f}")


if __name__ == "__main__":
    main()
