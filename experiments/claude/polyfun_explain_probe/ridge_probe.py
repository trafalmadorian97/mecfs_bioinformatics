"""
THROWAWAY spike. Feasibility check for the polyfun-explainability surrogate:
fit ridge regression of the polyfun precomputed SNPVAR (snpvar_bin) on the
baseline-LF 2.2.UKB annotations, to recover per-annotation weights gamma_c.

This version uses the PRODUCTION-FEASIBLE streaming path so it doubles as a
memory-feasibility demonstration: it never holds the full design matrix, only
one chromosome at a time plus the 187x187 Gram matrix.

  Pass A (train chroms): accumulate per-column mean/std sufficient stats.
  Pass B (train chroms): accumulate G = Xs^T Xs, b = Xs^T y, yty, sum y, n.
  Solve (G + alpha*I) gamma = b.  Train R^2 from sufficient stats; test R^2 by
  streaming the held-out chromosome.

Annotation source: baselineLF2.2.UKB.<chr>.annot.gz (LDSC text, CHR/BP/SNP/CM +
187 annotations). Join to snpvar_meta on rsid (SNP).

Run: pixi r python experiments/claude/polyfun_explain_probe/ridge_probe.py
"""

from pathlib import Path

import numpy as np
import polars as pl
import pyarrow.parquet as pq

HERE = Path(__file__).parent
ANNOT_DIR = HERE / "annot_gz"
META_FILES = [HERE / "snpvar_meta.chr1_7.parquet", HERE / "snpvar_meta.chr8_22.parquet"]
NON_ANNOT = ["CHR", "BP", "SNP", "CM"]


def load_meta() -> pl.DataFrame:
    frames = [pl.from_pandas(pq.read_table(f).to_pandas()) for f in META_FILES if f.exists()]
    # dedup rsid (a few multiallelic dups) so the join is 1:1
    return pl.concat(frames).select(["SNP", "snpvar_bin"]).unique(subset="SNP")


def annot_cols_of(path: Path) -> list[str]:
    header = pl.read_csv(path, separator="\t", n_rows=0).columns
    return [c for c in header if c not in NON_ANNOT]


def load_chr(path: Path, annot_cols: list[str], meta: pl.DataFrame):
    df = pl.read_csv(path, separator="\t", infer_schema_length=None)
    df = df.with_columns([pl.col(c).cast(pl.Float32) for c in annot_cols])
    merged = df.join(meta, on="SNP", how="inner")
    X = merged.select(annot_cols).to_numpy().astype(np.float32)
    y = merged.select("snpvar_bin").to_numpy().ravel().astype(np.float64)
    return X, y, df.height, merged.height


def main() -> None:
    meta = load_meta()
    print(f"meta rows (deduped): {meta.height}")

    files = sorted(ANNOT_DIR.glob("*.annot.gz"),
                   key=lambda f: int(f.name.replace(".annot.gz", "").split(".")[-1]))
    chrom_of = {f: int(f.name.replace(".annot.gz", "").split(".")[-1]) for f in files}
    annot_cols = annot_cols_of(files[0])
    p = len(annot_cols)
    print(f"annotation files: {[f.name for f in files]}  (p={p} annotations)")

    test_file = files[-1]
    train_files = files[:-1]
    print(f"train chroms={[chrom_of[f] for f in train_files]}  test chrom={chrom_of[test_file]}")

    # ---- Pass A: mean/std over train ----
    n = 0
    sx = np.zeros(p)
    sxx = np.zeros(p)
    print("\n=== pass A (mean/std) ===")
    for f in train_files:
        X, y, nraw, njoin = load_chr(f, annot_cols, meta)
        print(f"  {f.name}: annot={nraw} joined={njoin} (rate {njoin/nraw:.3f})")
        Xd = X.astype(np.float64)
        n += Xd.shape[0]
        sx += Xd.sum(0)
        sxx += (Xd * Xd).sum(0)
        del X, Xd
    mu = sx / n
    var = sxx / n - mu**2
    sd = np.sqrt(np.maximum(var, 0))
    sd[sd == 0] = 1.0

    # ---- Pass B: Gram accumulation over train ----
    G = np.zeros((p, p))
    b = np.zeros(p)
    yty = 0.0
    sy = 0.0
    ny = 0
    print("=== pass B (Gram) ===")
    for f in train_files:
        X, y, _, _ = load_chr(f, annot_cols, meta)
        Xs = (X.astype(np.float64) - mu) / sd
        G += Xs.T @ Xs
        b += Xs.T @ y
        yty += float(y @ y)
        sy += float(y.sum())
        ny += y.shape[0]
        del X, Xs, y
    ss_tot_tr = yty - sy * sy / ny
    print(f"train n={ny}  Gram={G.shape} ({G.nbytes/1e6:.1f} MB, independent of n)")

    # test chromosome (small enough to hold)
    Xte, yte, _, _ = load_chr(test_file, annot_cols, meta)
    Xte_s = (Xte.astype(np.float64) - mu) / sd
    ss_tot_te = float(((yte - yte.mean()) ** 2).sum())

    print("\n=== ridge: cross-chromosome R^2 ===")
    I = np.eye(p)
    best = None
    for alpha in [0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 1e5]:
        gamma = np.linalg.solve(G + alpha * I, b)
        # train R^2 from sufficient stats: SS_res = yty - 2 g.b + g.G.g
        ss_res_tr = yty - 2 * gamma @ b + gamma @ G @ gamma
        r2_tr = 1 - ss_res_tr / ss_tot_tr
        pred_te = Xte_s @ gamma
        r2_te = 1 - float(((yte - pred_te) ** 2).sum()) / ss_tot_te
        print(f"alpha={alpha:9.1f}  train R^2={r2_tr:.4f}  TEST R^2={r2_te:.4f}")
        if best is None or r2_te > best[1]:
            best = (alpha, r2_te, gamma)

    alpha, r2_te, gamma = best
    print(f"\nbest alpha={alpha}  test R^2={r2_te:.4f}")
    order = np.argsort(-np.abs(gamma))
    print("\n=== top 25 annotations by |standardized gamma_c| ===")
    for i in order[:25]:
        print(f"  {gamma[i]:+.4e}  {annot_cols[i]}")


if __name__ == "__main__":
    main()
