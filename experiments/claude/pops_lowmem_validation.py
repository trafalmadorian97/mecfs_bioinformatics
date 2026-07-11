"""
Validate the low-peak-memory POPs reimplementation against stock POPs on real data.

The chr22 example subset already showed machine-precision agreement across the full
pipeline; this harness confirms the same on the real DecodeME and IBD runs and
demonstrates the payoff: IBD selects ~14,564 features uncapped and OOM-kills the box
under stock POPs, but the low-memory path runs it uncapped.

Runs (each writes run.log with /usr/bin/time -v peak RSS):
  run_decode        low-mem on DecodeME (compare to stock uncapped in pops_feature_cap/cap_none)
  run_ibd_capped    low-mem on IBD capped 9000 (compare to stock in pops_ibd)
  run_ibd_uncapped  low-mem on IBD uncapped (stock OOMs here; the payoff run)
  analyze           print fidelity + peak-RSS table

    pixi r python experiments/claude/pops_lowmem_validation.py run_decode
    pixi r python experiments/claude/pops_lowmem_validation.py analyze
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ENV_PYTHON = REPO / ".pixi/envs/default/bin/python"
ASSET_STORE = REPO / "assets" / "base_asset_store"
SOURCE = ASSET_STORE / "reference_data/pops/source/extracted/pops_source_extracted"
LOWMEM_SCRIPT = REPO / "mecfs_bio/build_system/task/pops/lowmem/pops_lowmem.py"
MUNGED_PREFIX = ASSET_STORE / "other_files/pops_features_munged/pops_features"
GENE_ANNOT = SOURCE / "example/data/utils/gene_annot_jun10.txt"
CONTROL = SOURCE / "example/data/utils/features_jul17_control.txt"

DECODE_MAGMA = (
    ASSET_STORE / "gwas/ME_CFS/DecodeME/processed/magma"
    "/decode_me_gwas_1_build_37_magma_pops_gene_analysis/gene_analysis_output"
)
IBD_MAGMA = (
    ASSET_STORE / "gwas/inflammatory_bowl_disease/liu_et_al_2023_ibd_meta/processed/magma"
    "/liu_et_al_ibd_2023_eur_37_magma_pops_gene_analysis/gene_analysis_output"
)

OUT = REPO / "experiments/claude/pops_lowmem"
# Stock reference outputs produced earlier.
STOCK_DECODE = REPO / "experiments/claude/pops_feature_cap/cap_none/pops_results"
STOCK_IBD_CAPPED = REPO / "experiments/claude/pops_ibd/pops_results"

RUNS = {
    "run_decode": ("decode", DECODE_MAGMA, None),
    "run_ibd_capped": ("ibd_capped", IBD_MAGMA, "9000"),
    "run_ibd_uncapped": ("ibd_uncapped", IBD_MAGMA, None),
}


def _num_chunks() -> int:
    return len(list(MUNGED_PREFIX.parent.glob("pops_features.cols.*.txt")))


def run(name: str) -> None:
    label, magma, cap = RUNS[name]
    out_dir = OUT / label
    out_dir.mkdir(parents=True, exist_ok=True)
    args = [
        "/usr/bin/time", "-v", str(ENV_PYTHON), str(LOWMEM_SCRIPT),
        "--pops_source_dir", str(SOURCE),
        "--gene_annot_path", str(GENE_ANNOT),
        "--feature_mat_prefix", str(MUNGED_PREFIX),
        "--num_feature_chunks", str(_num_chunks()),
        "--magma_prefix", str(magma),
        "--control_features_path", str(CONTROL),
        "--out_prefix", str(out_dir / "pops_results"),
        "--verbose",
    ]
    if cap is not None:
        args += ["--feature_selection_max_num", cap]
    log_path = out_dir / "run.log"
    print(f"running low-mem {label} -> {log_path}")
    with log_path.open("w", encoding="utf-8") as log:
        proc = subprocess.run(args, stdout=log, stderr=subprocess.STDOUT)
    print(f"exit={proc.returncode}")


def _preds(prefix: Path) -> pd.Series:
    return pd.read_csv(str(prefix) + ".preds", sep="\t").set_index("ENSGID")["PoPS_Score"]


def _coefs(prefix: Path) -> pd.Series:
    return pd.read_csv(str(prefix) + ".coefs", sep="\t").set_index("parameter")["beta"]


def _peak_rss_gb(log_path: Path) -> float | None:
    if not log_path.exists():
        return None
    for line in log_path.read_text().splitlines():
        if "Maximum resident set size" in line:
            return int(line.split(":")[-1].strip()) / 1e6
    return None


def _compare(label: str, stock_prefix: Path, low_prefix: Path) -> None:
    if not Path(str(low_prefix) + ".preds").exists():
        print(f"[{label}] low-mem .preds missing (run not finished?)")
        return
    s, low = _preds(stock_prefix), _preds(low_prefix)
    j = pd.concat([s, low], axis=1, join="inner")
    j.columns = ["s", "l"]
    sc, lc = _coefs(stock_prefix), _coefs(low_prefix)
    meta = ["METHOD", "SELECTED_CV_ALPHA", "BEST_CV_SCORE"]
    sb = sc[~sc.index.isin(meta)].astype(float)
    lb = lc[~lc.index.isin(meta)].astype(float)
    jb = pd.concat([sb, lb], axis=1, join="inner")
    jb.columns = ["s", "l"]
    print(f"[{label}]")
    print(f"  preds: n={len(j)} max_abs_diff={np.abs(j.s - j.l).max():.2e} "
          f"spearman={j.s.corr(j.l, method='spearman'):.6f}")
    print(f"  alpha: stock={sc.get('SELECTED_CV_ALPHA')} lowmem={lc.get('SELECTED_CV_ALPHA')}")
    print(f"  coefs: n={len(jb)} same_set={set(sb.index) == set(lb.index)} "
          f"max_abs_diff={np.abs(jb.s - jb.l).max():.2e}")
    print(f"  peak RSS: stock={_peak_rss_gb(stock_prefix.parent / 'run.log')} GB "
          f"lowmem={_peak_rss_gb(low_prefix.parent / 'run.log')} GB")


def analyze() -> None:
    _compare("DecodeME uncapped", STOCK_DECODE, OUT / "decode" / "pops_results")
    _compare("IBD capped 9000", STOCK_IBD_CAPPED, OUT / "ibd_capped" / "pops_results")
    # Uncapped IBD: stock OOMs, so compare low-mem uncapped to the capped runs instead.
    up = OUT / "ibd_uncapped" / "pops_results"
    if Path(str(up) + ".preds").exists():
        annot = pd.read_csv(GENE_ANNOT, sep="\t").set_index("ENSGID")["NAME"]
        u = _preds(up)
        print("[IBD uncapped low-mem] (stock OOMs at this feature count)")
        print(f"  n_genes={len(u)} peak RSS={_peak_rss_gb(up.parent / 'run.log')} GB")
        cap = _preds(STOCK_IBD_CAPPED)
        jj = pd.concat([u, cap], axis=1, join="inner").dropna()
        jj.columns = ["uncapped", "capped9000"]
        print(f"  vs stock-capped-9000: spearman={jj.uncapped.corr(jj.capped9000, method='spearman'):.4f}")
        top = u.sort_values(ascending=False).head(15)
        print("  top 15 genes: " + ", ".join(annot.get(g, "?") for g in top.index))
    else:
        print("[IBD uncapped low-mem] not run yet")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) >= 2 else ""
    if cmd in RUNS:
        run(cmd)
    elif cmd == "analyze":
        analyze()
    else:
        print(__doc__)
        sys.exit(1)
