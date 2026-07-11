"""
Experiment: sensitivity of DecodeME PoPS results to --feature_selection_max_num.

POPs selects every feature passing marginal p<0.05 by default (no cap). The paper
reports 2,512 to 26,155 features selected per trait, so an uncapped run can hold a
very wide feature matrix in memory during training (plus the RidgeCV SVD), which is
the remaining memory risk after the munge OOM is fixed by the streaming splitter.

This script runs pops.py directly against the cached munged features and the
DecodeME MAGMA ensembl gene analysis for several caps, recording peak RSS and the
number of features actually selected, then compares the resulting PoPS gene rankings
across caps (Spearman correlation + top-K overlap) to judge whether capping
materially changes the model.

Run one cap at a time (memory-bound, so keep runs sequential):
    pixi r python experiments/claude/pops_feature_cap_experiment.py run <cap|none>
Then analyze once the runs are done:
    pixi r python experiments/claude/pops_feature_cap_experiment.py analyze
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
# Invoke the pixi default-env interpreter directly (not via "pixi r") so
# /usr/bin/time measures the pops.py process itself rather than the pixi wrapper.
ENV_PYTHON = REPO / ".pixi/envs/default/bin/python"
ASSET_STORE = REPO / "assets" / "base_asset_store"
POPS_SOURCE = (
    ASSET_STORE / "reference_data/pops/source/extracted/pops_source_extracted"
)
MUNGED_DIR = ASSET_STORE / "other_files/pops_features_munged"
MUNGED_PREFIX = MUNGED_DIR / "pops_features"
MAGMA_PREFIX = (
    ASSET_STORE
    / "gwas/ME_CFS/DecodeME/processed/magma"
    / "decode_me_gwas_1_build_37_magma_pops_gene_analysis"
    / "gene_analysis_output"
)
GENE_ANNOT = POPS_SOURCE / "example/data/utils/gene_annot_jun10.txt"
CONTROL_FEATURES = POPS_SOURCE / "example/data/utils/features_jul17_control.txt"
OUT_ROOT = REPO / "experiments/claude/pops_feature_cap"

# "none" = uncapped (POPs default). For DecodeME the uncapped run selects 5,383
# features, so caps at or above that are non-binding (identical to uncapped); the
# informative caps are the binding ones below 5,383.
CAPS = ["none", "5000", "2500", "1000"]


def _num_feature_chunks() -> int:
    chunks = list(MUNGED_DIR.glob("pops_features.cols.*.txt"))
    assert chunks, f"No munged feature chunks found in {MUNGED_DIR}"
    return len(chunks)


def _out_prefix(cap: str) -> Path:
    return OUT_ROOT / f"cap_{cap}" / "pops_results"


def run_one(cap: str) -> None:
    out_prefix = _out_prefix(cap)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    args = [
        "/usr/bin/time",
        "-v",
        str(ENV_PYTHON),
        str(POPS_SOURCE / "pops.py"),
        "--gene_annot_path",
        str(GENE_ANNOT),
        "--feature_mat_prefix",
        str(MUNGED_PREFIX),
        "--num_feature_chunks",
        str(_num_feature_chunks()),
        "--magma_prefix",
        str(MAGMA_PREFIX),
        "--control_features_path",
        str(CONTROL_FEATURES),
        "--out_prefix",
        str(out_prefix),
        "--verbose",
    ]
    if cap != "none":
        args += ["--feature_selection_max_num", cap]
    log_path = out_prefix.parent / "run.log"
    print(f"[cap={cap}] running -> {log_path}")
    with log_path.open("w", encoding="utf-8") as log:
        proc = subprocess.run(args, stdout=log, stderr=subprocess.STDOUT)
    print(f"[cap={cap}] exit={proc.returncode}")


def _load_scores(cap: str) -> pd.Series:
    preds = pd.read_csv(str(_out_prefix(cap)) + ".preds", sep="\t")
    return preds.set_index("ENSGID")["PoPS_Score"]


def _num_selected(cap: str) -> int:
    marg = pd.read_csv(str(_out_prefix(cap)) + ".marginals", sep="\t", index_col=0)
    return int(marg["selected"].sum())


def analyze() -> None:
    available = [c for c in CAPS if Path(str(_out_prefix(c)) + ".preds").exists()]
    assert available, "No completed runs to analyze."
    print(f"Completed runs: {available}\n")
    for cap in available:
        print(f"  cap={cap:>6}  features_selected={_num_selected(cap)}")
    print()
    if "none" not in available:
        print("Uncapped run missing; skipping cross-cap comparison.")
        return
    base = _load_scores("none")
    top20_base = set(base.sort_values(ascending=False).head(20).index)
    print("Comparison vs uncapped (none):")
    print(f"  {'cap':>6}  {'spearman':>9}  {'top20_overlap':>13}")
    for cap in available:
        if cap == "none":
            continue
        scores = _load_scores(cap)
        joined = pd.concat([base, scores], axis=1, join="inner")
        spearman = joined.iloc[:, 0].corr(joined.iloc[:, 1], method="spearman")
        top20_cap = set(scores.sort_values(ascending=False).head(20).index)
        overlap = len(top20_base & top20_cap)
        print(f"  {cap:>6}  {spearman:>9.4f}  {overlap:>11}/20")


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "run":
        run_one(sys.argv[2])
    elif len(sys.argv) >= 2 and sys.argv[1] == "analyze":
        analyze()
    else:
        print(__doc__)
        sys.exit(1)
