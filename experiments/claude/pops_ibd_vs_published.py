"""
Sanity check: compare our IBD PoPS scores against the published POPs results.

Runs pops.py on the Liu et al 2023 IBD POPs-specific MAGMA (uncapped, default
settings, matching PopsRunTask) and compares the resulting per-gene PoPS scores to
the manuscript's IBD PoPS (PoPS_FullResults.txt.gz, trait == IBD, cohorts PASS and
UKB) via Spearman correlation and top-K gene overlap.

The underlying GWAS differ (ours is Liu 2023; the paper's IBD is PASS/UKB), so we do
not expect identical scores; a strong positive correlation and high top-gene overlap
indicate our pipeline reproduces the trait biology.

    pixi r python experiments/claude/pops_ibd_vs_published.py run
    pixi r python experiments/claude/pops_ibd_vs_published.py analyze
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ENV_PYTHON = REPO / ".pixi/envs/default/bin/python"
ASSET_STORE = REPO / "assets" / "base_asset_store"
POPS_SOURCE = ASSET_STORE / "reference_data/pops/source/extracted/pops_source_extracted"
MUNGED_DIR = ASSET_STORE / "other_files/pops_features_munged"
MUNGED_PREFIX = MUNGED_DIR / "pops_features"
MAGMA_PREFIX = (
    ASSET_STORE
    / "gwas/inflammatory_bowl_disease/liu_et_al_2023_ibd_meta/processed/magma"
    / "liu_et_al_ibd_2023_eur_37_magma_pops_gene_analysis"
    / "gene_analysis_output"
)
GENE_ANNOT = POPS_SOURCE / "example/data/utils/gene_annot_jun10.txt"
CONTROL_FEATURES = POPS_SOURCE / "example/data/utils/features_jul17_control.txt"
OUT_DIR = REPO / "experiments/claude/pops_ibd"
OUT_PREFIX = OUT_DIR / "pops_results"
PUBLISHED = REPO / "PoPS_FullResults.txt.gz"


def _num_feature_chunks() -> int:
    chunks = list(MUNGED_DIR.glob("pops_features.cols.*.txt"))
    assert chunks, f"No munged feature chunks found in {MUNGED_DIR}"
    return len(chunks)


def run(feature_selection_max_num: str | None = None) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
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
        str(OUT_PREFIX),
        "--verbose",
    ]
    # IBD selects ~14.5k features uncapped, which OOMs the 15GB box during training.
    # A cap keeps peak memory in bounds for the comparison; see run.log for the
    # actual peak RSS. Capping mildly perturbs scores (see the cap-sensitivity
    # experiment) but preserves the strongest genes, which drive the comparison.
    if feature_selection_max_num is not None:
        args += ["--feature_selection_max_num", feature_selection_max_num]
    log_path = OUT_DIR / "run.log"
    print(f"running IBD POPs -> {log_path}")
    with log_path.open("w", encoding="utf-8") as log:
        proc = subprocess.run(args, stdout=log, stderr=subprocess.STDOUT)
    print(f"exit={proc.returncode}")


def _our_scores() -> pd.Series:
    preds = pd.read_csv(str(OUT_PREFIX) + ".preds", sep="\t")
    return preds.set_index("ENSGID")["PoPS_Score"]


def _published_scores(cohort: str) -> pd.Series:
    df = pd.read_csv(PUBLISHED, sep="\t", compression="gzip")
    df = df[(df["trait"] == "IBD") & (df["cohort"] == cohort)]
    return df.set_index("ensgid")["pops_score"]


def analyze() -> None:
    ours = _our_scores()
    annot = pd.read_csv(GENE_ANNOT, sep="\t").set_index("ENSGID")["NAME"]
    print(f"Our IBD PoPS: {len(ours)} genes\n")
    print("Top 15 genes (ours):")
    top = ours.sort_values(ascending=False).head(15)
    for ensgid, score in top.items():
        print(f"  {annot.get(ensgid, '?'):<12} {ensgid}  {score:.3f}")
    print()
    for cohort in ["PASS", "UKB"]:
        pub = _published_scores(cohort)
        joined = pd.concat([ours, pub], axis=1, join="inner").dropna()
        joined.columns = ["ours", "pub"]
        spearman = joined["ours"].corr(joined["pub"], method="spearman")
        pearson = joined["ours"].corr(joined["pub"], method="pearson")
        for k in (20, 50, 100):
            our_top = set(ours.sort_values(ascending=False).head(k).index)
            pub_top = set(pub.sort_values(ascending=False).head(k).index)
            overlap = len(our_top & pub_top)
            print(
                f"vs IBD/{cohort}: n={len(joined)} spearman={spearman:.3f} "
                f"pearson={pearson:.3f} top{k}_overlap={overlap}/{k}"
            )
        print()


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "run":
        run(sys.argv[2] if len(sys.argv) >= 3 else None)
    elif len(sys.argv) >= 2 and sys.argv[1] == "analyze":
        analyze()
    else:
        print(__doc__)
        sys.exit(1)
