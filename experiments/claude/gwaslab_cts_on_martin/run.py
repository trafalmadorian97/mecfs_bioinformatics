"""
Experiment: run gwaslab's S-LDSC cell-type regression on MARTIN's munged sumstats.

GOAL -- isolate whether the Peter-vs-Martin cell-type p-value gap is due to the input SNP set or
the S-LDSC implementation (gwaslab's port vs original LDSC 1.0.1).

We feed Martin's munged `.sumstats.gz` (columns SNP, A1, A2, Z, N -- his exact ~1,157,905-SNP set)
directly into gwaslab's cts regression, BYPASSING gwaslab's internal `_get_hapmap3` re-restriction
(which would drop his ~5k unique SNPs and collapse the set back toward Peter's). This is the
faithful analog of what `ldsc.py` does: merge the munged sumstats with the LD-score files and
regress. We call the lower-level `_estimate_h2_cts_by_ldsc` directly to skip `_get_hapmap3`.

INTERPRETATION (compare the headline tissue p-value):
  gwaslab(Martin input) ~= Peter's gwaslab value  -> the input SNP set does NOT explain the gap
                                                     -> the gap is the implementation (gwaslab vs LDSC).
  gwaslab(Martin input) shifts toward "less significant" -> the input SNP set matters.

We run with baseline v1.2 (Peter's) for direct comparison to Peter, and v1.0 (Martin's likely
baseline) for a Martin-replica. Peter reference values (fresh gwaslab 4.1.7 from CI run
26691178699): multitissue A08.186.211.Brain = 2.334e-7, Brain_Cortex = 2.347e-6;
gtex_brain Brain_Cortex ~= 0.000106.
"""

import sys
import shutil
import tempfile
from pathlib import Path

import pandas as pd
from gwaslab.info.g_Log import Log
from gwaslab.util.util_ex_ldsc import _estimate_h2_cts_by_ldsc

from mecfs_bio.build_system.task.gwaslab.gwaslab_cell_analysis_by_sldsc import (
    prepend_path_to_ldcts_file,
)

REPO = Path(__file__).resolve().parents[3]
MARTIN_SUMSTATS = REPO / "gwas_1.regenie.filtered.rsids.munged.txt.sumstats.gz"
# Peter's HapMap3-filtered sumstats (has CHR/POS per rsID); used only to supply CHR/POS columns
# that gwaslab's cts column-check decorator requires (they are not used in the regression itself).
PETER_CSV = REPO / "martin_exchange_5/hapmap_filtered_sumstats.csv"
OUT_DIR = Path(__file__).resolve().parent / "results"

_REF = REPO / "assets/base_asset_store/reference_data/linkage_disequilibrium_scores"
_AUX = (
    REPO
    / "assets/base_asset_store/reference_data/linkage_disequilibrium_score_regression_aux_data"
)

# Un-partitioned regression weights (hm3 no-MHC) -- shared across runs.
W_LD_CHR = str(
    _AUX
    / "LDSCORE_1000G_Phase/extracted/s_ldsc_hapmap3_regression_weights_extracted"
    / "1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.@"
)

# 'baseline' partitioned-model references (the conditioning model).
BASELINES = {
    "v1.2": str(
        _REF
        / "LDSCORE_1000G_Phase3_baseline/extracted/base_model_partitioned_ld_scores_extracted"
        / "baseline_v1.2/baseline.@"
    ),
    "v1.0": str(
        _REF
        / "LDSCORE_1000G_Phase3_baseline_v1_0/extracted/base_model_partitioned_ld_scores_v1_0_extracted"
        / "1000G_EUR_Phase3_baseline/baseline.@"
    ),
}

# Cell-type-specific LD-score datasets: (extracted dir, .ldcts index filename, headline tissue).
CTS_DATASETS = {
    "gtex_brain": (
        _REF
        / "LDSCORE_LDSC_SEG/extracted/partitioned_model_gtex_brain_ld_scores_extracted",
        "GTEx_brain.ldcts",
        "Brain_Cortex",
    ),
    "multitissue": (
        _REF
        / "LDSCORE_LDSC_SEG/extracted/multi_tissue_gene_expression_partitioned_ld_scores_extracted",
        "Multi_tissue_gene_expr.ldcts",
        "A08.186.211.Brain",
    ),
}


def load_martin_sumstats() -> pd.DataFrame:
    """
    Martin's munged file: SNP A1 A2 Z N. Drop the merge-alleles padding rows (NaN Z).

    gwaslab's cts entry point is decorated with a presence-check for CHR/POS/EA/NEA (none of which
    the cts regression actually uses -- it merges on SNP and regresses chi^2 = Z^2). We add EA/NEA
    from A1/A2 and bring CHR/POS in from Peter's HapMap3 CSV so the check passes. Martin-unique SNPs
    (absent from Peter) simply get NaN CHR/POS, which is fine (column presence is all that matters).
    """
    df = pd.read_csv(MARTIN_SUMSTATS, sep=r"\s+").dropna(subset=["Z"])
    df["EA"] = df["A1"]
    df["NEA"] = df["A2"]
    chr_pos = pd.read_csv(PETER_CSV)[["rsID", "CHR", "POS"]]
    df = df.merge(chr_pos, left_on="SNP", right_on="rsID", how="left")
    return df


def load_peter_as_munged() -> pd.DataFrame:
    """Peter's HapMap3 CSV mapped to the same SNP/A1/A2/Z/N schema (diagnostic control)."""
    df = pd.read_csv(PETER_CSV)
    df = df.rename(columns={"rsID": "SNP", "EA": "A1", "NEA": "A2"})
    df["Z"] = df["BETA"] / df["SE"]
    df["EA"] = df["A1"]
    df["NEA"] = df["A2"]
    return df[["SNP", "A1", "A2", "Z", "N", "EA", "NEA", "CHR", "POS"]]


def run_one(
    baseline_label: str,
    dataset_label: str,
    df: pd.DataFrame | None = None,
    tag: str = "martin_input",
) -> pd.DataFrame:
    cts_dir, index_name, headline = CTS_DATASETS[dataset_label]
    if df is None:
        df = load_martin_sumstats()
    with tempfile.TemporaryDirectory() as tmp:
        # The .ldcts index stores paths relative to its own directory; copy the cts dir and
        # prepend an absolute prefix (same workaround as CellAnalysisByLDSCTask).
        local_cts = Path(tmp) / "cts"
        shutil.copytree(cts_dir, local_cts)
        index_path = local_cts / index_name
        prepend_path_to_ldcts_file(
            index_path, prefix_path=local_cts, output_path=index_path
        )
        summary = _estimate_h2_cts_by_ldsc(
            df,
            log=Log(),
            verbose=True,
            ref_ld_chr=BASELINES[baseline_label],
            ref_ld_chr_cts=str(index_path),
            w_ld_chr=W_LD_CHR,
        )
    summary = summary.sort_values("Coefficient_P_value").reset_index(drop=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / f"{tag}_{dataset_label}_baseline_{baseline_label}.csv"
    summary.to_csv(out_path, index=False)
    head_p = summary.loc[summary["Name"] == headline, "Coefficient_P_value"]
    print(
        f"\n==> [{tag} | {dataset_label} / baseline {baseline_label}] n_snps_input={len(df)}  "
        f"{headline} p = {head_p.iloc[0]:.4e}  (saved {out_path.name})"
    )
    print(summary.head(6).to_string(index=False))
    return summary


def _peter_full_annovar_rsids() -> set[str]:
    """rsIDs in Peter's full (pre-HapMap3) annovar-assigned sumstats (~8.8M)."""
    full = pd.read_parquet(
        REPO
        / "assets/base_asset_store/gwas/ME_CFS/DecodeME/processed/decode_me_gwas_1_assign_rsids_via_dbsnp150.parquet",
        columns=["rsid"],
    )
    return set(full["rsid"])


def run_followup_2(dataset_label: str = "gtex_brain", baseline_label: str = "v1.2"):
    """
    Split Martin's ~5,318 unique SNPs (those collapsing the cts result) by their upstream cause:

      A = present in Peter's FULL annovar data but dropped by gwaslab's HapMap3 snplist
          -> isolates the gwaslab-HapMap3-snplist choice (db150 subset vs canonical w_hm3)
      B = absent from Peter's data entirely
          -> isolates upstream rsID-assignment disagreement

    Run cts on (common + A) and (common + B) to see which group drives the collapse.
    Reference points: common-only -> ~1.07e-4 (matches Peter); common + A + B (full) -> 0.33.
    """
    m = load_martin_sumstats()
    peter_hm3 = set(pd.read_csv(PETER_CSV)["rsID"])
    common = m[m["SNP"].isin(peter_hm3)]
    unique = m[~m["SNP"].isin(peter_hm3)]
    in_peter_full = unique["SNP"].isin(_peter_full_annovar_rsids())
    group_a = unique[in_peter_full]  # real SNPs gwaslab's HapMap3 list dropped
    group_b = unique[~in_peter_full]  # upstream rsID-assignment artifacts
    print(
        f"common={len(common)}  A(gwaslab-HM3-list dropped, real)={len(group_a)}  "
        f"B(rsID-assignment artifact)={len(group_b)}"
    )
    run_one(
        baseline_label,
        dataset_label,
        df=pd.concat([common, group_a]),
        tag="common_plus_A_gwaslab_hm3_list",
    )
    run_one(
        baseline_label,
        dataset_label,
        df=pd.concat([common, group_b]),
        tag="common_plus_B_rsid_assignment",
    )


def main():

    if sys.argv[1:2] == ["followup2"]:
        run_followup_2(*sys.argv[2:])
        return

    # default: validate on the small/fast gtex_brain first, then the headline multitissue.
    jobs = sys.argv[1:] or [
        "gtex_brain:v1.2",
        "gtex_brain:v1.0",
        "multitissue:v1.2",
        "multitissue:v1.0",
    ]
    for job in jobs:
        dataset_label, baseline_label = job.split(":")
        run_one(baseline_label=baseline_label, dataset_label=dataset_label)


if __name__ == "__main__":
    main()
