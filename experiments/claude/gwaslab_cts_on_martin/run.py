"""
Experiment: run gwaslab's S-LDSC cell-type regression on controlled variations of Martin's AND
Peter's summary statistics, through an IDENTICAL harness, to localize the Peter-vs-Martin gap.

Each run is identified by three names: (dataframe_name, dataset_label, baseline_label). Together
they uniquely determine both the input that is run and the output path
`results/{dataframe_name}_{dataset_label}_baseline_{baseline_label}.csv`, so nothing is ever
silently overwritten and the experiment is fully controlled.

Method: feed a munged sumstats DataFrame (columns SNP, A1, A2, Z, N, ...) straight into gwaslab's
cts regression, BYPASSING gwaslab's internal `_get_hapmap3` re-restriction (call the lower-level
`_estimate_h2_cts_by_ldsc` directly). This is the faithful analog of what `ldsc.py` does: merge the
munged sumstats with the LD-score files and regress chi^2 = Z^2.

KEY RESULT (gtex_brain, baseline v1.2, Brain_Cortex headline):
  peter                                    -> 1.06e-4   (== Peter's documented value; harness control)
  martin (full)                            -> 0.33      (collapse)
  martin_restricted_to_peter_snps          -> 1.07e-4   (== Peter; given the same SNPs, gwaslab agrees)
  common_snps_plus_present_in_peter_annovar -> 1.3e-4   (benign: Martin-unique SNPs Peter HAS but
                                                          gwaslab's HapMap3 list dropped)
  common_snps_plus_absent_from_peter_annovar -> 0.33    (the collapse: Martin-unique SNPs absent from
                                                          Peter's annovar output)
So the driver is the SNP set -- specifically the ~2,600 SNPs absent from Peter's annovar output, NOT
the implementation, baseline, or gwaslab HapMap3 list.
"""

import shutil
import sys
import tempfile
from functools import lru_cache
from pathlib import Path

import pandas as pd
from gwaslab.info.g_Log import Log
from gwaslab.util.util_ex_ldsc import _estimate_h2_cts_by_ldsc

from mecfs_bio.build_system.task.gwaslab.gwaslab_cell_analysis_by_sldsc import (
    prepend_path_to_ldcts_file,
)

REPO = Path(__file__).resolve().parents[3]
MARTIN_SUMSTATS = REPO / "gwas_1.regenie.filtered.rsids.munged.txt.sumstats.gz"
# Peter's HapMap3-filtered sumstats (his exact pre-S-LDSC SNP set, with CHR/POS/EA/NEA/BETA/SE).
PETER_CSV = REPO / "martin_exchange_5/hapmap_filtered_sumstats.csv"
# Peter's full pre-HapMap3 annovar-assigned sumstats (~8.8M), used to classify Martin-unique SNPs.
PETER_FULL_ANNOVAR = (
    REPO
    / "assets/base_asset_store/gwas/ME_CFS/DecodeME/processed/decode_me_gwas_1_assign_rsids_via_dbsnp150.parquet"
)
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


# --------------------------------------------------------------------------------------------------
# Base loaders (cached, so the full experiment reuses them across runs)
# --------------------------------------------------------------------------------------------------
@lru_cache(maxsize=1)
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


@lru_cache(maxsize=1)
def load_peter_as_munged() -> pd.DataFrame:
    """Peter's HapMap3 CSV mapped to the same SNP/A1/A2/Z/N schema (the harness control)."""
    df = pd.read_csv(PETER_CSV)
    df = df.rename(columns={"rsID": "SNP", "EA": "A1", "NEA": "A2"})
    df["Z"] = df["BETA"] / df["SE"]
    df["EA"] = df["A1"]
    df["NEA"] = df["A2"]
    return df[["SNP", "A1", "A2", "Z", "N", "EA", "NEA", "CHR", "POS"]]


@lru_cache(maxsize=1)
def _peter_hm3_rsids() -> frozenset[str]:
    return frozenset(pd.read_csv(PETER_CSV)["rsID"])


@lru_cache(maxsize=1)
def _peter_full_annovar_rsids() -> frozenset[str]:
    return frozenset(pd.read_parquet(PETER_FULL_ANNOVAR, columns=["rsid"])["rsid"])


# --------------------------------------------------------------------------------------------------
# Dataframe registry: name -> fetcher. Each fetcher returns one controlled input variation.
# --------------------------------------------------------------------------------------------------
def _martin_common() -> pd.DataFrame:
    m = load_martin_sumstats()
    return m[m["SNP"].isin(_peter_hm3_rsids())]


def _martin_unique() -> pd.DataFrame:
    m = load_martin_sumstats()
    return m[~m["SNP"].isin(_peter_hm3_rsids())]


def _fetch_martin() -> pd.DataFrame:
    return load_martin_sumstats()


def _fetch_peter() -> pd.DataFrame:
    return load_peter_as_munged()


def _fetch_martin_restricted_to_peter_snps() -> pd.DataFrame:
    return _martin_common()


def _fetch_common_snps_plus_present_in_peter_annovar() -> pd.DataFrame:
    # common + Martin-unique SNPs that ARE present in Peter's full annovar data (they're absent only
    # from the HapMap3-filtered comparison set because gwaslab's HapMap3 list dropped them). Benign.
    uniq = _martin_unique()
    group_a = uniq[uniq["SNP"].isin(_peter_full_annovar_rsids())]
    return pd.concat([_martin_common(), group_a])


def _fetch_common_snps_plus_absent_from_peter_annovar() -> pd.DataFrame:
    # common + Martin-unique SNPs ABSENT from Peter's annovar output entirely. These drive the
    # collapse (found to be mostly DecodeME-QC-failed variants Martin retained; see
    # experiments/claude/group_b_root_cause).
    uniq = _martin_unique()
    group_b = uniq[~uniq["SNP"].isin(_peter_full_annovar_rsids())]
    return pd.concat([_martin_common(), group_b])


DATAFRAMES = {
    "peter": _fetch_peter,
    "martin": _fetch_martin,
    "martin_restricted_to_peter_snps": _fetch_martin_restricted_to_peter_snps,
    "common_snps_plus_present_in_peter_annovar": _fetch_common_snps_plus_present_in_peter_annovar,
    "common_snps_plus_absent_from_peter_annovar": _fetch_common_snps_plus_absent_from_peter_annovar,
}


# --------------------------------------------------------------------------------------------------
# The harness
# --------------------------------------------------------------------------------------------------
def run_one(
    dataframe_name: str, dataset_label: str, baseline_label: str
) -> pd.DataFrame:
    """
    Run gwaslab cts on DATAFRAMES[dataframe_name] with the given cts dataset and baseline model.

    The triple (dataframe_name, dataset_label, baseline_label) uniquely determines the input and
    the output path, so runs never collide.
    """
    cts_dir, index_name, headline = CTS_DATASETS[dataset_label]
    df = DATAFRAMES[dataframe_name]()
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
    out_path = OUT_DIR / f"{dataframe_name}_{dataset_label}_baseline_{baseline_label}.csv"
    summary.to_csv(out_path, index=False)
    head_p = summary.loc[summary["Name"] == headline, "Coefficient_P_value"].iloc[0]
    print(
        f"\n==> [{dataframe_name} | {dataset_label} / baseline {baseline_label}] "
        f"n_snps_input={len(df)}  {headline} p = {head_p:.4e}  (saved {out_path.name})"
    )
    return summary


def summary_to_markdown(table: pd.DataFrame) -> str:
    """Render the (dataframe x baseline) headline-p-value summary table as markdown."""
    return table.to_markdown(floatfmt=".3e")


def run_full_experiment(
    dataset_label: str = "gtex_brain",
    baseline_labels: tuple[str, ...] = ("v1.2", "v1.0"),
) -> pd.DataFrame:
    """
    Run EVERY dataframe variation (Peter's + Martin's) x every baseline through the same harness,
    then print a single controlled summary table of the headline-tissue p-value.
    """
    _, _, headline = CTS_DATASETS[dataset_label]
    rows = []
    for dataframe_name in DATAFRAMES:
        for baseline_label in baseline_labels:
            summary = run_one(dataframe_name, dataset_label, baseline_label)
            p = summary.loc[summary["Name"] == headline, "Coefficient_P_value"].iloc[0]
            rows.append(
                {
                    "dataframe": dataframe_name,
                    "baseline": baseline_label,
                    f"{headline}_p": p,
                }
            )
    table = pd.DataFrame(rows).pivot(
        index="dataframe", columns="baseline", values=f"{headline}_p"
    )
    print(f"\n### {dataset_label}: {headline} p-value\n")
    print(summary_to_markdown(table))
    table.to_csv(OUT_DIR / f"SUMMARY_{dataset_label}.csv")
    return table


def main():
    args = sys.argv[1:]
    if not args or args[0] == "full":
        run_full_experiment(*(args[1:]))
    else:
        # run.py <dataframe_name> <dataset_label> <baseline_label>
        run_one(args[0], args[1], args[2])


if __name__ == "__main__":
    main()
