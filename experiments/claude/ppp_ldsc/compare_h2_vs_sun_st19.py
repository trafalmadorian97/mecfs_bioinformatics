"""
Compare our LDSC per-protein heritabilities against Sun et al. 2023 ST19.

ST19 gives a SNP-based heritability DECOMPOSITION per protein (proteins with pQTLs only):
    - "pQTL component": cis pQTLs / trans pQTLs / all pQTLs  (fine-mapped lead-variant variance)
    - "Polygenic component"                                  (residual polygenic term)
    - "Total heritability (TH)" = all pQTLs + polygenic

Our estimator is genome-wide LD-score regression with a free intercept, run on ~1.2M
HapMap3 SNPs, with GenomicSEM's per-trait chi^2 filter max(0.001*N, 80). At N ~ 34k the
threshold is ~80, so the very large cis-pQTL peaks (chi^2 up to thousands) are removed
before regression. LDSC with a free intercept also does not attribute single large-effect
loci to h^2. So a priori:
    - our h^2 should track Sun's POLYGENIC component most closely;
    - it should track TOTAL heritability poorly, because TH is dominated by the cis-pQTL
      variance that our filter/estimator deliberately does not capture;
    - cis_excluded should, if anything, be even more polygenic-like.

This script joins on OID and reports Pearson + Spearman correlations and summary stats for
each pairing, plus a couple of illustrative cis-dominated proteins.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import openpyxl
import polars as pl
from scipy.stats import pearsonr, spearmanr

REPO = Path(__file__).resolve().parents[3]
H2_PATH = (
    REPO
    / "assets/base_asset_store/gwas/ukbb_ppp/ppp_heritability/analysis"
    / "ppp_heritability_hapmap_3.parquet"
)
XLSX = (
    REPO
    / "assets/base_asset_store/reference_data/pqtl_data/sun_et_al_2023/raw"
    / "sun_et_al_2023_supplementary.xlsx"
)

# ST19 column offsets (0-indexed) under the merged header on row 2/3.
_C_PROTEIN_ID = 0
_C_CIS_PQTL = 1
_C_TRANS_PQTL = 2
_C_ALL_PQTL = 3
_C_POLYGENIC = 4
_C_TOTAL_H2 = 5


def _oid_from_protein_id(protein_id: str) -> str | None:
    """SIRPA:P78324:OID20304:v1 -> OID20304."""
    for part in protein_id.split(":"):
        if part.startswith("OID"):
            return part
    return None


def read_st19(path: Path) -> pl.DataFrame:
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb["ST19"]
    rows: list[dict] = []
    for row in ws.iter_rows(values_only=True):
        pid = row[_C_PROTEIN_ID]
        if not (isinstance(pid, str) and ":OID" in pid):
            continue
        oid = _oid_from_protein_id(pid)
        if oid is None:
            continue
        rows.append(
            {
                "oid": oid,
                "sun_cis_pqtl": row[_C_CIS_PQTL],
                "sun_trans_pqtl": row[_C_TRANS_PQTL],
                "sun_all_pqtl": row[_C_ALL_PQTL],
                "sun_polygenic": row[_C_POLYGENIC],
                "sun_total_h2": row[_C_TOTAL_H2],
            }
        )
    return pl.DataFrame(rows)


def _report(label: str, x: np.ndarray, y: np.ndarray) -> None:
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    pr = pearsonr(x, y)
    sr = spearmanr(x, y)
    print(
        f"{label:38s} n={x.size:4d}  "
        f"Pearson r={pr.statistic:+.3f} (p={pr.pvalue:.1e})  "
        f"Spearman rho={sr.statistic:+.3f} (p={sr.pvalue:.1e})"
    )


def main() -> None:
    h2 = pl.read_parquet(H2_PATH)
    st19 = read_st19(XLSX)
    print(f"our h2 rows: {h2.height}  (proteins: {h2['oid'].n_unique()})")
    print(f"ST19 rows:   {st19.height}")

    allv = h2.filter(pl.col("variant_set") == "all_variants").select(
        "oid", "gene", pl.col("h2").alias("our_h2_all"), "n_bar"
    )
    cisx = h2.filter(pl.col("variant_set") == "cis_excluded").select(
        "oid", pl.col("h2").alias("our_h2_cisx")
    )

    joined = allv.join(st19, on="oid", how="inner").join(cisx, on="oid", how="left")
    print(f"joined on OID: {joined.height} proteins\n")

    g = joined.to_pandas()
    print("=== our all_variants LDSC h^2 vs Sun ST19 ===")
    _report("all_variants vs Total heritability", g["our_h2_all"].values, g["sun_total_h2"].values)
    _report("all_variants vs Polygenic component", g["our_h2_all"].values, g["sun_polygenic"].values)
    _report("all_variants vs all-pQTL component", g["our_h2_all"].values, g["sun_all_pqtl"].values)
    _report("all_variants vs trans-pQTL component", g["our_h2_all"].values, g["sun_trans_pqtl"].values)
    print("\n=== our cis_excluded LDSC h^2 vs Sun ST19 ===")
    _report("cis_excluded vs Total heritability", g["our_h2_cisx"].values, g["sun_total_h2"].values)
    _report("cis_excluded vs Polygenic component", g["our_h2_cisx"].values, g["sun_polygenic"].values)
    _report("cis_excluded vs trans-pQTL component", g["our_h2_cisx"].values, g["sun_trans_pqtl"].values)
    _report("cis_excluded vs (polygenic+trans)", g["our_h2_cisx"].values, (g["sun_polygenic"] + g["sun_trans_pqtl"]).values)

    print("\n=== distributions ===")
    for col in ["our_h2_all", "our_h2_cisx", "sun_total_h2", "sun_polygenic", "sun_all_pqtl", "sun_trans_pqtl"]:
        v = g[col].values.astype(float)
        v = v[np.isfinite(v)]
        print(f"{col:16s} median={np.median(v):.4f}  mean={np.mean(v):.4f}  p10={np.percentile(v,10):.4f}  p90={np.percentile(v,90):.4f}")

    print("\n=== most cis-dominated proteins (Sun cis pQTL), our h^2 for contrast ===")
    top = joined.sort("sun_cis_pqtl", descending=True).head(8)
    print(top.select("gene", "sun_cis_pqtl", "sun_total_h2", "sun_polygenic", "our_h2_all", "our_h2_cisx"))

    print("\n=== highest polygenic-component proteins ===")
    topp = joined.sort("sun_polygenic", descending=True).head(8)
    print(topp.select("gene", "sun_polygenic", "sun_total_h2", "our_h2_all", "our_h2_cisx"))


if __name__ == "__main__":
    main()
