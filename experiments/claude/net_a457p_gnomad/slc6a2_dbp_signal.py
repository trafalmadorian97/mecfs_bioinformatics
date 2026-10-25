"""Characterize the diastolic blood pressure (Keaton et al.) signal at SLC6A2 found by slc6a2_region_lookup_repo_gwas.py.

- Top DBP SNPs in SLC6A2 +/- 50 kb (GRCh37) from the repo asset.
- Whether the top SNPs are SLC6A2 cis-eQTLs in GTEx v8 (direction of the DBP-raising allele on expression).
- 1000 Genomes EUR LD (Ensembl REST) between the DBP lead and nearby candidates (rs5569, eQTL leads).
- Open Targets GWAS credible sets containing the DBP lead (other blood-pressure GWAS).

Run: pixi r python experiments/claude/net_a457p_gnomad/slc6a2_dbp_signal.py
"""

import json
from pathlib import Path

import polars as pl
import requests

from experiments.claude.net_a457p_gnomad.slc6a2_region_lookup_repo_gwas import (
    NEAR_FLANK,
    SLC6A2_END_37,
    SLC6A2_START_37,
    main_checkout_meta_to_path,
    region_snps,
)
from mecfs_bio.assets.gwas.blood_pressure.keaton_et_al_diastolic.analysis.keaton_dbp_standard_analysis import (
    KEATON_DBP_STANDARD_ANALYSIS,
)

OUT_DIR = Path(__file__).parent
GRCH37_TO_38_OFFSET = -33_912
"""Constant within this region (checked against rs9937434: b37 55,693,855, b38 55,659,943)."""
ENSEMBL_LD = "https://grch37.rest.ensembl.org/ld/human/pairwise/{a}/{b}?population_name=1000GENOMES:phase_3:EUR"
GTEX_EQTL = "https://gtexportal.org/api/v2/association/singleTissueEqtl"
OT_API = "https://api.platform.opentargets.org/api/v4/graphql"
LD_PARTNERS = (
    "rs5569",
    "rs9937434",
    "rs17841327",
    "rs187714",
    "rs55830311",
    "rs3785151",
)
N_TOP = 15


def top_dbp_snps() -> pl.DataFrame:
    table_path = main_checkout_meta_to_path()(
        KEATON_DBP_STANDARD_ANALYSIS.magma_tasks.parquet_file_task.meta
    )
    snps = region_snps(table_path)
    return (
        snps.filter(
            pl.col("POS").is_between(
                SLC6A2_START_37 - NEAR_FLANK, SLC6A2_END_37 + NEAR_FLANK
            )
        )
        .sort("P")
        .head(N_TOP)
    )


def gtex_eqtls_by_rsid() -> dict[str, list[dict]]:
    response = requests.get(
        GTEX_EQTL,
        params={"gencodeId": "ENSG00000103546.18", "datasetId": "gtex_v8"},
        timeout=300,
    )
    response.raise_for_status()
    by_rsid: dict[str, list[dict]] = {}
    for row in response.json()["data"]:
        by_rsid.setdefault(row["snpId"], []).append(row)
    return by_rsid


def ld_r2(a: str, b: str) -> float | None:
    response = requests.get(
        ENSEMBL_LD.format(a=a, b=b),
        headers={"Content-Type": "application/json"},
        timeout=120,
    )
    if not response.ok:
        return None
    rows = response.json()
    return float(rows[0]["r2"]) if rows else None


OT_QUERY = """
query($variantIds: [String!]!) {
  credibleSets(variantIds: $variantIds, studyTypes: [gwas], page: {index: 0, size: 50}) {
    count
    rows {
      pValueMantissa pValueExponent beta
      variant { id rsIds }
      study { traitFromSource publicationFirstAuthor nSamples }
      locus(variantIds: $variantIds) { rows { posteriorProbability } }
      l2GPredictions(page: {index: 0, size: 2}) { rows { score target { approvedSymbol } } }
    }
  }
}
"""


def ot_credible_sets(variant_id_38: str) -> list[dict]:
    response = requests.post(
        OT_API,
        json={"query": OT_QUERY, "variables": {"variantIds": [variant_id_38]}},
        timeout=300,
    )
    assert response.ok, response.text[:1000]
    return response.json()["data"]["credibleSets"]["rows"]


def main() -> None:
    top = top_dbp_snps()
    eqtls = gtex_eqtls_by_rsid()
    lead = top.row(0, named=True)
    print(f"DBP lead: {lead}")
    print("\nTop DBP SNPs (SLC6A2 +/- 50 kb, GRCh37) and their GTEx v8 SLC6A2 eQTLs:")
    for row in top.iter_rows(named=True):
        hits = eqtls.get(row["rsID"], [])
        eqtl_text = "; ".join(
            f"{h['tissueSiteDetailId']} NES={h['nes']:.2f} (alt {h['variantId'].split('_')[3]}) p={h['pValue']:.1e}"
            for h in hits
        )
        print(
            f"  {row['rsID']:12} pos37={row['POS']} EA={row['EA']} NEA={row['NEA']} EAF={row['EAF']} "
            f"beta={row['BETA']} p={row['P']:.2e} | eQTL: {eqtl_text or 'none'}"
        )
    print("\n1000G EUR r2 with DBP lead:")
    ld = {partner: ld_r2(lead["rsID"], partner) for partner in LD_PARTNERS}
    for partner, r2 in ld.items():
        print(f"  {lead['rsID']} - {partner}: r2={r2}")
    variant_id_38 = f"16_{lead['POS'] + GRCH37_TO_38_OFFSET}_{lead['NEA']}_{lead['EA']}"
    sets = ot_credible_sets(variant_id_38)
    print(
        f"\nOpen Targets GWAS credible sets containing {variant_id_38} ({lead['rsID']}): {len(sets)}"
    )
    for cs in sorted(sets, key=lambda c: c["pValueExponent"]):
        pip = max(
            (m["posteriorProbability"] for m in cs["locus"]["rows"]), default=None
        )
        l2g = [
            f"{p['target']['approvedSymbol']}:{p['score']:.2f}"
            for p in cs["l2GPredictions"]["rows"]
        ]
        print(
            f"  {cs['study']['traitFromSource'][:60]:60} {cs['study']['publicationFirstAuthor']} N={cs['study']['nSamples']} "
            f"lead={cs['variant']['id']} p={cs['pValueMantissa']:.1f}e{cs['pValueExponent']} PIP(lead DBP SNP)={pip} L2G={l2g}"
        )
    (OUT_DIR / "slc6a2_dbp_signal.json").write_text(
        json.dumps(
            {
                "top_snps": top.to_dicts(),
                "ld_r2_with_lead": ld,
                "ot_credible_sets": sets,
            },
            indent=1,
        )
    )


if __name__ == "__main__":
    main()
