"""Print a compact summary of query_gnomad.json (output of query_gnomad.py).

Run: pixi r python experiments/claude/net_a457p_gnomad/summarize.py
"""

import json
from pathlib import Path

RESULTS_PATH = Path(__file__).parent / "query_gnomad.json"


def main() -> None:
    results = json.loads(RESULTS_PATH.read_text())
    print("A457P present in gnomAD v4:", results["a457p_in_gnomad_r4"])
    for variant_id, detail in results["codon_457_detail_r4"].items():
        variant = detail["variant"]
        hgvsp = variant["transcript_consequences"][0]["hgvsp"]
        print(f"\n== {variant_id} {hgvsp} joint AC={variant['joint']['ac']}")
        for block in ("exome", "genome"):
            if variant[block]:
                carriers = {
                    p["id"]: p["ac"]
                    for p in variant[block]["populations"]
                    if p["ac"] > 0
                    and "_" not in p["id"]
                    and p["id"] not in ("XX", "XY")
                }
                print(f"  {block} AC by ancestry: {carriers}")
        print(
            "  predictors:",
            {p["id"]: p["value"] for p in variant["in_silico_predictors"]},
        )
        clinvar = detail["clinvar_variant"]
        if clinvar:
            print(
                "  ClinVar:",
                clinvar["clinical_significance"],
                "|",
                clinvar["review_status"],
                "| VCV",
                clinvar["clinvar_variation_id"],
            )
            for submission in clinvar["submissions"]:
                print(
                    "    -",
                    submission["submitter_name"],
                    submission["clinical_significance"],
                    [c["name"] for c in submission["conditions"]],
                    submission["last_evaluated"],
                )
    print("\n== A457P dbNSFP (myvariant.info)")
    print(json.dumps(results["a457p_dbnsfp"], indent=1))


if __name__ == "__main__":
    main()
