"""List SLC6A2 missense/LoF variants with gnomAD v4 joint allele frequency at or above MIN_AF.

These are the coding variants with enough carriers for single-variant PheWAS lookups.

Run: pixi r python experiments/claude/net_a457p_gnomad/common_coding_variants.py
"""

import requests

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
MIN_AF = 1e-3
CODING_CONSEQUENCES = {
    "missense_variant",
    "stop_gained",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "start_lost",
}

QUERY = """
query {
  gene(gene_symbol: "SLC6A2", reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id rsids hgvsp consequence
      joint { ac an homozygote_count }
    }
  }
}
"""


def main() -> None:
    response = requests.post(GNOMAD_API, json={"query": QUERY}, timeout=300)
    response.raise_for_status()
    variants = response.json()["data"]["gene"]["variants"]
    rows = []
    for v in variants:
        joint = v["joint"]
        if v["consequence"] not in CODING_CONSEQUENCES or not joint or not joint["an"]:
            continue
        af = joint["ac"] / joint["an"]
        if af >= MIN_AF:
            rows.append((af, v))
    for af, v in sorted(rows, key=lambda r: -r[0]):
        print(
            f"{v['variant_id']:22} {v['rsids']} {v['hgvsp']:14} {v['consequence']:18} "
            f"AF={af:.4g} hom={v['joint']['homozygote_count']}"
        )


if __name__ == "__main__":
    main()
