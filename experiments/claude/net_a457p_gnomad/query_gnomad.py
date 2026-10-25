"""Query the gnomAD GraphQL API for the SLC6A2 (NET) A457P variant from Shannon et al. 2000.

Finds p.Ala457Pro among SLC6A2 variants in gnomAD v4 and v2, then fetches the
variant-level detail (frequencies, in silico predictors, ClinVar) and gene constraint.

Run: pixi r python experiments/claude/net_a457p_gnomad/query_gnomad.py
"""

import json

import requests

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
GENE = "SLC6A2"
TARGET_HGVSP = "p.Ala457Pro"


def run_query(query: str, variables: dict[str, str]) -> dict:
    response = requests.post(
        GNOMAD_API, json={"query": query, "variables": variables}, timeout=120
    )
    response.raise_for_status()
    payload = response.json()
    assert "errors" not in payload, payload["errors"]
    return payload["data"]


GENE_VARIANTS_QUERY = """
query($gene: String!, $dataset: DatasetId!, $ref: ReferenceGenomeId!) {
  gene(gene_symbol: $gene, reference_genome: $ref) {
    canonical_transcript_id
    mane_select_transcript { ensembl_id refseq_id }
    gnomad_constraint { exp_lof obs_lof oe_lof oe_lof_lower oe_lof_upper pLI lof_z
                        exp_mis obs_mis oe_mis oe_mis_lower oe_mis_upper mis_z syn_z }
    variants(dataset: $dataset) {
      variant_id rsids hgvsp hgvsc consequence transcript_id
      exome { ac an homozygote_count }
      genome { ac an homozygote_count }
    }
  }
}
"""

VARIANT_QUERY = """
query($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id rsids caid
    exome { ac an homozygote_count hemizygote_count filters
            populations { id ac an homozygote_count } }
    genome { ac an homozygote_count filters
             populations { id ac an homozygote_count } }
    joint { ac an homozygote_count filters }
    in_silico_predictors { id value flags }
    transcript_consequences { gene_symbol transcript_id hgvsc hgvsp major_consequence
                              is_canonical is_mane_select polyphen_prediction sift_prediction }
  }
  clinvar_variant(variant_id: $variantId, reference_genome: GRCh38) {
    clinical_significance review_status clinvar_variation_id gold_stars
    submissions { clinical_significance conditions { name } last_evaluated review_status submitter_name }
  }
}
"""


def main() -> None:
    results: dict = {}
    for dataset, ref in [("gnomad_r4", "GRCh38"), ("gnomad_r2_1", "GRCh37")]:
        gene = run_query(
            GENE_VARIANTS_QUERY, {"gene": GENE, "dataset": dataset, "ref": ref}
        )["gene"]
        variants = gene["variants"]
        results[dataset] = {
            "canonical": gene["canonical_transcript_id"],
            "mane": gene["mane_select_transcript"],
            "constraint": gene["gnomad_constraint"],
            "n_variants": len(variants),
            "codon_457_variants": [
                v for v in variants if v["hgvsp"] and "Ala457" in v["hgvsp"]
            ],
        }
    results["a457p_in_gnomad_r4"] = any(
        TARGET_HGVSP in v["hgvsp"] for v in results["gnomad_r4"]["codon_457_variants"]
    )
    results["codon_457_detail_r4"] = {
        v["variant_id"]: run_query(
            VARIANT_QUERY, {"variantId": v["variant_id"], "dataset": "gnomad_r4"}
        )
        for v in results["gnomad_r4"]["codon_457_variants"]
        if v["consequence"] == "missense_variant"
    }
    results["a457p_dbnsfp"] = fetch_dbnsfp(A457P_HG38_HGVS)
    print(json.dumps(results, indent=2))


A457P_HG38_HGVS = "chr16:g.55698005G>C"
"""c.1369G>C on NM_001172501 (MANE); GCx -> CCx, Ala457Pro."""

DBNSFP_FIELDS = (
    "dbnsfp.aa,dbnsfp.cadd,dbnsfp.revel,dbnsfp.alphamissense,dbnsfp.primateai,"
    "dbnsfp.esm1b,dbnsfp.polyphen2,dbnsfp.sift,dbnsfp.mutationtaster,dbnsfp.metarnn,"
    "dbnsfp.bayesdel,dbnsfp.gerp++,dbnsfp.phylop,dbnsfp.phastcons,dbnsfp.clinpred,"
    "dbnsfp.mvp,dbnsfp.varity,cadd.phred,clinvar,dbsnp.rsid,gnomad_exome,gnomad_genome"
)


def fetch_dbnsfp(hgvs: str) -> dict:
    response = requests.get(
        f"https://myvariant.info/v1/variant/{hgvs}",
        params={"assembly": "hg38", "fields": DBNSFP_FIELDS},
        timeout=120,
    )
    response.raise_for_status()
    return response.json()


if __name__ == "__main__":
    main()
