"""Simple PheWAS-style survey of phenotypes linked to SLC6A2 (NET) via coding and regulatory variation.

Sources:
- Open Targets Platform GraphQL: credible sets attached to SLC6A2 (eQTL/pQTL/sQTL sets whose qtlGeneId is
  SLC6A2, plus GWAS sets that L2G assigns to it), their colocalisation with GWAS credible sets, gene-burden
  evidence (UK Biobank exome rare-variant collapsing tests from Genebass, AstraZeneca PheWAS, Regeneron),
  and mouse knockout phenotypes.
- GTEx v8 API: significant single-tissue cis-eQTLs for SLC6A2.

Run: pixi r python experiments/claude/net_a457p_gnomad/phewas_slc6a2.py > experiments/claude/net_a457p_gnomad/phewas_slc6a2.json
"""

import json

import requests

OT_API = "https://api.platform.opentargets.org/api/v4/graphql"
GTEX_API = "https://gtexportal.org/api/v2"
SLC6A2_ENSG = "ENSG00000103546"
SLC6A2_GENCODE = "ENSG00000103546.18"
PAGE_SIZE = 10


def ot_query(query: str, variables: dict) -> dict:
    response = requests.post(
        OT_API, json={"query": query, "variables": variables}, timeout=300
    )
    assert response.ok, response.text[:2000]
    payload = response.json()
    assert "errors" not in payload, payload["errors"]
    return payload["data"]


CREDIBLE_SETS_QUERY = """
query($id: String!, $index: Int!, $size: Int!) {
  target(ensemblId: $id) {
    credibleSets(page: {index: $index, size: $size}) {
      count
      rows {
        studyLocusId studyType qtlGeneId isTransQtl confidence
        pValueMantissa pValueExponent beta
        variant { id rsIds mostSevereConsequence { label } }
        study { id projectId traitFromSource biosample { biosampleName } }
        l2GPredictions(page: {index: 0, size: 3}) { rows { score target { approvedSymbol } } }
        colocalisation(studyTypes: [gwas], page: {index: 0, size: 50}) {
          count
          rows {
            colocalisationMethod h4 clpp betaRatioSignAverage
            otherStudyLocus {
              studyLocusId pValueMantissa pValueExponent beta
              variant { id rsIds }
              study { id projectId traitFromSource publicationFirstAuthor nSamples }
            }
          }
        }
      }
    }
  }
}
"""

ASSOCIATIONS_QUERY = """
query($id: String!, $index: Int!, $size: Int!) {
  target(ensemblId: $id) {
    mousePhenotypes { modelPhenotypeLabel modelPhenotypeClasses { label }
                      biologicalModels { allelicComposition geneticBackground } }
    associatedDiseases(page: {index: $index, size: $size}) {
      count
      rows { score disease { id name } datasourceScores { id score } }
    }
  }
}
"""

EVIDENCE_QUERY = """
query($id: String!, $efoIds: [String!]!) {
  target(ensemblId: $id) {
    evidences(efoIds: $efoIds, datasourceIds: ["gene_burden", "gwas_credible_sets", "eva", "orphanet",
                                                "gene2phenotype", "clingen", "uniprot_variants"], size: 500) {
      count
      rows {
        datasourceId disease { id name } diseaseFromSource score
        projectId cohortId statisticalMethod statisticalMethodOverview
        pValueMantissa pValueExponent beta oddsRatio
        studyCasesWithQualifyingVariants studyCases studySampleSize
        variant { id rsIds } variantRsId clinicalSignificances literature
        credibleSet { studyLocusId pValueMantissa pValueExponent beta effectAlleleFrequencyFromSource confidence
                      study { traitFromSource projectId publicationFirstAuthor nSamples nCases }
                      variant { id rsIds } }
      }
    }
  }
}
"""


def fetch_paged(query: str, field: str) -> list[dict]:
    rows: list[dict] = []
    index = 0
    while True:
        page = ot_query(query, {"id": SLC6A2_ENSG, "index": index, "size": PAGE_SIZE})[
            "target"
        ]
        rows.extend(page[field]["rows"])
        if len(rows) >= page[field]["count"] or not page[field]["rows"]:
            return rows
        index += 1


def fetch_associations() -> dict:
    first = ot_query(ASSOCIATIONS_QUERY, {"id": SLC6A2_ENSG, "index": 0, "size": 0})[
        "target"
    ]
    return {
        "mouse_phenotypes": first["mousePhenotypes"],
        "associated_diseases": fetch_paged(ASSOCIATIONS_QUERY, "associatedDiseases"),
    }


def fetch_genetic_evidence(efo_ids: list[str]) -> list[dict]:
    return ot_query(EVIDENCE_QUERY, {"id": SLC6A2_ENSG, "efoIds": efo_ids})["target"][
        "evidences"
    ]["rows"]


def fetch_gtex_eqtls() -> list[dict]:
    response = requests.get(
        f"{GTEX_API}/association/singleTissueEqtl",
        params={"gencodeId": SLC6A2_GENCODE, "datasetId": "gtex_v8"},
        timeout=300,
    )
    response.raise_for_status()
    return response.json()["data"]


def main() -> None:
    associations = fetch_associations()
    efo_ids = [row["disease"]["id"] for row in associations["associated_diseases"]]
    results = {
        "credible_sets": fetch_paged(CREDIBLE_SETS_QUERY, "credibleSets"),
        **associations,
        "genetic_evidence": fetch_genetic_evidence(efo_ids),
        "gtex_v8_eqtls": fetch_gtex_eqtls(),
    }
    print(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
