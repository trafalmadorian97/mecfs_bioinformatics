"""Full-results PheWAS lookups for SLC6A2 (NET): gene-level rare-variant tests and single-variant lookups.

a) Coding variation:
   - Genebass (UK Biobank 394k exomes, SAIGE-GENE): gene-level SKAT-O / burden / SKAT p-values for SLC6A2 across
     ~4,500 phenotypes, for the pLoF and missense|LC burden sets, with synonymous as a negative control.
   - FinnGen R12 (~500k Finns, GRCh38) single-variant PheWAS for the SLC6A2 coding variants with gnomAD v4 AF >= 0.1%
     (from common_coding_variants.py).
b) Expression:
   - FinnGen R12 single-variant PheWAS for lead SLC6A2 cis-eQTL / transcript-usage QTL variants (GTEx v8/v10 via the
     GTEx API and Open Targets credible sets; see phewas_slc6a2.py).

Writes genebass_slc6a2.parquet and finngen_slc6a2_lookups.parquet next to this script.

Run: pixi r python experiments/claude/net_a457p_gnomad/phewas_lookups.py
"""

from pathlib import Path

import attrs
import polars as pl
import requests

OUT_DIR = Path(__file__).parent
GENEBASS_API = "https://main.genebass.org/api"
FINNGEN_API = "https://r12.finngen.fi/api"
OT_API = "https://api.platform.opentargets.org/api/v4/graphql"
SLC6A2_ENSG = "ENSG00000103546"
GENEBASS_BURDEN_SETS = ("pLoF", "missense|LC", "synonymous")
GENEBASS_PHENO_KEY = ["trait_type", "phenocode", "pheno_sex", "coding", "modifier"]
TIMEOUT_SECONDS = 900


@attrs.frozen(slots=True)
class LookupVariant:
    variant_id: str
    """GRCh38 chrom-pos-ref-alt."""
    rsid: str
    role: str


LOOKUP_VARIANTS = (
    LookupVariant(
        "16-55659943-G-C",
        "rs9937434",
        "eQTL lead: adrenal gland, skin (GTEx v8); alt lowers expression",
    ),
    LookupVariant(
        "16-55660341-A-C",
        "rs17841327",
        "eQTL lead: skin (GTEx v10, TwinsUK); alt lowers expression",
    ),
    LookupVariant(
        "16-55672587-A-G",
        "rs187714",
        "eQTL lead: testis (GTEx v10); colocalises with birth weight",
    ),
    LookupVariant(
        "16-55648759-A-C",
        "rs55830311",
        "eQTL lead: esophagus (GTEx v10); colocalises with height, lean mass",
    ),
    LookupVariant(
        "16-55703744-G-A",
        "rs42460",
        "tuQTL lead: testis; colocalises with placental weight",
    ),
    LookupVariant(
        "16-55678607-G-C",
        "rs3785151",
        "GWAS lead: MVP sensorineural hearing loss (L2G SLC6A2 0.83); weak skin eQTL",
    ),
    LookupVariant("16-55669586-C-T", "rs1805065", "coding: p.Thr99Ile, gnomAD AF 1.7%"),
    LookupVariant(
        "16-55685231-G-A", "rs1805066", "coding: p.Val245Ile, gnomAD AF 0.9%"
    ),
    LookupVariant(
        "16-55671926-G-A", "rs11568343", "coding: p.Cys27Tyr, gnomAD AF 0.8%"
    ),
    LookupVariant(
        "16-55671848-T-G", "rs552629380", "coding: p.Met1? start lost, gnomAD AF 0.13%"
    ),
    LookupVariant(
        "16-55691982-C-G", "rs45564432", "coding: p.Thr283Arg, gnomAD AF 0.1%"
    ),
)


def get_json(url: str) -> dict | list:
    response = requests.get(url, timeout=TIMEOUT_SECONDS)
    response.raise_for_status()
    return response.json()


def records_to_frame(records: list[dict]) -> pl.DataFrame:
    """Genebass encodes missing values as the string NA; map them to null before schema inference."""
    cleaned = [{k: (None if v == "NA" else v) for k, v in r.items()} for r in records]
    return pl.DataFrame(cleaned, infer_schema_length=None)


def fetch_genebass() -> pl.DataFrame:
    phenotypes = records_to_frame(get_json(f"{GENEBASS_API}/phenotypes")).select(
        *GENEBASS_PHENO_KEY,
        "description",
        "coding_description",
        "category",
        "n_cases",
        "n_controls",
    )
    results = pl.concat(
        [
            records_to_frame(
                get_json(f"{GENEBASS_API}/phewas/{SLC6A2_ENSG}?burdenSet={burden_set}")[
                    "phewas"
                ]
            ).select(
                *GENEBASS_PHENO_KEY,
                "annotation",
                "Pvalue",
                "Pvalue_Burden",
                "Pvalue_SKAT",
                "BETA_Burden",
                "total_variants",
            )
            for burden_set in GENEBASS_BURDEN_SETS
        ],
        how="vertical_relaxed",
    )
    return results.join(phenotypes, on=GENEBASS_PHENO_KEY, how="left")


def fetch_finngen_variant(variant: LookupVariant) -> pl.DataFrame | None:
    """None when the variant is not in the FinnGen release (the API returns 404)."""
    response = requests.get(
        f"{FINNGEN_API}/variant/{variant.variant_id}", timeout=TIMEOUT_SECONDS
    )
    if response.status_code == 404:
        print(f"not in FinnGen R12: {variant.rsid} {variant.variant_id}")
        return None
    response.raise_for_status()
    payload = response.json()
    assert payload["results"], f"no FinnGen results for {variant.variant_id}"
    return (
        records_to_frame(payload["results"])
        .select(
            "phenocode",
            "phenostring",
            "category",
            "beta",
            "sebeta",
            "pval",
            "maf",
            "n_case",
            "n_control",
        )
        .with_columns(
            pl.lit(variant.variant_id).alias("variant_id"),
            pl.lit(variant.rsid).alias("rsid"),
            pl.lit(variant.role).alias("role"),
        )
    )


OT_CREDIBLE_SETS_QUERY = """
query($variantIds: [String!]!, $index: Int!) {
  credibleSets(variantIds: $variantIds, studyTypes: [gwas], page: {index: $index, size: 50}) {
    count
    rows {
      studyLocusId pValueMantissa pValueExponent beta confidence
      variant { id rsIds }
      study { id projectId traitFromSource publicationFirstAuthor nSamples nCases }
      locus(variantIds: $variantIds) { rows { posteriorProbability variant { id } } }
      l2GPredictions(page: {index: 0, size: 1}) { rows { score target { approvedSymbol } } }
    }
  }
}
"""


def fetch_ot_credible_sets(variant: LookupVariant) -> pl.DataFrame:
    """GWAS credible sets (any study in Open Targets) whose 95% set contains the variant."""
    ot_variant_id = variant.variant_id.replace("-", "_")
    rows: list[dict] = []
    index = 0
    while True:
        response = requests.post(
            OT_API,
            json={
                "query": OT_CREDIBLE_SETS_QUERY,
                "variables": {"variantIds": [ot_variant_id], "index": index},
            },
            timeout=TIMEOUT_SECONDS,
        )
        assert response.ok, response.text[:2000]
        page = response.json()["data"]["credibleSets"]
        rows.extend(page["rows"])
        if len(rows) >= page["count"] or not page["rows"]:
            break
        index += 1
    return pl.DataFrame(
        [
            {
                "variant_id": variant.variant_id,
                "rsid": variant.rsid,
                "role": variant.role,
                "trait": r["study"]["traitFromSource"],
                "study_id": r["study"]["id"],
                "project": r["study"]["projectId"],
                "first_author": r["study"]["publicationFirstAuthor"],
                "n_samples": r["study"]["nSamples"],
                "lead_variant": r["variant"]["id"],
                "lead_p": r["pValueMantissa"] * 10.0 ** r["pValueExponent"],
                "lead_beta": r["beta"],
                "confidence": r["confidence"],
                "lookup_variant_pip": max(
                    (m["posteriorProbability"] for m in r["locus"]["rows"]),
                    default=None,
                ),
                "l2g_top_gene": next(
                    (
                        p["target"]["approvedSymbol"]
                        for p in r["l2GPredictions"]["rows"]
                    ),
                    None,
                ),
            }
            for r in rows
        ],
        infer_schema_length=None,
    )


def main() -> None:
    genebass = fetch_genebass()
    genebass.write_parquet(OUT_DIR / "genebass_slc6a2.parquet")
    print("genebass rows:", genebass.height)
    finngen_frames = [fetch_finngen_variant(v) for v in LOOKUP_VARIANTS]
    finngen = pl.concat(
        [f for f in finngen_frames if f is not None], how="vertical_relaxed"
    )
    finngen.write_parquet(OUT_DIR / "finngen_slc6a2_lookups.parquet")
    print("finngen rows:", finngen.height)
    ot_sets = pl.concat(
        [fetch_ot_credible_sets(v) for v in LOOKUP_VARIANTS], how="diagonal_relaxed"
    )
    ot_sets.write_parquet(OUT_DIR / "ot_gwas_credible_sets_containing_lookups.parquet")
    print("open targets credible-set rows:", ot_sets.height)


if __name__ == "__main__":
    main()
