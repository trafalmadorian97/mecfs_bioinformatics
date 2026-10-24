"""Print a compact summary of phewas_slc6a2.json (output of phewas_slc6a2.py).

Run: pixi r python experiments/claude/net_a457p_gnomad/summarize_phewas.py
"""

import json
from collections import Counter
from pathlib import Path

RESULTS_PATH = Path(__file__).parent / "phewas_slc6a2.json"
H4_THRESHOLD = 0.8
CLPP_THRESHOLD = 0.01


def p_value(row: dict) -> str:
    return f"{row['pValueMantissa']:.2g}e{row['pValueExponent']}"


def is_strong_coloc(coloc: dict) -> bool:
    if coloc["h4"] is not None:
        return coloc["h4"] >= H4_THRESHOLD
    return coloc["clpp"] is not None and coloc["clpp"] >= CLPP_THRESHOLD


def describe_credible_set(cs: dict) -> str:
    study = cs["study"]
    biosample = (study.get("biosample") or {}).get("biosampleName")
    variant = cs["variant"]
    consequence = (variant.get("mostSevereConsequence") or {}).get("label")
    l2g = [
        f"{r['target']['approvedSymbol']}:{r['score']:.2f}"
        for r in cs["l2GPredictions"]["rows"]
    ]
    return (
        f"{cs['studyType']:6} {study['projectId']:<14} {str(biosample or study['traitFromSource'])[:45]:45} "
        f"lead={variant['id']} {variant['rsIds']} {consequence} p={p_value(cs)} beta={cs['beta']} "
        f"qtlGene={cs['qtlGeneId']} trans={cs['isTransQtl']} L2G={l2g}"
    )


def print_credible_sets(credible_sets: list[dict]) -> None:
    print("Credible sets by type:", Counter(cs["studyType"] for cs in credible_sets))
    for cs in sorted(
        credible_sets, key=lambda c: (c["studyType"], c["pValueExponent"])
    ):
        print(describe_credible_set(cs))
        for coloc in cs["colocalisation"]["rows"]:
            if not is_strong_coloc(coloc):
                continue
            other = coloc["otherStudyLocus"]
            print(
                f"     COLOC {coloc['colocalisationMethod']} h4={coloc['h4']} clpp={coloc['clpp']} "
                f"sign={coloc['betaRatioSignAverage']} -> {other['study']['traitFromSource']} "
                f"({other['study']['projectId']}, {other['study']['publicationFirstAuthor']}) "
                f"lead={other['variant']['id']} p={p_value(other)} beta={other['beta']}"
            )


def print_genetic_evidence(evidence: list[dict]) -> None:
    print(
        "\nGenetic evidence rows by datasource:",
        Counter(e["datasourceId"] for e in evidence),
    )
    for e in sorted(
        evidence, key=lambda r: (r["datasourceId"], r["pValueExponent"] or 0)
    ):
        p = p_value(e) if e["pValueMantissa"] is not None else None
        if e["datasourceId"] == "gwas_credible_sets":
            cs = e["credibleSet"]
            study = cs["study"]
            print(
                f"  gwas L2G={e['score']:.2f} {study['traitFromSource'][:50]:50} {study['publicationFirstAuthor']} "
                f"N={study['nSamples']} cases={study['nCases']} lead={cs['variant']['id']} {cs['variant']['rsIds']} "
                f"p={p_value(cs)} beta={cs['beta']} EAF={cs['effectAlleleFrequencyFromSource']} conf={cs['confidence']}"
            )
            continue
        print(
            f"  {e['datasourceId']:14} {e['disease']['name'][:40]:40} src={str(e['diseaseFromSource'])[:40]:40} "
            f"proj={e['projectId']} method={e['statisticalMethod']} p={p} beta={e['beta']} OR={e['oddsRatio']} "
            f"carriers={e['studyCasesWithQualifyingVariants']}/{e['studyCases']} N={e['studySampleSize']} "
            f"var={e['variantRsId']} clinsig={e['clinicalSignificances']}"
        )


def print_gtex(eqtls: list[dict]) -> None:
    print(
        "\nGTEx v8 significant SLC6A2 cis-eQTLs per tissue (n variants, top variant):"
    )
    by_tissue: dict[str, list[dict]] = {}
    for row in eqtls:
        by_tissue.setdefault(row["tissueSiteDetailId"], []).append(row)
    for tissue, rows in sorted(by_tissue.items(), key=lambda kv: -len(kv[1])):
        top = min(rows, key=lambda r: r["pValue"])
        print(
            f"  {tissue:40} n={len(rows):4} top={top['snpId']} {top['variantId']} p={top['pValue']:.2g} NES={top['nes']:.2f}"
        )


def main() -> None:
    results = json.loads(RESULTS_PATH.read_text())
    print_credible_sets(results["credible_sets"])
    print_genetic_evidence(results["genetic_evidence"])
    print_gtex(results["gtex_v8_eqtls"])
    print("\nMouse phenotypes:")
    for mp in results["mouse_phenotypes"]:
        print(
            f"  {mp['modelPhenotypeLabel']} | {[m['allelicComposition'] for m in mp['biologicalModels']][:2]}"
        )
    print("\nTop associated diseases (overall score, datasource scores):")
    for row in sorted(results["associated_diseases"], key=lambda r: -r["score"])[:40]:
        sources = {d["id"]: round(d["score"], 2) for d in row["datasourceScores"]}
        print(f"  {row['score']:.3f} {row['disease']['name'][:50]:50} {sources}")


if __name__ == "__main__":
    main()
