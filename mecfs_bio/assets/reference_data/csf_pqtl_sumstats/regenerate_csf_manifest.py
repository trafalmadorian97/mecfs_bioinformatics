"""
Regenerate the Western et al. 2024 CSF pQTL summary-statistics manifest.

The manifest is a static list mapping every published CSF aptamer to the GWAS
Catalog study accession whose GWAS-SSF file holds its summary statistics. It is the
CSF analogue of regenerate_ppp_manifest.py.

Two upstream sources are combined:

  1. The GWAS Catalog REST API, listing the 7,008 studies of PMID 39528825, each
     with an accession and a human-readable trait string.
  2. aptamer_info.xlsx (from the paper's Box deposit), the SomaScan 7k analyte
     table: 7,584 rows, of which the 576 with a non-null "Step Removed" are the
     unpublished aptamers, leaving exactly 7,008.

The trait string is the only link between the two. The Catalog assigns each study a
UNIQUE trait, having itself disambiguated shared target names by appending
"(analyte X####.##)". A four-rule resolver recovers a complete accession -> analyte
bijection; the resolver asserts the bijection, so a future upstream re-curation that
breaks the mapping fails loudly here rather than silently at build time.

Run via:
  pixi r python mecfs_bio/assets/reference_data/csf_pqtl_sumstats/regenerate_csf_manifest.py

Downloads no summary statistics; it fetches only study metadata (a few MB of JSON)
and the 826 KB analyte table. Re-run to pick up any upstream re-curation, then diff
the result against the committed manifest.
"""

import csv
import json
import re
import urllib.request
from pathlib import Path

import pandas as pd

_DIR = Path(__file__).parent
MANIFEST_PATH = _DIR / "csf_aptamer_manifest.csv"

PUBMED_ID = "39528825"
# GWAS Catalog REST endpoint listing every study of a publication, paged.
_REST_SEARCH_URL = (
    "https://www.ebi.ac.uk/gwas/rest/api/studies/search/"
    f"findByPublicationIdPubmedId?pubmedId={PUBMED_ID}&size=500"
)
EXPECTED_STUDY_COUNT = 7008

# aptamer_info.xlsx in the paper's Box deposit (shared folder CSF_Soma7K_pQTL). The
# file id is stable; the download endpoint needs only the shared_name and file id, no
# login. If Box ever renames the deposit, re-derive both from the folder page.
_BOX_SHARED_NAME = "3kms8fjz228qw7z8pti9ar3ivhyofewn"
_BOX_APTAMER_INFO_FILE_ID = "f_1748814237193"
_BOX_DOWNLOAD_URL = (
    "https://wustl.app.box.com/index.php?rm=box_download_shared_file"
    f"&shared_name={_BOX_SHARED_NAME}&file_id={_BOX_APTAMER_INFO_FILE_ID}"
)

# aptamer_info.xlsx columns.
_ANALYTES_COL = "Analytes"  # aptamer primary key, e.g. X13681.173
_SEQID_COL = "SeqId"  # e.g. 13681-173
_UNIPROT_COL = "UniProt"
_ENTREZ_SYMBOL_COL = "EntrezGeneSymbol"
_TARGET_FULL_NAME_COL = "TargetFullName"
_STEP_REMOVED_COL = "Step Removed"  # non-null => aptamer not published

# Rule 1: the Catalog disambiguates shared target names by appending this to the trait.
_ANALYTE_IN_TRAIT_PATTERN = re.compile(r"\(analyte (X[\d.]+)\)")

# Rule 4: three Casein kinase II traits the first three rules cannot resolve (the
# trait names do not match any TargetFullName + " levels"). Mapping verified by hand.
_TRAIT_OVERRIDES = {
    "Casein kinase II subunit alpha-2 levels": "X13681.173",
    "Casein kinase II alpha-1: beta heterotetramer levels": "X5225.50",
    "Casein kinase II alpha-2: beta heterotetramer levels": "X5226.36",
}

MANIFEST_COLUMNS = [
    "analyte",
    "seq_id",
    "uniprot",
    "entrez_gene_symbol",
    "target_full_name",
    "accession",
]


def fetch_studies() -> list[tuple[str, str]]:
    """Page the GWAS Catalog REST API, returning (accession, trait) for every study."""
    studies: list[tuple[str, str]] = []
    url: str | None = _REST_SEARCH_URL
    while url is not None:
        with urllib.request.urlopen(url) as response:
            payload = json.load(response)
        for study in payload["_embedded"]["studies"]:
            trait = study["diseaseTrait"]["trait"]
            studies.append((study["accessionId"], trait))
        url = payload["_links"].get("next", {}).get("href")
    return studies


def download_aptamer_info(dest: Path) -> Path:
    """Download aptamer_info.xlsx from the Box deposit."""
    urllib.request.urlretrieve(_BOX_DOWNLOAD_URL, dest)
    return dest


def published_aptamers(aptamer_info_path: Path) -> pd.DataFrame:
    """The 7,008 published aptamers (Step Removed is null), indexed by analyte."""
    frame = pd.read_excel(aptamer_info_path, sheet_name="Sheet1")
    published = frame[frame[_STEP_REMOVED_COL].isna()].copy()
    assert len(published) == EXPECTED_STUDY_COUNT, (
        f"expected {EXPECTED_STUDY_COUNT} published aptamers, got {len(published)}"
    )
    return published.set_index(_ANALYTES_COL, drop=False)


def build_resolver(published: pd.DataFrame):
    """Return a function trait -> analyte, closed over the published analyte table.

    Four rules, tried in order:
      1. an explicit "(analyte X####.##)" in the trait (the Catalog's own tiebreak);
      2. an exact "<TargetFullName> levels" match to a uniquely-named aptamer;
      3. the same, case-insensitively;
      4. a hardcoded override for three Casein kinase II traits.
    """
    exact: dict[str, list[str]] = {}
    lower: dict[str, list[str]] = {}
    for analyte, target in zip(
        published[_ANALYTES_COL], published[_TARGET_FULL_NAME_COL]
    ):
        key = f"{target} levels"
        exact.setdefault(key, []).append(analyte)
        lower.setdefault(key.lower(), []).append(analyte)

    def resolve(trait: str) -> tuple[str, str]:
        if trait in _TRAIT_OVERRIDES:
            return _TRAIT_OVERRIDES[trait], "override"
        match = _ANALYTE_IN_TRAIT_PATTERN.search(trait)
        if match is not None:
            return match.group(1), "analyte_in_trait"
        if len(exact.get(trait, [])) == 1:
            return exact[trait][0], "exact_name"
        if len(lower.get(trait.lower(), [])) == 1:
            return lower[trait.lower()][0], "case_insensitive_name"
        raise KeyError(trait)

    return resolve


def build_manifest_rows(
    studies: list[tuple[str, str]], published: pd.DataFrame
) -> list[dict]:
    """Resolve every study to an aptamer and assert a complete bijection."""
    resolve = build_resolver(published)
    published_analytes = set(published[_ANALYTES_COL])

    rows: list[dict] = []
    rule_counts: dict[str, int] = {}
    unresolved: list[tuple[str, str]] = []
    used: dict[str, str] = {}  # analyte -> accession that claimed it
    for accession, trait in studies:
        try:
            analyte, rule = resolve(trait)
        except KeyError:
            unresolved.append((accession, trait))
            continue
        assert analyte in published_analytes, (
            f"{accession}: resolved analyte {analyte} is not a published aptamer"
        )
        assert analyte not in used, (
            f"analyte {analyte} claimed by both {used[analyte]} and {accession}"
        )
        used[analyte] = accession
        rule_counts[rule] = rule_counts.get(rule, 0) + 1
        info = published.loc[analyte]
        rows.append(
            {
                "analyte": analyte,
                "seq_id": info[_SEQID_COL],
                "uniprot": info[_UNIPROT_COL],
                "entrez_gene_symbol": info[_ENTREZ_SYMBOL_COL],
                "target_full_name": info[_TARGET_FULL_NAME_COL],
                "accession": accession,
            }
        )

    assert not unresolved, f"{len(unresolved)} traits did not resolve: {unresolved[:5]}"
    missing = published_analytes - set(used)
    assert not missing, (
        f"{len(missing)} published aptamers were never hit: {list(missing)[:5]}"
    )
    assert len(rows) == EXPECTED_STUDY_COUNT, (
        f"expected {EXPECTED_STUDY_COUNT} rows, got {len(rows)}"
    )
    print("resolved by rule:", rule_counts)
    return rows


def write_manifest(rows: list[dict], manifest_path: Path) -> None:
    rows_sorted = sorted(rows, key=lambda row: row["analyte"])
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS)
        writer.writeheader()
        writer.writerows(rows_sorted)
    print(f"wrote {len(rows_sorted)} rows to {manifest_path}")


def main() -> None:
    aptamer_info_path = _DIR / "aptamer_info.xlsx"
    print("fetching study list from GWAS Catalog ...")
    studies = fetch_studies()
    assert len(studies) == EXPECTED_STUDY_COUNT, (
        f"expected {EXPECTED_STUDY_COUNT} studies, got {len(studies)}"
    )
    print(f"got {len(studies)} studies")
    print("downloading aptamer_info.xlsx from Box ...")
    download_aptamer_info(aptamer_info_path)
    published = published_aptamers(aptamer_info_path)
    rows = build_manifest_rows(studies, published)
    write_manifest(rows, MANIFEST_PATH)


if __name__ == "__main__":
    main()
