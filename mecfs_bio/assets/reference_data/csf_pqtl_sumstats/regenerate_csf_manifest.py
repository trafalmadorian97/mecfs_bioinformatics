"""
Regenerate the Western et al. 2024 CSF pQTL summary-statistics manifest.

The manifest is a static list mapping every published CSF aptamer to the GWAS
Catalog study whose GWAS-SSF file holds its summary statistics. It is the CSF
analogue of regenerate_ppp_manifest.py.

Glossary (the terms this script maps between):

  aptamer      One SomaScan affinity reagent: a modified-DNA binder that measures one
               protein target. One row of the aptamer table. 7,584 exist; 7,008 are
               published (the rest have a non-null "Step Removed").
  analyte      An aptamer's stable identifier string, e.g. "X13681.173" (the
               "Analytes" column). The aptamer primary key -- gene symbol is not
               unique across aptamers.
  seq_id       The same aptamer identifier in dash form, e.g. "13681-173" (the
               "SeqId" column).
  target name  The human-readable protein name an aptamer targets, e.g.
               "Beta-crystallin B2" (the "TargetFullName" column). NOT unique:
               several aptamers can share one target name.
  aptamer      aptamer_info.xlsx, the SomaScan 7k table with one row per aptamer.
    table      From the paper's Box deposit.
  study        One GWAS Catalog deposit: the genome-wide association results for a
               single aptamer measured in CSF. There are 7,008, one per published
               aptamer.
  accession    A study's stable GWAS Catalog identifier, e.g. "GCST90421540".
  trait        The human-readable label the GWAS Catalog gives a study, e.g.
               "Beta-crystallin B2 levels". Unique per study; the Catalog
               disambiguates aptamers that share a target name by appending
               "(analyte X####.##)".

The script maps each study (identified by its accession) to the aptamer it measured
(identified by its analyte), using the trait as the only key the two sources share.

Two upstream sources are combined:

  1. The GWAS Catalog REST API, listing the 7,008 studies of PMID 39528825, each
     with its accession and trait.
  2. The aptamer table (aptamer_info.xlsx), which supplies each aptamer's analyte,
     seq_id, target name and publication status: 7,584 rows, of which the 576 with
     a non-null "Step Removed" are the unpublished aptamers, leaving exactly 7,008.

The trait is the only link between the two. The Catalog assigns each study a UNIQUE
trait, having itself disambiguated shared target names as noted above. A four-rule
resolver recovers a complete accession -> analyte bijection; the resolver asserts the
bijection, so a future upstream re-curation that breaks the mapping fails loudly here
rather than silently at build time.

Run via:
  pixi r python mecfs_bio/assets/reference_data/csf_pqtl_sumstats/regenerate_csf_manifest.py

Downloads no summary statistics; it fetches only study metadata (a few MB of JSON),
the 826 KB aptamer table, and each study's tiny md5sum.txt (7,008 small requests, run
in parallel) to record the .tsv.gz md5 the build later verifies. Re-run to pick up any
upstream re-curation, then diff the result against the committed manifest.
"""

import functools
import json
import re
import tempfile
import urllib.request
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.error import URLError

import pandas as pd

from mecfs_bio.build_system.task.csf_database.gwas_catalog_url import (
    gwas_catalog_bucket,
)
from mecfs_bio.constants.csf_database_constants import GcstAccession
from mecfs_bio.util.download.verify import verify_hash
from mecfs_bio.util.retry.call_with_retries import call_with_retries

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
# MD5 of aptamer_info.xlsx (826,640 bytes) as served by the endpoint above. Pinning it
# turns a silent upstream swap into a loud failure: a changed table would otherwise be
# accepted and could shift the manifest. Re-run to update the hash if Box re-deposits.
_APTAMER_INFO_MD5 = "4913a42d671f3e35d4997eab9d3e670f"

# aptamer_info.xlsx columns.
_ANALYTES_COL = "Analytes"  # aptamer primary key, e.g. X13681.173
_SEQID_COL = "SeqId"  # e.g. 13681-173
_UNIPROT_COL = "UniProt"
_ENTREZ_SYMBOL_COL = "EntrezGeneSymbol"
_TARGET_FULL_NAME_COL = "TargetFullName"
_STEP_REMOVED_COL = "Step Removed"  # non-null => aptamer not published

# Rule 1: the Catalog disambiguates shared target names by appending this to the trait.
_ANALYTE_IN_TRAIT_PATTERN = re.compile(r"\(analyte (X[\d.]+)\)")

# Rule 4: three Casein kinase II traits the first three rules cannot resolve, because
# the Catalog writes the two catalytic subunits as "alpha-1" / "alpha-2" while the
# aptamer table writes them "alpha" / "alpha'" -- so the trait carries no
# "(analyte ...)" tiebreak (rule 1 misses) and does not match any target name +
# " levels" (rules 2-3 miss).
#
# The pairing is not a free choice: the subunit identity is fixed by UniProt --
# alpha-1 = CSNK2A1 (P68400), alpha-2 = CSNK2A2 (P19784) -- which pins each trait to
# exactly one aptamer (the alpha-2 monomer is the sole published CSNK2A2 monomer; the
# two heterotetramers differ only in CSNK2A1 vs CSNK2A2). To re-verify a row, open its
# study page and read the trait, then confirm the analyte's subunit in the aptamer
# table against UniProt:
#   https://www.ebi.ac.uk/gwas/studies/GCST90422297  "...subunit alpha-2 levels"
#   https://www.ebi.ac.uk/gwas/studies/GCST90426283  "...alpha-1: beta heterotetramer levels"
#   https://www.ebi.ac.uk/gwas/studies/GCST90426284  "...alpha-2: beta heterotetramer levels"
#   https://www.uniprot.org/uniprotkb/P68400  CSNK2A1 (alpha-1)
#   https://www.uniprot.org/uniprotkb/P19784  CSNK2A2 (alpha-2)
_TRAIT_OVERRIDES = {
    "Casein kinase II subunit alpha-2 levels": "X13681.173",  # CSNK2A2 monomer
    "Casein kinase II alpha-1: beta heterotetramer levels": "X5225.50",  # CSNK2A1|CSNK2B
    "Casein kinase II alpha-2: beta heterotetramer levels": "X5226.36",  # CSNK2A2|CSNK2B
}

# GWAS Catalog publishes a per-study md5sum.txt next to the .tsv.gz, with one line
# "<md5> <filename>" per file. We pin the md5 of the exact .tsv.gz each aptamer task
# downloads, so a corrupted or re-deposited file fails the build's download check.
_FTP_ROOT = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
# Bounded concurrency for the 7,008 tiny md5sum.txt fetches; polite to the EBI FTP,
# which refuses connections under a heavier burst.
_MD5_FETCH_WORKERS = 8

# Output manifest columns.
ANALYTE_MANIFEST_COLUMN = "analyte"
SEQ_ID_MANIFEST_COLUMN = "seq_id"
UNIPROT_MANIFEST_COLUMN = "uniprot"
ENTREZ_SYMBOL_MANIFEST_COLUMN = "entrez_gene_symbol"
TARGET_FULL_NAME_MANIFEST_COLUMN = "target_full_name"
ACCESSION_MANIFEST_COLUMN = "accession"
MD5_MANIFEST_COLUMN = "md5"

MANIFEST_COLUMNS = [
    ANALYTE_MANIFEST_COLUMN,
    SEQ_ID_MANIFEST_COLUMN,
    UNIPROT_MANIFEST_COLUMN,
    ENTREZ_SYMBOL_MANIFEST_COLUMN,
    TARGET_FULL_NAME_MANIFEST_COLUMN,
    ACCESSION_MANIFEST_COLUMN,
    MD5_MANIFEST_COLUMN,
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
    """Download aptamer_info.xlsx from the Box deposit and verify its MD5."""
    urllib.request.urlretrieve(_BOX_DOWNLOAD_URL, dest)
    verify_hash(dest, _APTAMER_INFO_MD5)
    return dest


def _read_md5sum_file(url: str) -> str:
    with urllib.request.urlopen(url) as response:
        return response.read().decode()


def _fetch_sumstats_md5(accession: str) -> str:
    """The published md5 of a study's {accession}.tsv.gz, read from its md5sum.txt.

    The FTP layout mirrors gwas_catalog_sumstats_url: the md5sum.txt sits beside the
    .tsv.gz in the study's bucket directory, with one "<md5> <filename>" line per file.
    The EBI FTP intermittently refuses connections under load, so the read is retried.
    """
    bucket = gwas_catalog_bucket(GcstAccession(accession))
    url = f"{_FTP_ROOT}/{bucket}/{accession}/md5sum.txt"
    filename = f"{accession}.tsv.gz"
    text = call_with_retries(
        functools.partial(_read_md5sum_file, url), retry_on=(URLError,)
    )
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[1] == filename:
            return parts[0]
    raise AssertionError(f"{filename} not listed in {url}")


def fetch_sumstats_md5s(accessions: list[str]) -> dict[str, str]:
    """Fetch the .tsv.gz md5 for every accession, in parallel over the tiny requests."""
    with ThreadPoolExecutor(max_workers=_MD5_FETCH_WORKERS) as pool:
        md5s = list(pool.map(_fetch_sumstats_md5, accessions))
    return dict(zip(accessions, md5s))


def published_aptamers(aptamer_info_path: Path) -> pd.DataFrame:
    """The 7,008 published aptamers (Step Removed is null), indexed by analyte."""
    frame = pd.read_excel(aptamer_info_path, sheet_name="Sheet1")
    published = frame[frame[_STEP_REMOVED_COL].isna()].copy()
    assert len(published) == EXPECTED_STUDY_COUNT, (
        f"expected {EXPECTED_STUDY_COUNT} published aptamers, got {len(published)}"
    )
    return published.set_index(_ANALYTES_COL, drop=False)


def build_resolver(
    published: pd.DataFrame,
) -> Callable[[str], tuple[str, str]]:
    """Return a resolver mapping a trait to (analyte, rule), closed over the
    published aptamers. The second element names the rule that matched, for reporting.

    Four rules. Rule 4 (the hardcoded override) is a short-circuit checked FIRST; the
    other three are then tried in order 1 -> 2 -> 3:
      1. an explicit "(analyte X####.##)" in the trait (the Catalog's own tiebreak);
      2. an exact "<target name> levels" match to a uniquely-named aptamer;
      3. the same, case-insensitively;
      4. a hardcoded override for three Casein kinase II traits (see _TRAIT_OVERRIDES).
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
    studies: list[tuple[str, str]],
    published: pd.DataFrame,
    md5_by_accession: dict[str, str],
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
        # 15 published aptamers have no EntrezGeneSymbol in the table; fall back to
        # their UniProt id (all 15 have one) so the symbol -- which feeds the slim
        # file's asset_id and project -- stays meaningful and non-empty.
        symbol = info[_ENTREZ_SYMBOL_COL]
        if pd.isna(symbol):
            symbol = info[_UNIPROT_COL]
        rows.append(
            {
                ANALYTE_MANIFEST_COLUMN: analyte,
                SEQ_ID_MANIFEST_COLUMN: info[_SEQID_COL],
                UNIPROT_MANIFEST_COLUMN: info[_UNIPROT_COL],
                ENTREZ_SYMBOL_MANIFEST_COLUMN: symbol,
                TARGET_FULL_NAME_MANIFEST_COLUMN: info[_TARGET_FULL_NAME_COL],
                ACCESSION_MANIFEST_COLUMN: accession,
                MD5_MANIFEST_COLUMN: md5_by_accession[accession],
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
    rows_sorted = sorted(rows, key=lambda row: row[ANALYTE_MANIFEST_COLUMN])
    frame = pd.DataFrame.from_records(rows_sorted, columns=MANIFEST_COLUMNS)
    # Every cell must be populated -- in particular the UniProt fallback in
    # build_manifest_rows must have filled every missing gene symbol, so no NaN
    # reaches the CSV to be serialized as a literal "nan" or an empty asset_id.
    nan_cells = frame.isna().sum()
    assert not frame.isna().to_numpy().any(), f"manifest has NaN cells:\n{nan_cells}"
    frame.to_csv(manifest_path, index=False)
    print(f"wrote {len(rows_sorted)} rows to {manifest_path}")


def main() -> None:
    print("fetching study list from GWAS Catalog ...")
    studies = fetch_studies()
    assert len(studies) == EXPECTED_STUDY_COUNT, (
        f"expected {EXPECTED_STUDY_COUNT} studies, got {len(studies)}"
    )
    print(f"got {len(studies)} studies")
    print("fetching per-study sumstats md5s from the GWAS Catalog FTP ...")
    md5_by_accession = fetch_sumstats_md5s([accession for accession, _ in studies])
    print("downloading aptamer_info.xlsx from Box ...")
    # The xlsx is a transient input, not a committed asset: download it to a temp dir
    # so a re-run leaves no untracked binary in the tracked source tree. Only the CSV
    # manifest is durable.
    with tempfile.TemporaryDirectory() as tmp:
        aptamer_info_path = download_aptamer_info(Path(tmp) / "aptamer_info.xlsx")
        published = published_aptamers(aptamer_info_path)
        rows = build_manifest_rows(studies, published, md5_by_accession)
    write_manifest(rows, MANIFEST_PATH)


if __name__ == "__main__":
    main()
