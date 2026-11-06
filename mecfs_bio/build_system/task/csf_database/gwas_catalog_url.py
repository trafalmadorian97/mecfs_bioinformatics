"""
Build the EBI GWAS Catalog FTP download URL for a summary-statistics accession.

The Catalog groups per-study files into buckets 1,000 accessions wide, named
"GCST{lower}-GCST{upper}" where the range is inclusive and starts at ...001. For
example GCST90421540 lives under GCST90421001-GCST90422000. The bucket is derived
arithmetically from the accession's numeric part rather than hardcoded, so the whole
GCST90421001-GCST90428xxx span of the Western deposit is covered by one rule.
"""

import re

from mecfs_bio.constants.csf_database_constants import GcstAccession

_FTP_ROOT = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
_ACCESSION_PATTERN = re.compile(r"^GCST(\d+)$")
_BUCKET_WIDTH = 1000


def _accession_number(accession: GcstAccession) -> int:
    match = _ACCESSION_PATTERN.match(accession)
    assert match is not None, f"not a GCST accession: {accession!r}"
    return int(match.group(1))


def gwas_catalog_bucket(accession: GcstAccession) -> str:
    """The "GCST{lower}-GCST{upper}" directory a study's files live under.

    Buckets are 1,000 wide and start at ...001, so GCST90421001-GCST90422000 holds
    accessions ...001 through ...000 of the next thousand inclusive.
    """
    number = _accession_number(accession)
    # ...001..000 maps to the same bucket; shift by one so the exact multiple of 1000
    # falls in the lower bucket, matching the Catalog's ...001-...000 ranges.
    lower = ((number - 1) // _BUCKET_WIDTH) * _BUCKET_WIDTH + 1
    upper = lower + _BUCKET_WIDTH - 1
    width = len(str(number))
    return f"GCST{lower:0{width}d}-GCST{upper:0{width}d}"


def gwas_catalog_sumstats_url(accession: GcstAccession) -> str:
    """The full https URL of a study's harmonized-into-SSF .tsv.gz summary statistics."""
    bucket = gwas_catalog_bucket(accession)
    return f"{_FTP_ROOT}/{bucket}/{accession}/{accession}.tsv.gz"
