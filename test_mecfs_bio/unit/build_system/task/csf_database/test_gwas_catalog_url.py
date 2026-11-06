import pytest

from mecfs_bio.build_system.task.csf_database.gwas_catalog_url import (
    gwas_catalog_bucket,
    gwas_catalog_sumstats_url,
)
from mecfs_bio.constants.csf_database_constants import GcstAccession


def test_gwas_catalog_bucket_boundaries():
    # Buckets are 1,000 wide and run ...001 through ...000 of the next thousand,
    # inclusive: an exact multiple of 1000 is the TOP of its bucket, and +1 opens the
    # next one. These three cases pin both edges.
    assert (
        gwas_catalog_bucket(GcstAccession("GCST90421000"))
        == "GCST90420001-GCST90421000"
    )
    assert (
        gwas_catalog_bucket(GcstAccession("GCST90421001"))
        == "GCST90421001-GCST90422000"
    )
    assert (
        gwas_catalog_bucket(GcstAccession("GCST90422000"))
        == "GCST90421001-GCST90422000"
    )


def test_gwas_catalog_sumstats_url():
    assert gwas_catalog_sumstats_url(GcstAccession("GCST90421540")) == (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90421001-GCST90422000/GCST90421540/GCST90421540.tsv.gz"
    )


def test_gwas_catalog_bucket_rejects_non_accession():
    with pytest.raises(AssertionError):
        gwas_catalog_bucket(GcstAccession("not-an-accession"))
