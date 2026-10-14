"""
The 1000 Genomes phase 3 EUR reference VCF (hg19), self-hosted and pinned by checksum.

This is gwaslab's 1kg_eur_hg19 panel: multi-allelic variants are decomposed, variants
are normalized, and INFO/AF carries the EUR allele frequency. Records are indexed for
chromosomes 1-22 and X only.

It is rehosted rather than fetched from gwaslab's own URL, for the same reasons as
eur_hg38_30x_vcf.py: gwaslab republishes its panels in place, and by September 2026
the gwaslab link for this file served a "File Deleted" page. The copy hosted here was
taken from a local gwaslab download whose md5 matched gwaslab's catalogue entry.

The checksum below is what pins the content.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

THOUSAND_GENOMES_EUR_HG19_VCF = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="thousand_genomes_eur_hg19_vcf",
        group="thousand_genomes",
        sub_group="eur_hg19",
        sub_folder=PurePath("raw"),
        extension=".vcf.gz",
    ),
    url=(
        "https://www.dropbox.com/scl/fi/q3wndz2wockuqto5pq0cc/"
        "EUR.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz"
        "?rlkey=zq8fzrjtbbl20r1wr74p4tqa1&dl=1"
    ),
    md5_hash="2c78cb84cb1f90b576510decc45e5b9b",
)
