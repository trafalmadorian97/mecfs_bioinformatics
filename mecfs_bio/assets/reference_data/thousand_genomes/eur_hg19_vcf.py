"""
The 1000 Genomes phase 3 EUR reference VCF (hg19), as distributed by gwaslab.

Multi-allelic variants are decomposed, variants are normalized, and INFO/AF carries
the EUR allele frequency. Records are indexed for chromosomes 1-22 and X only.

gwaslab republishes its panels in place (its hg38 panel silently gained chrX), so the
checksum is what pins the content. If the URL stops serving this checksum, rehost the
file as was done for eur_hg38_30x_vcf.py.
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
        "https://www.dropbox.com/scl/fi/sxnjd37t7e677wlluzeit/"
        "EUR.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz"
        "?rlkey=9om2qu3tfjnq8nxj0tgwp8bd0&dl=1"
    ),
    md5_hash="2c78cb84cb1f90b576510decc45e5b9b",
)
