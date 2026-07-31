"""
The 1000 Genomes 30x EUR reference VCF (hg38), self-hosted and pinned by checksum.

This is gwaslab's EUR panel (multi-allelic variants decomposed, variants normalized,
INFO/AF carrying EUR allele frequency), rehosted rather than fetched from gwaslab's
own URL. gwaslab republishes its panels in place: the copy served at its URL gained
chrX at some point after October 2025, silently changing the variant set without any
version marker. Anything built on top of it would have been irreproducible.

That matters more than usual downstream, because the PPP variant index derived from
this file fixes the row order that every per-protein beta/se file is written
against. Those files carry no variant identifiers, so a source that shifts under the
index does not produce a detectably wrong answer -- it produces silently misaligned
floats across thousands of files.

The checksum below is what pins it: the URL is ours, but only the hash makes the
content immutable.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

THOUSAND_GENOMES_EUR_HG38_30X_VCF = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="thousand_genomes_eur_hg38_30x_vcf",
        group="thousand_genomes",
        sub_group="eur_hg38_30x",
        sub_folder=PurePath("raw"),
        extension=".vcf.gz",
    ),
    url=(
        "https://www.dropbox.com/scl/fi/qtrh2d6u6p8405ug7a6fc/"
        "EUR.ALL.split_norm_af.1kg_30x.hg38.vcf.gz"
        "?rlkey=z25adzkqjurvloak46fl93i8h&dl=1"
    ),
    md5_hash="9e1aa56d8b843fef724b67e23932c95b",
)
