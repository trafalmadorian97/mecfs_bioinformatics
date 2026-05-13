"""
As mentioned in the multivariate GWAS tutorial here: https://github.com/GenomicSEM/GenomicSEM/wiki/5.-Multivariate-GWAS
Genomic SEM requires thousand genome reference data.
This Task downloads this data, re-hosted on Dropbox.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

GENOMES1K_REFERENCE_FOR_GENOMIC_SEM = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="genomic_sem_reference",
        sub_group="genomes_1k",
        sub_folder=PurePath("raw"),
        id=AssetId("1k_genomes_for_genomic_sem"),
        extension=".gz",
        filename="reference.1000G.maf.0.005.txt",
    ),
    url="https://www.dropbox.com/scl/fi/pp3iftty01ekhilk19gmx/reference.1000G.maf.0.005.txt.gz?rlkey=kgxww2tf7118j582jx46dq0rd&dl=1",
    md5_hash="b710949ecf4b7ff56c395a12174261b6",
)
