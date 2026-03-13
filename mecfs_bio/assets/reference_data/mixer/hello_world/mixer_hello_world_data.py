"""
MiXeR hello-world reference data: small sample data for fast system tests.

Downloads the mixer_hello_world.tar.gz tarball from the comorment/mixer repository,
extracts it, and generates .ld files from the included PLINK data (chr21 and chr22).

The prepared directory contains everything needed to run the cross-trait MiXeR
hello-world example: .bed/.bim/.fam, .ld files, trait1.sumstats.gz, trait2.sumstats.gz,
annotation files, and gene ontology files.

Source: https://github.com/comorment/mixer/tree/main/reference/mixer_hello_world
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import ReferenceFileMeta
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.build_system.task.extract_tar_gzip_task import ExtractTarGzipTask
from mecfs_bio.build_system.task.mixer.mixer_task import MixerLDGenerationTask

_HELLO_WORLD_TARBALL_URL = (
    "https://raw.githubusercontent.com/comorment/mixer/"
    "refs/heads/main/reference/mixer_hello_world/mixer_hello_world.tar.gz"
)

# Download the hello-world tarball
MIXER_HELLO_WORLD_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="mixer",
        sub_group="hello_world",
        sub_folder=PurePath("raw"),
        extension=".gz",
        id="mixer_hello_world_raw",
    ),
    url=_HELLO_WORLD_TARBALL_URL,
    md5_hash=None,
)

# Extract the tarball
MIXER_HELLO_WORLD_EXTRACTED = ExtractTarGzipTask.create(
    asset_id="mixer_hello_world_extracted",
    source_task=MIXER_HELLO_WORLD_RAW,
)

# Prepare reference data: copy all files and generate .ld files from PLINK data
MIXER_HELLO_WORLD_PREPARED = MixerLDGenerationTask(
    meta=SimpleDirectoryMeta(id="mixer_hello_world_prepared"),
    plink_data_task=MIXER_HELLO_WORLD_EXTRACTED,
    chromosomes=(21, 22),
    bfile_prefix_pattern="g1000_eur_hm3_chr{chr}",
)
