"""
LAVA vignette data: small sample data bundled with the LAVA R package.

These are used for fast system tests that don't require downloading
the full g1000 EUR LD reference (~300MB+).

The vignette data includes:
- A small g1000_test LD reference in PLINK format (.bed/.bim/.fam)
- A test locus definition file (42 loci)
- Summary statistics for several phenotypes (bmi, neuro, depression, etc.)
- A sample overlap file

Source: https://github.com/josefin-werme/LAVA/tree/main/vignettes/data
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameWhiteSpaceSepTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.build_system.task.download_files_into_directory_task import (
    DownloadEntry,
    DownloadFilesIntoDirectoryTask,
)

_VIGNETTE_RAW_URL = "https://raw.githubusercontent.com/josefin-werme/LAVA/refs/heads/main/vignettes/data"

# LD reference: g1000_test .bed/.bim/.fam (PLINK format)
LAVA_VIGNETTE_LD_REF = DownloadFilesIntoDirectoryTask(
    meta=SimpleDirectoryMeta(id="lava_vignette_g1000_test_ld_ref"),
    entries=[
        DownloadEntry(
            url=f"{_VIGNETTE_RAW_URL}/g1000_test.bed",
            filename="g1000_test.bed",
            md5hash=None,
            post_download_command=None,
        ),
        DownloadEntry(
            url=f"{_VIGNETTE_RAW_URL}/g1000_test.bim",
            filename="g1000_test.bim",
            md5hash=None,
            post_download_command=None,
        ),
        DownloadEntry(
            url=f"{_VIGNETTE_RAW_URL}/g1000_test.fam",
            filename="g1000_test.fam",
            md5hash=None,
            post_download_command=None,
        ),
    ],
)

# Locus definition file (42 test loci)
LAVA_VIGNETTE_LOCI = DownloadFileTask(
    ReferenceFileMeta(
        group="lava_vignette",
        sub_group="data",
        sub_folder=PurePath("raw"),
        id="lava_vignette_test_loci",
        filename="test",
        extension=".loci",
    ),
    url=f"{_VIGNETTE_RAW_URL}/test.loci",
    md5_hash=None,
)

# BMI summary statistics (continuous phenotype)
# Columns: CHR BP SNPID_UKB A1 A2 MAF BETA SE STAT P NMISS INFO_UKB
LAVA_VIGNETTE_BMI_SUMSTATS = DownloadFileTask(
    ReferenceFileMeta(
        group="lava_vignette",
        sub_group="data",
        sub_folder=PurePath("raw"),
        id="lava_vignette_bmi_sumstats",
        filename="bmi.sumstats",
        extension=".txt",
        read_spec=DataFrameReadSpec(DataFrameWhiteSpaceSepTextFormat(comment_code="#")),
    ),
    url=f"{_VIGNETTE_RAW_URL}/bmi.sumstats.txt",
    md5_hash=None,
)

# Neuroticism summary statistics (continuous phenotype)
# Columns: CHR POS SNP A1 A2 EAF BETA SE Z P N_analyzed INFO
LAVA_VIGNETTE_NEURO_SUMSTATS = DownloadFileTask(
    ReferenceFileMeta(
        group="lava_vignette",
        sub_group="data",
        sub_folder=PurePath("raw"),
        id="lava_vignette_neuro_sumstats",
        filename="neuro.sumstats",
        extension=".txt",
        read_spec=DataFrameReadSpec(DataFrameWhiteSpaceSepTextFormat(comment_code="#")),
    ),
    url=f"{_VIGNETTE_RAW_URL}/neuro.sumstats.txt",
    md5_hash=None,
)

# Depression summary statistics (binary phenotype: 170756 cases, 329443 controls)
# Columns: SNP A1 A2 Freq LogOR StdErrLogOR P N OR
LAVA_VIGNETTE_DEPRESSION_SUMSTATS = DownloadFileTask(
    ReferenceFileMeta(
        group="lava_vignette",
        sub_group="data",
        sub_folder=PurePath("raw"),
        id="lava_vignette_depression_sumstats",
        filename="depression.sumstats",
        extension=".txt",
        read_spec=DataFrameReadSpec(DataFrameWhiteSpaceSepTextFormat(comment_code="#")),
    ),
    url=f"{_VIGNETTE_RAW_URL}/depression.sumstats.txt",
    md5_hash=None,
)
