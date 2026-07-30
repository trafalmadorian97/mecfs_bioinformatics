"""
Task to download gene-level summary data from the Astra-Zenica rare variant collapsing analysis of the UK biobank trait

120010#Ever had chronic Fatigue Syndrome or Myalgic Encephalomyelitis (M.E.)

I believe this trait reflects responses to the UK Biobank Pain Survey.  If so, it is a higher quality way to define a
UK Biobank ME/CFS cohort than some of the alternatives.

See: Samms, Gemma L., and Chris P. Ponting. "Defining a High-Quality Myalgic Encephalomyelitis/Chronic Fatigue Syndrome cohort in UK Biobank." NIHR open research 5 (2025): 39.


The AZ Phewas homepage is  https://azphewas.com/

The methodology is described in

Wang, Quanli, et al. "Rare variant contribution to human disease in 281,104 UK Biobank exomes."
Nature 597.7877 (2021): 527-532.


I have re-uploaded the csv table from azphewas.com to my dropbox



"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MECFS_AZ_PHEWAS = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id="mecfs_az_phewas_gene_level_summary",
        trait="ME_CFS",
        project="az_phewas",
        sub_dir="raw",
        project_path=PurePath(
            "120010_Ever_had_chronic_Fatigue_Syndrome_or_Myalgic_Encephalomy_GL_0513T221814Z_AZ_PheWAS_Portal.csv"
        ),
        read_spec=DataFrameReadSpec(
            DataFrameTextFormat(separator=","),
        ),
    ),
    url="https://www.dropbox.com/scl/fi/bio5hw7q07u1cnb4mjy17/120010_Ever_had_chronic_Fatigue_Syndrome_or_Myalgic_Encephalomy_GL_0513T221814Z_AZ_PheWAS_Portal.csv?rlkey=90lwlfl7tz37eugz7tel3rbcb&dl=1",
    md5_hash="8726146c22f0a8efbe10a3e24a3e0e4c",
)
