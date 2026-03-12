from pathlib import PurePath

from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec, DataFrameTextFormat
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

KAMITAKI_ET_AL_HHV7_DNA_RAW=DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id="kamitaki_et_al_2025_hhb7_raw",
        trait="hhb7_dna",
        project="kamataki_et_al",
        sub_dir="raw",
        project_path=PurePath("HHV7_invnorm.bgen.MAF.stats.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://data.broadinstitute.org/lohlab/virome_summary_statistics/blood/HHV7_invnorm.bgen.MAF.stats.gz",
    md5_hash=None,
)