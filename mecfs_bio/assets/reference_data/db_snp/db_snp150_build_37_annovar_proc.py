from pathlib import PurePath

from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

"""
See: https://annovar.openbioinformatics.org/en/latest/user-guide/download/#additional-databases
"""
DB_SNP150_ANNOVAR_PROC = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="db_snp_reference_data",
        sub_group="build_37",
        sub_folder=PurePath("annovar"),
        asset_id="db_snp150_annovar_proc",
        filename="hg19_avsnp150",
        extension=".txt",
        read_spec=DataFrameReadSpec(
            DataFrameTextFormat(
                separator=r"\t",
                # comment_char="#",
                has_header=False,
                column_names=[
                    "CHROM",
                    "POS",
                    "END",
                    "REF",
                    "ALT",
                    "rsid",
                ],
            )
        ),
    ),
    url="https://www.dropbox.com/scl/fi/rvonq6jk3o88mtzc83mrn/hg19_avsnp150.txt?rlkey=gitowj4jrw2wjzbx2uyqi2xxp&dl=1",
    md5_hash="c17e8f96e9b36041455069be9c459555",
)
