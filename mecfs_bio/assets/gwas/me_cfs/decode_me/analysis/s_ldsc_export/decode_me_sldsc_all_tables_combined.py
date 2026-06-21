from wcmatch.pathlib import PurePath

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.s_ldsc_export.decode_me_sldsc_all_tables import (
    DECODE_ME_SLDSC_ALL_TABLES,
)
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.task.copy_files_into_directory_task import (
    CopyFilesIntoDirectoryTask,
    CopySource,
)

DECODE_ME_SLDSC_ALL_TABLES_COMBINED = CopyFilesIntoDirectoryTask(
    sources=[CopySource(item, suffix=".csv") for item in DECODE_ME_SLDSC_ALL_TABLES],
    meta=ResultDirectoryMeta(
        id="decode_me_sldsc",
        trait="ME_CFS",
        project="DecodeME",
        sub_dir=PurePath("analysis/pooled_s_ldsc_results"),
    ),
)
