from pathlib import PurePath

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.build_38.build_38_s_ldsc_export.build_38_s_ldsc_all_tables import (
    BUILD_38_DECODE_ME_SLDSC_ALL_TABLES,
)
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.task.copy_files_into_directory_task import (
    CopyFilesIntoDirectoryTask,
    CopySource,
)
from mecfs_bio.build_system.task.multiple_testing_table_task import (
    MultipleTestingTableTask,
)
from mecfs_bio.build_system.task.zip_dir_task import ZipDirTask

BUILD_38_DECODE_ME_S_LDSC_ALL_TABLES_BONFERRONI = [
    MultipleTestingTableTask.create_from_result_table_task(
        alpha=0.05,
        method="bonferroni",
        asset_id=tsk.asset_id + "_bonferroni",
        p_value_column="Coefficient_P_value",
        source_task=tsk,
        apply_filter=False,
    )
    for tsk in BUILD_38_DECODE_ME_SLDSC_ALL_TABLES
]

BUILD_38_DECODE_ME_S_LDSC_ALL_TABLES_BONFERRONI_COMBINED = CopyFilesIntoDirectoryTask(
    sources=[
        CopySource(item, suffix=".csv")
        for item in BUILD_38_DECODE_ME_S_LDSC_ALL_TABLES_BONFERRONI
    ],
    meta=ResultDirectoryMeta(
        id="decode_me_sldsc_bonferroni_build_38",
        trait="ME_CFS",
        project="DecodeME",
        sub_dir=PurePath("analysis/pooled_s_ldsc_results"),
    ),
)
BUILD_38_DECODE_ME_S_LDSC_ALL_TABLES_BONFERRONI_ZIP = ZipDirTask.create(
    source_dir_task=BUILD_38_DECODE_ME_S_LDSC_ALL_TABLES_BONFERRONI_COMBINED,
    asset_id="decode_me_sldsc_bonferroni_archive_build_38",
)
