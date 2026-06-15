from pathlib import PurePath

from mecfs_bio.assets.gwas.multi_trait.genomic_sem.decode_me_minus_pain_subtraction_ols import (
    DECODE_ME_MINUS_PAIN_SUBTRACTION_OLS,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.task.copy_file_from_directory_task import (
    CopyFileFromDirectoryTask,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_config import (
    GWAS_RESULTS_SUBDIR,
    SUBTRACTION_R_FILENAME,
)

DECODE_ME_MINUS_PAIN_OLS_EXTRACT = (
    CopyFileFromDirectoryTask.create_from_gwas_data_in_directory(
        asset_id="decode_me_minus_pain_ols_extracted",
        source_directory_task=DECODE_ME_MINUS_PAIN_SUBTRACTION_OLS,
        path_inside_directory=PurePath(GWAS_RESULTS_SUBDIR) / SUBTRACTION_R_FILENAME,
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        trait="me_cfs_minus_pain_ols",
        project="decode_me_minus_johnston",
        sub_dir="raw",
        extension=".parquet",
    )
)
