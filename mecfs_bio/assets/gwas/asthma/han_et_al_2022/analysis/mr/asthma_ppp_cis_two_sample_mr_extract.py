from pathlib import PurePath

from mecfs_bio.assets.gwas.asthma.han_et_al_2022.analysis.asthma_ppp_two_sample_mr_cis import (
    HAN_2022_ASTHMA_TSMR,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.copy_file_from_directory_task import (
    CopyFileFromDirectoryTask,
)
from mecfs_bio.build_system.task.two_sample_mr_task import MAIN_RESULT_DF_PATH

HAN_2022_ASTHMA_CIS_PPP_TSMR_EXTRACT = CopyFileFromDirectoryTask.create_result_table(
    asset_id="two_sample_mr_ppp_cis_asthma_han_extract_results",
    source_directory_task=HAN_2022_ASTHMA_TSMR,
    path_inside_directory=PurePath(MAIN_RESULT_DF_PATH),
    extension=".csv",
    read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
)
