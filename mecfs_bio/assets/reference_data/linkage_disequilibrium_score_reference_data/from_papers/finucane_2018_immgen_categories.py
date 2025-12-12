"""
Categories for mouse immune cell types from the immgen project

source:
Suplementary Table 10 in
Supplementary material to
Finucane, Hilary K., et al. "Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types." Nature genetics 50.4 (2018): 621-629.

"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

FICUANE_2018_IMMGEN_CATEGORIES = DownloadFileTask(
    meta=ReferenceFileMeta(
        asset_id="ficuane_2018_immgen_categories",
        group="linkage_disequilibrium_score_regression_aux_data",
        sub_group="data_from_papers",
        sub_folder=PurePath("ficuane_2018"),
        filename="ficuane_2018_immgen_categories.csv",
        extension=".csv",
        read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
    ),
    url="https://www.dropbox.com/scl/fi/587e0eypohr2azdlt0x54/Finucane_2018_Immgen_categories.csv?rlkey=nmjm7nrrpinqpm78q9wrbfft0&dl=1",
    md5_hash="6a1d3c65ace87be0be6f56bd493eb3db",
)
