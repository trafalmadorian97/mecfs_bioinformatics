"""
Categories for tissue and cell types from the GTEx and Franke lab datasets

source:
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

FICUANE_2018_FRANKE_GTEX_CATEGORIES = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="ficuane_2018_combined_gtex_franke_categories",
        group="linkage_disequilibrium_score_regression_aux_data",
        sub_group="data_from_papers",
        sub_folder=PurePath("ficuane_2018"),
        filename="ficuane_2018_combined_gtex_franke_categories.csv",
        extension=".csv",
        read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
    ),
    url="https://www.dropbox.com/scl/fi/2cr5ipddzkde6rj33tcxk/Finucane_2018_combined_gtex_franke.csv?rlkey=w5knnradpjugl2l1orml9j0g9&dl=1",
    md5_hash="3f121ae8936e380828a14dc76c2fd309",
)
