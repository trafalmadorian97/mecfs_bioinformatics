"""
Categories for tissue and cell types from the roadmap project used in Finucane et al.

source:
My judgement

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

ROADMAP_CELL_TYPE_CATEGORIES_FOR_LDSC = DownloadFileTask(
    meta=ReferenceFileMeta(
        asset_id="roadmap_tissue_type_categories",
        group="linkage_disequilibrium_score_regression_aux_data",
        sub_group="custom_data",
        sub_folder=PurePath("roadmap_categories"),
        filename="roadmap_cell_type_categories",
        extension=".csv",
        read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
    ),
    url="https://www.dropbox.com/scl/fi/cz1yty9sj8vsqnddarpqb/cell_type_series_categorized.csv?rlkey=upn71pljutjywj7dz3fyb785e&st=azq3dty6&dl=1",
    md5_hash="f1f871fa446ef19ef6473860b9f12741",
)
