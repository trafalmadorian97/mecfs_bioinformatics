import polars as pl
from mecfs_bio.assets.gwas.imaging_derived_heart_phenotypes.pirruccello_et_al_2022.raw.raw_right_heart_data import \
    PIRRUCCELLO_RAW_RIGHT_HEART_DATA
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec, DataFrameTextFormat
from mecfs_bio.build_system.task.extraction_one_file_from_zip_task import ExtractFromZipTask

PIRRUCCELLO_EXTRACTED_RVEF_DATA=ExtractFromZipTask.create_from_zipped_gwas_data(
PIRRUCCELLO_RAW_RIGHT_HEART_DATA,
    asset_id="pirruccello_et_al_2022_extracted_rvef_data",
    file_to_extract="invnorm_RVEF.tsv.gz",
    sub_dir="extracted",
    read_spec=DataFrameReadSpec(
        DataFrameTextFormat("\t",
                            schema_overrides={
                                "GENPOS":pl.String()
                            })
    )
)