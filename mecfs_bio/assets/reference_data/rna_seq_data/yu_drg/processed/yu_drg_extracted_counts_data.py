from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.raw.yu_drg_raw_counts_rdata import YU_DRG_SRC1_RDATA
from mecfs_bio.build_system.task.extract_dataframe_from_rdata_task import ExtractDataFrameFromRDataTask

YU_DRG_EXTRACTED_COUNTS=ExtractDataFrameFromRDataTask.create(
    asset_id="yu_drg_extracted_counts",
    rdata_file_task=YU_DRG_SRC1_RDATA,
    r_dataframe_name="HS.counts.raw",
    r_package_list=['Seurat','base']
)