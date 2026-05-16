from mecfs_bio.assets.gwas.me_cfs.decode_me.raw_gwas_data.decode_me_quality_control_snps import DECODE_ME_QC_SNPS
from mecfs_bio.assets.gwas.me_cfs.decode_me_male.raw.decode_me_male_raw import DECODE_ME_GWAS_MALE_TASK
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec, DataFrameParquetFormat
from mecfs_bio.build_system.task.filter_snps_task import FilterSNPsTask

DECODE_ME_MALE_FILTER_SNPS_TASK = FilterSNPsTask(
    raw_gwas_task=DECODE_ME_GWAS_MALE_TASK,
    snp_list_task=DECODE_ME_QC_SNPS,
    meta=FilteredGWASDataMeta(
        id=AssetId("decode_me_male_filtered_for_quality_control"),
        trait="ME_CFS",
        project="DecodeME",
        sub_dir="processed",
        read_spec=DataFrameReadSpec(format=DataFrameParquetFormat()),
    ),
)
