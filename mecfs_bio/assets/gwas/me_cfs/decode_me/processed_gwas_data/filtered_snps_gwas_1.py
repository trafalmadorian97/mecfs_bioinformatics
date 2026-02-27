"""
Task to filter genetic variants from DecodeME GWAS 1 to keep only those passing quality control.

This is consistent with the README here: https://osf.io/rgqs3/files/axp4k.
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.raw_gwas_data.decode_me_gwas_1 import (
    DECODE_ME_GWAS_1_TASK,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.raw_gwas_data.decode_me_quality_control_snps import (
    DECODE_ME_QC_SNPS,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.task.filter_snps_task import FilterSNPsTask

DECODE_ME_FILTER_SNPS_GWAS_1_TASK = FilterSNPsTask(
    raw_gwas_task=DECODE_ME_GWAS_1_TASK,
    snp_list_task=DECODE_ME_QC_SNPS,
    meta=FilteredGWASDataMeta(
        id=AssetId("decode_me_gwas_1_filtered_for_quality_control"),
        trait="ME_CFS",
        project="DecodeME",
        sub_dir="processed",
        read_spec=DataFrameReadSpec(format=DataFrameParquetFormat()),
    ),
)
