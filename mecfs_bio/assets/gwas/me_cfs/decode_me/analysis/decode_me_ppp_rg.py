from mecfs_bio.asset_generator.ukbb_ppp_rg_asset_generator import generate_ppp_rg_assets
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.build_38.decode_me_gwas_1_assign_rsids_build_38 import (
    DECODE_ME_GWAS_1_ASSIGN_RSIDS_BUILD_38,
)
from mecfs_bio.build_system.task.dataframe_output import ParquetWriteOptions
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_rg_task import PppRgConfig

DECODE_ME_ORIG_PARQUET = GwasLabSumstatsToTableTask.create_from_source_task(
    source_tsk=DECODE_ME_GWAS_1_ASSIGN_RSIDS_BUILD_38,  # DECODE_ME_GWAS_1_SUMSTATS_MINIMAL_FILTERING,
    asset_id="DecodeME_keep_version_parquet_table_from_sumstats",
    sub_dir="processed",
    # pipe=CompositePipe([ComputeBetaPipe(), ComputeSEPipe()]),
    write_options=ParquetWriteOptions(
        compression="zstd", byte_stream_split_floats=True
    ),
)


DECODE_ME_PPP_RG_CIS_EXCLUDED = generate_ppp_rg_assets(
    "decode_me_ppp_rg_cis_excluded",
    trait_task=DECODE_ME_ORIG_PARQUET,
    config=PppRgConfig(variant_set="cis_excluded"),
)

DECODE_ME_PPP_RG_CIS_INCLUDED = generate_ppp_rg_assets(
    "decode_me_ppp_rg_cis_included",
    trait_task=DECODE_ME_ORIG_PARQUET,
    config=PppRgConfig(variant_set="all_variants"),
)
