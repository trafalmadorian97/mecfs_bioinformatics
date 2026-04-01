from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.raw_gwas_data_harmonized.bentham_2015_harmonized_build_38 import (
    BENTHAM_2015_HARMONIZED_BUILD_38,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
)

"""
Noted some missing data here
For instance the variant rs887369 has a p value but no beta.
"""
BENTHAM_2015_SUMSTATS_MINIMAL_PROCESSING = GWASLabCreateSumstatsTask(
    df_source_task=BENTHAM_2015_HARMONIZED_BUILD_38,
    target_asset_id=AssetId(
        "bentham_2015_harmonized_build_38_sumstats_minimal_processing"
    ),
    basic_check=True,
    genome_build="infer",
    fmt="gwascatalog",
    drop_col_list=["effect_allele_frequency"],  # this is empty
)
