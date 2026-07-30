from mecfs_bio.asset_generator.ukbb_ppp_rg_asset_generator import generate_ppp_rg_assets
from mecfs_bio.assets.gwas.asthma.han_et_al_2022.analysis.asthma_standard_analysis import (
    HAN_ASTHMA_STANDARD_ANALYSIS,
)
from mecfs_bio.build_system.task.dataframe_output import ParquetWriteOptions
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_rg_task import PppRgConfig

ASTHMA_PARQUET = GwasLabSumstatsToTableTask.create_from_source_task(
    source_tsk=HAN_ASTHMA_STANDARD_ANALYSIS.labeled_lead_variant_tasks.raw_sumstats_task,
    asset_id="asthma_keep_version_parquet_table_from_sumstats",
    sub_dir="processed",
    pipe=CompositePipe([ComputeBetaPipe(), ComputeSEPipe()]),
    write_options=ParquetWriteOptions(
        compression="zstd", byte_stream_split_floats=True
    ),
)


ASTHMA_PPP_RG_CIS_EXCLUDED = generate_ppp_rg_assets(
    "han_asthma_ppp_rg_cis_excluded",
    trait_task=ASTHMA_PARQUET,
    config=PppRgConfig(
        variant_set="cis_excluded", trait_total_sample_size=393859
    ),  # sample size is from sumstats file
)
