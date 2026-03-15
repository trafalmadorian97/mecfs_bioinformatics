"""
Asset generator to model the genetic architecture of ME/CFS using univariate MiXeR.
"""

from mecfs_bio.asset_generator.mixer_asset_generator import (
    univariate_mixer_asset_generator,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.auxiliary.prevalance_info import (
    DECODE_ME_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.assets.reference_data.mixer.processed.mixer_g1000_plink_eur_extracted import (
    MIXER_G1000_PLINK_EXTRACTED,
)
from mecfs_bio.build_system.task.mixer.mixer_task import MixerDataSource
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL

DECODE_ME_UNIVARIATE_MIXER = univariate_mixer_asset_generator(
    base_name="decode_me",
    trait_1_source=MixerDataSource(
        DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.join_task,
        alias="DecodeME",
        pipe=RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
        sample_info=DECODE_ME_PREVALENCE_INFO,
    ),
    reference_data_directory_task=MIXER_G1000_PLINK_EXTRACTED,
    name_in_plot="DecodeME",
    threads=6,
    # reps=list(range(1, 13)),
)
