from mecfs_bio.asset_generator.mixer_asset_generator import (
    univariate_mixer_asset_generator,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import (
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.auxiliary.sample_info import (
    JOHNSTON_ET_AL_SAMPLE_INFO,
)
from mecfs_bio.assets.reference_data.mixer.processed.mixer_g1000_plink_eur_extracted import (
    MIXER_G1000_PLINK_EXTRACTED,
)
from mecfs_bio.build_system.task.mixer.mixer_task import MixerDataSource

JOHNSTON_ET_AL_UNIVARIATE_MIXER = univariate_mixer_asset_generator(
    base_name="johnston_et_al_pain",
    trait_1_source=MixerDataSource(
        JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
        alias="Multsite_pain",
        sample_info=JOHNSTON_ET_AL_SAMPLE_INFO,
    ),
    reference_data_directory_task=MIXER_G1000_PLINK_EXTRACTED,
    name_in_plot="Multisite_pain",
    threads=5,
    # reps=list(range(1, 13)),
)
