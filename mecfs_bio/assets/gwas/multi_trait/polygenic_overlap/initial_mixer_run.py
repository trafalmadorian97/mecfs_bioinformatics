from mecfs_bio.assets.gwas.me_cfs.decode_me.auxiliary.prevalance_info import DECODE_ME_PREVALENCE_INFO
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import \
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import \
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.auxiliary.sample_info import JONHSTON_ET_AL_SAMPLE_INFO
from mecfs_bio.assets.reference_data.mixer.processed.mixer_g1000_plink_eur_extracted import MIXER_G1000_PLINK_EXTRACTED
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import QuantPhenotype
from mecfs_bio.build_system.task.mixer_bivariate_task import MixerTask, MixerDataSource
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL

MULTISITE_PAIN_DECODE_ME_INITIAL_MIXER=MixerTask.create(
    asset_id="initial_mixer_run",
    trait_1_source=MixerDataSource(

        DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.join_task,
        alias="DecodeME",
        pipe=RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
        sample_info=DECODE_ME_PREVALENCE_INFO
    ),
    trait_2_source=MixerDataSource(
        JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
        alias="Multisite_pain",
        sample_info=JONHSTON_ET_AL_SAMPLE_INFO
    ),
    ref_data_directory_task=MIXER_G1000_PLINK_EXTRACTED,
    # rep_file_pattern=None,
    num_reps=1,
    chr_args="1-4"
)