from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import \
    concrete_standard_analysis_generator_no_rsid
from mecfs_bio.assets.gwas.human_herpesvirus_7_dna.kamitaki_et_al_2025.raw.kamitaki_et_al_2025_raw import \
    KAMITAKI_ET_AL_HHV7_DNA_RAW
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import GWASLabColumnSpecifiers
from mecfs_bio.build_system.task.pipes.filter_rows_by_info_score import FilterRowsByInfoScorePipe


KAMITAKI_ET_AL_STANARD_ANALYSIS=concrete_standard_analysis_generator_no_rsid(
    base_name="kamitaki_et_al_2025_hhv7",
    raw_gwas_data_task=KAMITAKI_ET_AL_HHV7_DNA_RAW,
    fmt=GWASLabColumnSpecifiers(
        snpid="SNP",
        chrom="CHR",
        pos="BP",
        ea="ALLELE1",
        nea="ALLELE0",
        eaf="A1FREQ",
        p="P_BOLT_LMM_INF",
        beta="BETA",
        se="SE",
        OR=None,
        rsid=None,
        info="INFO",
    ),
    pre_pipe_before_rsid_assignment=FilterRowsByInfoScorePipe(min_score=0.8, info_col="INFO"),
    sample_size=490_401,

)