from mecfs_bio.asset_generator.concrete_standard_analysis_generator import \
    concrete_standard_analysis_generator_assumme_already_has_rsid, concrete_standard_analysis_generator_no_rsid
from mecfs_bio.assets.gwas.imaging_derived_heart_phenotypes.pirruccello_et_al_2022.processed.extracted_rvef_data import \
    PIRRUCCELLO_EXTRACTED_RVEF_DATA
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import GWASLabColumnSpecifiers

RVEF_STANDARD_ANALYSIS = concrete_standard_analysis_generator_assumme_already_has_rsid(
    base_name="pirruccello_et_al_2022_rvef",
    raw_gwas_data_task=PIRRUCCELLO_EXTRACTED_RVEF_DATA,
    fmt=GWASLabColumnSpecifiers(
        rsid="SNP",
        chrom="CHR",
        pos="BP",
        ea="ALLELE1",
        nea="ALLELE0",
        eaf="A1FREQ",
        p="P_LINREG",
        beta="BETA",
        se="SE",
        OR=None,
        snpid=None,
        info="INFO"

    ),
    sample_size=41_135
)

RVEF_STANDARD_ANALYSIS_ASSIGN_RSID = concrete_standard_analysis_generator_no_rsid(
    base_name="pirruccello_et_al_2022_rvef",
    raw_gwas_data_task=PIRRUCCELLO_EXTRACTED_RVEF_DATA,
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
        info="INFO"

    ),
    sample_size=41_135
)
