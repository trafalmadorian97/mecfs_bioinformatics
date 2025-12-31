from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import \
    concrete_standard_analysis_generator_no_rsid
from mecfs_bio.assets.gwas.imaging_derived_heart_phenotypes.pirruccello_et_al_2022.processed.extracted_ra_fac import \
    PIRRUCCELLO_EXTRACTED_RA_FAC_DATA
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import GWASLabColumnSpecifiers

RA_FAC_STANDARD_ANALYSIS_ASSIGN_RSID = concrete_standard_analysis_generator_no_rsid(
    base_name="pirruccello_et_al_2022_ra_fac",
    raw_gwas_data_task=PIRRUCCELLO_EXTRACTED_RA_FAC_DATA,
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
    sample_size=41_135,
    include_master_gene_lists=False
)
