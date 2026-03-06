from mecfs_bio.assets.gwas.ankylosing_spondylitis.ukbb.processed.ukbb_eur_ank_spond_filtered import (
    FILTERED_UKBB_ANK_SPOND,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
    GWASLabCreateSumstatsTask,
)

UK_BIOBANK_ANKYLOSING_SPONDYLITIS_SUMSTATS = GWASLabCreateSumstatsTask(
    df_source_task=FILTERED_UKBB_ANK_SPOND,
    asset_id="uk_biobank_ank_spond_sumstats",
    fmt=GWASLabColumnSpecifiers(
        chrom="chromosome",
        pos="base_pair_location",
        ea="effect_allele",
        nea="other_allele",
        beta="beta",
        se="standard_error",
        eaf="effect_allele_frequency",
        p="p_value",
        rsid="rsid",
    ),
    liftover_to="19",
    genome_build="infer",
    basic_check=True,
)
