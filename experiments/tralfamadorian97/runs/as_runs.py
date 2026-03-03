from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.ankylosing_spondylitis.finngen.processed.finngen_ank_spond_sumstats import \
    FINGNEN_ANK_SPOND_SUMSTATS
from mecfs_bio.assets.gwas.ankylosing_spondylitis.finngen.raw.raw_finngen_ank_spond_data import \
    FINNGEN_ANKYLOSING_SPONDYLITIS_DATA_RAW
from mecfs_bio.assets.gwas.ankylosing_spondylitis.million_veterans.processed.mv_eur_ank_spond_sumstats import \
    MILLION_VETERAN_ANKYLOSING_SPONDYLITIS_SUMSTATS
from mecfs_bio.assets.gwas.ankylosing_spondylitis.million_veterans.processed.mv_eur_ank_spond_sumstats_dump_to_parquet import \
    MILLION_VETERANS_ANK_SPOND_SUMSTATS_37_DUMP_TO_PARQUET
from mecfs_bio.assets.gwas.ankylosing_spondylitis.million_veterans.raw.raw_mv_eur_ank_spond_data import \
    MILLION_VETERAN_ANKYLOSING_SPONDYLITIS_EUR_DATA_RAW
from mecfs_bio.assets.gwas.ankylosing_spondylitis.ukbb.processed.ukbb_ank_spond_sumstats import \
    UK_BIOBANK_ANKYLOSING_SPONDYLITIS_SUMSTATS
from mecfs_bio.assets.gwas.ankylosing_spondylitis.ukbb.processed.ukbb_eur_ank_spond_filtered import \
    FILTERED_UKBB_ANK_SPOND
from mecfs_bio.assets.gwas.ankylosing_spondylitis.ukbb.processed.ukbb_eur_ank_spond_parquet import \
    UKBB_ANK_SPOND_PARQUET
from mecfs_bio.assets.gwas.ankylosing_spondylitis.ukbb.raw.raw_ukbb_eur_ank_spond_data import \
    UK_BIOBANK_ANKYLOSING_SPONDYLITIS_EUR_DATA_RAW


def go():
    DEFAULT_RUNNER.run(
        (
            [
                # MILLION_VETERAN_ANKYLOSING_SPONDYLITIS_SUMSTATS,
                # FINGNEN_ANK_SPOND_SUMSTATS
                # MILLION_VETERAN_ANKYLOSING_SPONDYLITIS_EUR_DATA_RAW,
                # FINNGEN_ANKYLOSING_SPONDYLITIS_DATA_RAW
             # UK_BIOBANK_ANKYLOSING_SPONDYLITIS_EUR_DATA_RAW
             #    FILTERED_UKBB_ANK_SPOND
             #    UKBB_ANK_SPOND_PARQUET
             #    MILLION_VETERANS_ANK_SPOND_SUMSTATS_37_DUMP_TO_PARQUET
             #    FILTERED_UKBB_ANK_SPOND
                UK_BIOBANK_ANKYLOSING_SPONDYLITIS_SUMSTATS
             ]
        ),

        incremental_save=True,
        must_rebuild_transitive=[
        ]

    )

if __name__ == '__main__':
    go()
