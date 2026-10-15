from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.ankylosing_spondylitis.multi_study.fixed_effects_meta_analysis.analysis.ank_spond_meta_analysis_manhattan import \
    ANK_SPOND_META_ANALYSIS_MANHATTAN
from mecfs_bio.assets.gwas.ankylosing_spondylitis.multi_study.fixed_effects_meta_analysis.analysis.standard_analysis import \
    ANK_SPOND_FIXED_EFFECTS_STANDARD_ANALYSIS


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
             #    UK_BIOBANK_ANKYLOSING_SPONDYLITIS_SUMSTATS

                # FINGNEN_ANK_SPOND_SUMSTATS_HARMONIZED,
                # MILLION_VETERANS_ANK_SPOND_SUMSTATS_HARMONIZED,
                #
                # UKBB_ANK_SPOND_SUMSTATS_HARMONIZED
                # ANK_SPOND_FIXED_EFFECTS_META_ANALYSIS
                # FINGEN_ANK_SPOND_MANHATTAN
                # MV_ANK_SPOND_MANHATTAN
                ANK_SPOND_META_ANALYSIS_MANHATTAN
            ]+ANK_SPOND_FIXED_EFFECTS_STANDARD_ANALYSIS.get_terminal_tasks()
        ),

        incremental_save=True,
        must_rebuild_transitive=[
        ]

    )

if __name__ == '__main__':
    go()
