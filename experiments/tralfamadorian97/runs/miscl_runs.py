from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.mi_ldl_correlation import MI_LDL_CORRELATION
from mecfs_bio.assets.gwas.multi_trait.lcv.ldl_mi_lcv_analysis import LDL_MI_LCV_ANALYSIS, LDL_HARMONIZE_WITH_MI


def run_miscl_analysis():
    DEFAULT_RUNNER.run(
        # [MILLION_VETERAN_MI_EUR_DATA_RAW]+MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.get_terminal_tasks(),
        # [MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.manhattan_task],
        [
            # LDL_HARMONIZE_WITH_MI
            # LDL_MI_LCV_ANALYSIS
            # THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOL
         ]
        +MI_LDL_CORRELATION.terminal_tasks(),
    # [MI_EUR_MANHATTAN],
        incremental_save=True,
        must_rebuild_transitive=[
            # THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOL
        ]
        # must_rebuild_transitive=[MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.manhattan_task]


    )


if __name__ == "__main__":
    run_miscl_analysis()
