from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.analysis_results.bentham_2015_standard_analysis import \
    BENTHAM_LUPUS_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.raw_gwas_data.bentham_2015_raw_build_37 import \
    BENTHAM_2015_RAW_BUILD_37


def lupus_analysis():
    """
    Script to analyze the lupus GWAS of Bentham et al
    """
    DEFAULT_RUNNER.run(
        BENTHAM_LUPUS_STANDARD_ANALYSIS.get_terminal_tasks(),
       # [BENTHAM_2015_RAW_BUILD_37] ,
        incremental_save=True,
        must_rebuild_transitive=[],
    )


if __name__ == "__main__":
    lupus_analysis()
