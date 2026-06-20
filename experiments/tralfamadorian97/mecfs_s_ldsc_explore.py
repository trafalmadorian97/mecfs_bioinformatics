from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.s_ldsc.decode_me_sldsc import DECODE_ME_S_LDSC
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.s_ldsc.decode_me_sldsc_all_tables_combined import \
    DECODE_ME_SLDSC_ALL_TABLES_COMBINED


def go():
    # task = DECODE_ME_S_LDSC.partitioned_tasks["multi_tissue_chromatin"].add_categories_task_unwrap
    task_2=DECODE_ME_S_LDSC.partitioned_tasks["cahoy_cns"].cell_analysis_task
    DEFAULT_RUNNER.run(
        [
            DECODE_ME_SLDSC_ALL_TABLES_COMBINED
            # task_2
         ],
        must_rebuild_transitive=[DECODE_ME_SLDSC_ALL_TABLES_COMBINED]
    )

if __name__ == '__main__':
    go()