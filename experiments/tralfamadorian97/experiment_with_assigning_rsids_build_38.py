from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.build_38.build_38_s_ldsc_export.build_38_s_ldsc_all_tables_bonferroni import \
    BUILD_38_DECODE_ME_S_LDSC_ALL_TABLES_BONFERRONI_ZIP
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.build_38.decode_me_build_38_s_ldsc import DECODE_ME_S_LDSC_BUILD_38
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.build_38.decode_me_gwas_1_assign_rsids_build_38 import \
    DECODE_ME_GWAS_1_ASSIGN_RSIDS_BUILD_38
from mecfs_bio.build_system.scheduler.topological_scheduler import TopologicalSchedulerSettings


def go():
    DEFAULT_RUNNER.run(
        [
            # DECODE_ME_GWAS_1_ASSIGN_RSIDS_BUILD_38,
            BUILD_38_DECODE_ME_S_LDSC_ALL_TABLES_BONFERRONI_ZIP
        ]
        # +
        # DECODE_ME_S_LDSC_BUILD_38.get_terminal_tasks()
        ,
        incremental_save=True,
        settings=TopologicalSchedulerSettings(print_progress=True)
    )


if __name__ == '__main__':
    go()