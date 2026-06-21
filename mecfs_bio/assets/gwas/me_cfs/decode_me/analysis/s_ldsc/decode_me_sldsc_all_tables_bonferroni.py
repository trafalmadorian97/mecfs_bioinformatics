from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.s_ldsc.decode_me_sldsc import DECODE_ME_S_LDSC
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.s_ldsc.decode_me_sldsc_all_tables import DECODE_ME_SLDSC_ALL_TABLES
from mecfs_bio.build_system.task.multiple_testing_table_task import MultipleTestingTableTask

DECODE_ME_S_LDSC_ALL_TABLES_BONFERRONI = [

    MultipleTestingTableTask.create_from_result_table_task(
        alpha=0.05,
        method="bonferroni",
        asset_id=tsk.asset_id+"_bonferroni",
        p_value_column="Coefficient_P_value",
        source_task=tsk,
        apply_filter=False
    )
    for tsk in DECODE_ME_SLDSC_ALL_TABLES

]