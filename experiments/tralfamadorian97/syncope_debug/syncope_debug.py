import gwaslab

from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL, GWASLAB_CHROM_COL, GWASLAB_P_COL, GWASLAB_POS_COL, \
    GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL


def go():
    # keaton_pre_liftover = gwaslab.load_pickle("assets/base_asset_store/gwas/diastolic_blood_pressure/keaton_et_al_dbp/gwaslab_sumstats/keaton_dbp_raw_sumstats_task_no_liftover.pickle")
    # syncope_pre_liftover = gwaslab.load_pickle()
    syncope:gwaslab.Sumstats = gwaslab.load_pickle("assets/base_asset_store/gwas/syncope/aegisdottir/gwaslab_sumstats/aegisdottir_et_al_raw_sumstats_exploded_and_scaled.pickle")
    dup = syncope.data.duplicated(subset=[GWASLAB_CHROM_COL,GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL]).any()
    print("yo")
    syncope_hapmap=syncope.filter_hapmap3()
    dup_hap = syncope_hapmap.data.duplicated(subset=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL]).any()
    import pdb; pdb.set_trace()
    print("yo")


if __name__ == "__main__":
    go()