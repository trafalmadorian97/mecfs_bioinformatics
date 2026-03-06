import gwaslab




import polars as pl


def go():
    dat= pl.read_csv("assets/base_asset_store/gwas/ankylosing_spondylitis/million_veterans/raw/GCST90476232.h.tsv.gz",separator="\t")
    import pdb; pdb.set_trace()
    print("yo")

def check_harm_effect(
                      ):
    pre = gwaslab.load_pickle("assets/base_asset_store/gwas/ankylosing_spondylitis/million_veterans/gwaslab_sumstats/million_veterans_ank_spond_sumstats.pickle")
    post = gwaslab.load_pickle("assets/base_asset_store/gwas/ankylosing_spondylitis/million_veterans/gwaslab_sumstats/million_veterans_ank_spond_sumstats_37_harmonized.pickle")
    pre_filtered= pre.data.loc[pre.data["rsID"]=="rs34686676"]
    post_filtered = post.data.loc[post.data["rsID"]=="rs34686676"]
    import pdb; pdb.set_trace()
    print("yo")


if __name__ == "__main__":
    check_harm_effect()