from pathlib import Path
import gwaslab as gl
import pandas as pd


def go():
    pth = Path("assets/base_asset_store/gwas/ME_CFS/DecodeME/gwaslab_sumstats/decode_me_gwas_1_37_sumstats_rsids_from_annovar.pickle")
    sumstats=gl.load_pickle(str(pth))
    filtered= sumstats.filter_hapmap3()
    dat:pd.DataFrame=filtered.data
    dat=  dat.drop("SNPID", axis=1)
    dat.to_csv(Path("martin_exchange_5")/"hapmap_filtered_sumstats.csv", index=False)


if __name__ =="__main__":
    go()