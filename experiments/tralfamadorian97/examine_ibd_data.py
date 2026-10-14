import polars as pl



def go():
    df = pl.scan_csv("assets/base_asset_store/gwas/inflammatory_bowl_disease/liu_et_al_2023_ibd_meta/raw/sumstats/ibd_EAS_EUR_SiKJEF_meta_IBD.TBL.txt.gz",
                     separator="\t",
                      )
    head = df.head().collect()
    print(head)
    import pdb; pdb.set_trace()
    print("yo")


if __name__ == "__main__":
    go()
