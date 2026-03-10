import pandas as pd

from pathlib import Path


def go():
    locpath = Path("assets/base_asset_store/reference_data/lava_locus_file/default/raw/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
    num_rows_per_group=10
    df= pd.read_csv(locpath, sep=r"\s+")
    out_frames =[]
    for num,grp in df.groupby("CHR"):
        chrom_rows = len(grp)
        for i in range(0,  chrom_rows+num_rows_per_group, num_rows_per_group    ):
            to_fuse = grp.iloc[i:i+num_rows_per_group]
            if len(to_fuse)>0:
                out_frames.append(
                   pd.DataFrame(
                       {
                           "CHR":[num],
                           "START":[to_fuse["START"].min()],
                           "STOP": [to_fuse["STOP"].max()]
                       }
                   )
                )
    result =pd.concat(out_frames, ignore_index=True)
    import pdb; pdb.set_trace()
    print("yo")

if __name__ == "__main__":
    go()



