import polars as pl
from pathlib import Path

from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, \
    GWASLAB_NON_EFFECT_ALLELE_COL

finn_path = Path(
    "assets/base_asset_store/gwas/ankylosing_spondylitis/finngne/processed/finngen_ank_spond_harmonized_dump_to_parquet.parquet")

mv_path = Path("assets/base_asset_store/gwas/ankylosing_spondylitis/million_veterans/processed/million_veterans_eur_ank_spond_harmonized_dump_to_parquet.parquet")
def go():
    df = pl.read_parquet(finn_path)
    counts = df.group_by(
                [
                    GWASLAB_CHROM_COL,
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL
                ],
    ).len().sort(by="len",descending=True)
    import pdb; pdb.set_trace()
    print("yo")

def go_mv():
    df = pl.read_parquet(mv_path)
    counts = df.group_by(
        [
            GWASLAB_CHROM_COL,
            GWASLAB_POS_COL,
            GWASLAB_EFFECT_ALLELE_COL,
            GWASLAB_NON_EFFECT_ALLELE_COL
        ],
    ).len().sort(by="len",descending=True)
    import pdb; pdb.set_trace()
    print("yo")

if __name__ == "__main__":
    go_mv()


