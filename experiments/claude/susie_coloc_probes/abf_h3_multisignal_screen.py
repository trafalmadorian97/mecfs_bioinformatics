"""
Screens the v1 coloc.abf results for high-H3 pairs that may be explained by multiple causal
signals rather than truly distinct signals.

For each tested pair with PP_H3 > 0.9, looks up the trait's association at the protein's lead
variant position (hg38) in the trait's build-38 harmonized summary statistics. If the
protein's lead is itself genome-wide significant for the trait (P < 5e-8), but is not the
trait's lead, the single-causal-variant assumption may have forced H3 where a secondary trait
signal is shared with the protein. This is only a screen: the protein lead could instead just
tag the trait lead through LD.
"""

from pathlib import Path

import polars as pl

from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_build_38_harmonized import (
    LIU_ET_AL_2023_IBD_EUR_BUILD_38_HARMONIZED,
)
from mecfs_bio.assets.gwas.ldl.willer_et_al.processed.willer_ldl_build_38_harmonized import (
    WILLER_LDL_BUILD_38_HARMONIZED,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_build_38_harmonized import (
    DECODE_ME_GWAS_1_BUILD_38_HARMONIZED,
)
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.processed.ra_seropositive_build_38_harmonized import (
    SEROPOS_RA_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.metadata_to_path.simple_meta_to_path import (
    simple_meta_to_relative_path,
)
from mecfs_bio.build_system.task.base_task import Task

ASSET_ROOT = Path("assets/base_asset_store")
H3_THRESHOLD = 0.9
GWS_MLOG10P = -pl.lit(5e-8).log10()

_TRAITS: dict[str, Task] = {
    "willer_ldl_ppp_coloc": WILLER_LDL_BUILD_38_HARMONIZED,
    "decode_me_gwas_1_ppp_coloc": DECODE_ME_GWAS_1_BUILD_38_HARMONIZED,
    "decode_ra_seropositive_ppp_coloc": SEROPOS_RA_BUILD_38_HARMONIZED,
    "liu_et_al_2023_ibd_eur_ppp_coloc": LIU_ET_AL_2023_IBD_EUR_BUILD_38_HARMONIZED,
}


def _trait_mlog10p_by_pos(trait_task: Task) -> pl.DataFrame:
    path = ASSET_ROOT / simple_meta_to_relative_path(trait_task.meta)
    lf = scan_dataframe_asset(FileAsset(path), trait_task.meta).to_native()
    cols = lf.collect_schema().names()
    mlog = pl.col("MLOG10P") if "MLOG10P" in cols else -pl.col("P").log10()
    return (
        lf.select(
            pl.col("CHR").cast(pl.Utf8),
            pl.col("POS").cast(pl.Int64),
            mlog.alias("TRAIT_MLOG10P_AT_PROTEIN_LEAD"),
        )
        .group_by("CHR", "POS")
        .agg(pl.col("TRAIT_MLOG10P_AT_PROTEIN_LEAD").max())
        .collect()
    )


def main() -> None:
    pl.Config.set_tbl_rows(60)
    pl.Config.set_tbl_width_chars(220)
    for coloc_name, trait_task in _TRAITS.items():
        results_path = next(ASSET_ROOT.glob(f"gwas/*/*/analysis/{coloc_name}/raw_results.parquet"))
        high_h3 = pl.read_parquet(results_path).filter(pl.col("PP_H3") > H3_THRESHOLD)
        trait = _trait_mlog10p_by_pos(trait_task)
        joined = high_h3.with_columns(
            pl.col("CHR").cast(pl.Utf8), pl.col("PROTEIN_LEAD_POS").cast(pl.Int64)
        ).join(trait, left_on=["CHR", "PROTEIN_LEAD_POS"], right_on=["CHR", "POS"], how="left")
        flagged = joined.filter(
            (pl.col("TRAIT_MLOG10P_AT_PROTEIN_LEAD") > GWS_MLOG10P)
            & (pl.col("PROTEIN_LEAD_POS") != pl.col("LEAD_POS"))
        ).with_columns((pl.col("PROTEIN_LEAD_POS") - pl.col("LEAD_POS")).abs().alias("LEAD_DIST_BP"))
        print(f"\n{coloc_name}: {high_h3.height} pairs with PP_H3 > {H3_THRESHOLD}; "
              f"{joined['TRAIT_MLOG10P_AT_PROTEIN_LEAD'].is_null().sum()} protein leads absent from trait; "
              f"{flagged.height} flagged")
        print(
            flagged.select(
                "LOCUS_ID", "GENE", "CIS_STATUS", "LEAD_MLOG10P", "PROTEIN_LEAD_MLOG10P",
                "TRAIT_MLOG10P_AT_PROTEIN_LEAD", "LEAD_DIST_BP", "PP_H3", "PP_H4",
            ).sort("TRAIT_MLOG10P_AT_PROTEIN_LEAD", descending=True)
        )


if __name__ == "__main__":
    main()
