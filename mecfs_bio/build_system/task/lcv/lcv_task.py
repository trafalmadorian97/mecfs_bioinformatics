
"""
Task to apply the Latent Causal Variable technique of O'Connor and Price to attempt
to estimate the causal director between two genetically correlated traits.
The key output value is GCP: Genetic Causality Proportion.

Citation:
O’Connor, Luke J., and Alkes L. Price. "Distinguishing genetic correlation
from causation across 52 diseases and complex traits." Nature genetics 50.12 (2018): 1728-1734.

"""
from pathlib import Path

import narwhals
from attrs import frozen
import narwhals as nw
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import GWASLAB_BETA_COL, GWASLAB_SE_COL, GWASLAB_CHROM_COL, GWASLAB_POS_COL, \
    GWASLAB_RSID_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL

Z_SCORE_COL = "_z_score_"
Z_SCORE_1 = "_z_score_1_"
Z_SCORE_2 = "_z_score_2_"

@frozen
class LCVTask(Task):
    """
    Task to apply the Latent Causal Variable technique of O'Connor and Price to attempt
    to estimate the causal director between two genetically correlated traits.
    The key output value is GCP: Genetic Causality Proportion.

    Assume that the trat 1 and trait 2 datasets have already been harmonized via HarmonizeGWASWithReferenceViaRSIDTask
    """
    _meta: Meta
    trait_1_data: Task
    trait_2_data: Task
    consolidated_ld_scores: Task

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.trait_1_data, self.trait_2_data]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        trait_1_asset = fetch(self.trait_1_data.asset_id)
        trait_2_asset = fetch(self.trait_2_data.asset_id)
        ld_scores_asset = fetch(self.consolidated_ld_scores.asset_id)
        df_trait_1 = make_z_score_frame( scan_dataframe_asset(asset=trait_1_asset,meta=self.trait_1_data.meta))
        df_trait_2 = make_z_score_frame(scan_dataframe_asset(asset=trait_2_asset,meta=self.trait_2_data.meta))
        df_ld_scores = scan_dataframe_asset(asset=ld_scores_asset,meta=self.consolidated_ld_scores.meta)






def make_z_score_frame(
        df: narwhals.LazyFrame,
)->nw.LazyFrame:
    return df.with_columns(
        (nw.col(GWASLAB_BETA_COL)/ nw.col(GWASLAB_SE_COL)).alias(Z_SCORE_COL),
    ).select(
     GWASLAB_CHROM_COL,
     GWASLAB_POS_COL,
        GWASLAB_RSID_COL,
        Z_SCORE_COL
    )

def align_traits_and_ld(
        df_trait_1: narwhals.LazyFrame,
        df_trait_2: narwhals.LazyFrame,
        ld_scores: narwhals.LazyFrame,
) -> nw.DataFrame:
    df_trait_1=df_trait_1.rename({Z_SCORE_COL:Z_SCORE_1 })
    df_trait_2=df_trait_2.rename({Z_SCORE_COL:Z_SCORE_2 })
    joined = df_trait_1.join(df_trait_2, on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, GWASLAB_RSID_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL], how="inner").join(
        ld_scores, left_on=[
            GWASLAB_CHROM_COL,
            GWASLAB_POS_COL,
            GWASLAB_RSID_COL,
        ],
        right_on=[
            "CHR", "BP","SNP"
        ]
    ).collect()
    return joined
