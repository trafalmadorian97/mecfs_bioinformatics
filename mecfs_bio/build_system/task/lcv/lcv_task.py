
"""
Task to apply the Latent Causal Variable technique of O'Connor and Price to attempt
to estimate the causal director between two genetically correlated traits.
The key output value is GCP: Genetic Causality Proportion.

Citation:
O’Connor, Luke J., and Alkes L. Price. "Distinguishing genetic correlation
from causation across 52 diseases and complex traits." Nature genetics 50.12 (2018): 1728-1734.

"""
from pathlib import Path, PurePath

import narwhals
from attrs import frozen
import narwhals as nw
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec, DataFrameParquetFormat
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.consolidate_ld_scores_task import LD_SCORE_CHROM_COL, LD_SCORE_POS_COL, \
    LD_SCORE_RSID_COL, LD_SCORE_LD_SCORE_COL
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import MULTI_TRAIT
from mecfs_bio.build_system.task.lcv.lcv_core import run_lcv
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import GWASLAB_BETA_COL, GWASLAB_SE_COL, GWASLAB_CHROM_COL, GWASLAB_POS_COL, \
    GWASLAB_RSID_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL

Z_SCORE_COL = "_z_score_"
Z_SCORE_1 = "_z_score_1_"
Z_SCORE_2 = "_z_score_2_"

@frozen
class LCVConfig:
    chisq_exclude_factor_threshold:float


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
    config: LCVConfig

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
        aligned = align_traits_and_ld(
            df_trait_1=df_trait_1,
            df_trait_2=df_trait_2,
            ld_scores=df_ld_scores,
        )
        lcv_result = run_lcv(
            ld_scores=aligned[LD_SCORE_LD_SCORE_COL].to_numpy(),
            z1=aligned[Z_SCORE_1].to_numpy(),
            z2=aligned[Z_SCORE_2].to_numpy(),
            chisq_exclude_factor_threshold=self.config.chisq_exclude_factor_threshold,
        )
        result_df = lcv_result.to_df()
        out_path = scratch_dir/"result.parquet"
        result_df.write_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(cls,
               asset_id:str,
    trait_1_data: Task,
    trait_2_data: Task,
    consolidated_ld_scores: Task,
    config: LCVConfig,
    ):
        meta= ResultTableMeta(
            id=AssetId(asset_id),
            trait=MULTI_TRAIT,
            project="causal_analysis",
            sub_dir=PurePath("lcv"),
            read_spec=DataFrameReadSpec(
                DataFrameParquetFormat()
            ),
            extension=".parquet"

        )
        return cls(
            trait_1_data=trait_1_data,
            trait_2_data=trait_2_data,
            consolidated_ld_scores=consolidated_ld_scores,
            config=config,
            meta=meta,
        )




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
            LD_SCORE_CHROM_COL, LD_SCORE_POS_COL,LD_SCORE_RSID_COL
        ]
    ).collect().sort(
        by=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, GWASLAB_RSID_COL],
    )
    return joined
