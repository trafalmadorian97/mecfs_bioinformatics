from pathlib import Path
from typing import Sequence

import narwhals
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SE_COL,
)


@frozen
class CaseControlSampleInfo:
    cases: int
    controls: int


@frozen
class GwasSource:
    task: Task
    sample_info: CaseControlSampleInfo
    pipe: DataProcessingPipe = IdentityPipe()


_key_cols = [
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_BETA_COL,
    GWASLAB_SE_COL,
]


@frozen
class FixedEffectsMetaAnalysisTask(Task):
    """
    Task to perform a fixed effects meta analysis on non-overlapping GWAS of the same trait
    Assumes all alleles are expressed with respect to the forward strand
    """

    _meta: Meta
    sources: Sequence[GwasSource]

    def __attrs_post_init__(self):
        assert len(self.sources) > 1

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [item.task for item in self.sources]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        asset = fetch(self.sources[0].task.asset_id)
        df = scan_dataframe_asset(asset, self.sources[0].task.meta)
        df = _select_df_1_columns(df)
        df = df.rename(
            {
                GWASLAB_BETA_COL: GWASLAB_BETA_COL + "_0",
                GWASLAB_SE_COL: GWASLAB_SE_COL + "_0",
            }
        )
        beta_col_list = [GWASLAB_BETA_COL + "_0"]
        se_col_list = [GWASLAB_SE_COL + "_0"]

        i = 0
        for source in self.sources[1:]:
            i += 1
            asset = fetch(source.task.asset_id)
            source_df = scan_dataframe_asset(asset, source.task.meta)
            source_df = source_df.select(
                [
                    GWASLAB_CHROM_COL,
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL,
                    GWASLAB_BETA_COL,
                    GWASLAB_SE_COL,
                ]
            )
            source_df = source_df.rename(
                {
                    GWASLAB_BETA_COL: GWASLAB_BETA_COL + f"_{i}",
                    GWASLAB_SE_COL: GWASLAB_SE_COL + f"_{i}",
                }
            )
            beta_col_list.append(GWASLAB_BETA_COL + f"_{i}")
            se_col_list.append(GWASLAB_SE_COL + f"_{i}")
            df = df.join(
                source_df,
                on=[
                    GWASLAB_CHROM_COL,
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL,
                ],
            )

        meta_beta, meta_std = _fixed_effects_beta_se_cols(
            beta_cols=beta_col_list,
            se_cols=se_col_list,
        )
        df = df.with_columns(
            meta_beta.alias(GWASLAB_BETA_COL),
            meta_std.alias(GWASLAB_SE_COL),
        )
        out_path = scratch_dir / (self.asset_id + ".parquet")
        df.sink_parquet(
            out_path,
        )
        return FileAsset(out_path)


def _fixed_effects_beta_se_cols(
    beta_cols: list[str], se_cols: list[str]
) -> tuple[narwhals.Expr, narwhals.Expr]:
    numerator = narwhals.col(beta_cols[0]) / (narwhals.col(se_cols[0]) ** 2)
    denominator = 1 / (narwhals.col(se_cols[0]) ** 2)
    for beta, std in zip(beta_cols[1:], se_cols[1:]):
        numerator += narwhals.col(beta) / narwhals.col(std) ** 2
        denominator += 1 / narwhals.col(std) ** 2
    return numerator / denominator, (1 / denominator).sqrt()


def _select_df_1_columns(
    df: narwhals.LazyFrame,
) -> narwhals.LazyFrame:
    cols = [
        GWASLAB_CHROM_COL,
        GWASLAB_POS_COL,
        GWASLAB_EFFECT_ALLELE_COL,
        GWASLAB_NON_EFFECT_ALLELE_COL,
        GWASLAB_BETA_COL,
        GWASLAB_SE_COL,
    ]
    schema = df.collect_schema()
    if GWASLAB_RSID_COL in schema:
        cols.append(GWASLAB_RSID_COL)
    return df.select(cols)
