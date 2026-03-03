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
    GWASLAB_SE_COL, GWASLAB_EFFECTIVE_SAMPLE_SIZE,
)
import structlog
logger = structlog.get_logger()


@frozen
class CaseControlSampleInfo:
    cases: int
    controls: int
    def effective_sample_size(self)-> int:
        return int(4/( 1/self.cases + 1/self.controls  ))

SampleInfo = CaseControlSampleInfo # add other types of sample info later

@frozen
class GwasSource:
    task: Task
    sample_info: SampleInfo
    pipe: DataProcessingPipe = IdentityPipe()



@frozen
class FixedEffectsMetaAnalysisTask(Task):
    """
    Task to perform a fixed effects meta analysis on non-overlapping GWAS of the same trait
    Assumes all alleles are expressed with respect to the forward strand

    The variants present in the output dataframe will be equal to the intersection of the variants in the studies
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
        df = self.sources[0].pipe.process(scan_dataframe_asset(asset, self.sources[0].task.meta))
        _check_unique_variants(df)
        _check_nonzero_se(df)
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
            _check_unique_variants(source_df)
            _check_nonzero_se(source_df)
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
            df_forward_match = df.join(
                source_df,
                on=[
                    GWASLAB_CHROM_COL,
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL,
                ],
            ).with_columns(
                narwhals.lit(False).alias(f"flipped_{i}")
            )
            source_df_flipped = get_reversed(source_df, beta_col=GWASLAB_BETA_COL + f"_{i}")
            df_reverse_match = df.join(
                source_df_flipped,
                on=[
                    GWASLAB_CHROM_COL,
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL,
                ],
            ).with_columns(
                narwhals.lit(True).alias(f"flipped_{i}")
            )
            df =narwhals.concat([df_forward_match, df_reverse_match], how="vertical")

        meta_beta, meta_std = _fixed_effects_beta_se_cols(
            beta_cols=beta_col_list,
            se_cols=se_col_list,
        )
        df = df.with_columns(
            meta_beta.alias(GWASLAB_BETA_COL),
            meta_std.alias(GWASLAB_SE_COL),
        )
        out_path = scratch_dir / (self.asset_id + ".parquet")
        df = add_effective_sample_size_column(
            df, [item.sample_info for item in self.sources]
        )
        df.sink_parquet(
            out_path,
        )
        final  =narwhals.scan_parquet(out_path, backend="polars")
        report_flips(
           final,num_sources=len(self.sources),
        )
        report_ouput_size(final)
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


def add_effective_sample_size_column(
        out_df: narwhals.LazyFrame,
        sample_info: list[SampleInfo],
) -> narwhals.LazyFrame:
    effective_sample_size = sum(
        item.effective_sample_size() for item in sample_info
    )
    return out_df.with_columns(
        narwhals.lit(effective_sample_size).alias(GWASLAB_EFFECTIVE_SAMPLE_SIZE)
    )


def get_reversed(
    df: narwhals.LazyFrame,
        beta_col:str
    )   -> narwhals.LazyFrame:
    return df.with_columns(
        narwhals.col(GWASLAB_EFFECT_ALLELE_COL).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
        narwhals.col(GWASLAB_NON_EFFECT_ALLELE_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
        (-1*narwhals.col(beta_col)).alias(beta_col),
    )


def report_flips(
        df: narwhals.LazyFrame,
        num_sources: int
):
    sum_cols = [
        narwhals.col(f"flipped_{i}").sum().alias(f"num_flipped_{i}") for i in range(1, num_sources )
    ]
    result = df.select(
        *sum_cols,
    ).collect()
    logger.debug(f"Flipped alleles:\n{result}")

def report_ouput_size(
        df: narwhals.LazyFrame,
):
    l = df.select(narwhals.len()).collect().item()
    logger.debug(f"Final meta-analysis has {l} variants")



def _check_unique_variants(
        df: narwhals.LazyFrame
):
    assert  df.select(
        [GWASLAB_CHROM_COL,
         GWASLAB_POS_COL,
         GWASLAB_EFFECT_ALLELE_COL,
         GWASLAB_NON_EFFECT_ALLELE_COL,]
    ).unique().select(narwhals.len()).collect().item() == df.select(narwhals.len()).collect().item()


def _check_nonzero_se(
        df: narwhals.LazyFrame,
):
    assert df.select(
        (narwhals.col(GWASLAB_SE_COL)>0).all()
    ).collect().item()