"""
This task uses GWASLAB's wrapper around the original LDSC code to estimate genetic correlation by cross-trait
ld-score regression.

While LDSC is a useful and powerful technique, its results should be taken with a grain of salt.

Here is a discussion (which covers both single and cross-trait analysis) from chapter 23 of the Handbook of Statistical Genomics:
Balding, David J., Ida Moltke, and John Marioni, eds. Handbook of statistical genomics. John Wiley & Sons, 2019.
Chapter authors are S. Burgess, C.N. Foley and V. Zuber




"
*A criticism of LD score regression is that every analysis for each pair of traits uses the same LD scores as the dependent variable in the regression model (and as LD scores have been precomputed by its proponents, literally the same LD scores are used in the majority of applied analyses). This means that any influential points in the regression will affect not only one LD score regression analysis, but all such analyses. LD scores are also likely to be a ‘weak instrument’ in the language of Mendelian randomization, as they will only explain a small proportion of variance in the dependent variable. Additionally, due to the scale of the data, it is not possible to provide a visual representation of an LD score regression analysis. Standard regression diagnostics are rarely, if ever, performed. Finally, results from LD score regression are not always consistent with known causal relationships; for example, the method did not find evidence for a genetic correlation between LDL cholesterol and CHD risk that survived a multiple testing correction (Bulik-Sullivan et al., 2015). The method has utility in mapping the genetic distance between related phenotypes, such as determining how closely related different psychiatric disorders are in terms of their genetic predictors (Cross-Disorder Group of the Psychiatric Genomics Consortium, 2013). However, the reliance of the method on numerous linearity and independence assumptions, incorrect weighting in the linear regression model (correct weights would require computation of the Cholesky decomposition of a matrix with dimension equal to the number of genetic variants in the model – misspecified weights are recommended for use in practice), and lack of validation against known causal relationships mean that results from the method should not be treated too seriously as an assessment of causality.*
"


"""

from pathlib import Path, PurePath
from typing import Sequence

import gwaslab
import pandas as pd
import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_sumstats import read_sumstats
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import GenomeBuild
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_rsid import (
    complement_reverse_expr,
    match_flipped_reference_expr,
    match_reference_expr,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_RSID_COL,
)

logger = structlog.get_logger()


@frozen
class FilterSettings:
    """
    Options for SNP filtering for CT-LDSC
    """

    remove_indels: bool = True
    remove_palindromic: bool = True
    remove_hla: bool = True
    keep_only_hapmap: bool = True


def filter_sumstats(sumstats: gwaslab.Sumstats, settings: FilterSettings, build: str):
    """
    Performing filtering roughly consistent with the procedure described in the methods section of Bulik-Sullivan et al.
    """
    if settings.remove_indels:
        logger.debug("filtering indels")
        sumstats.filter_indel(inplace=True, mode="out")
    if settings.remove_palindromic:
        logger.debug("filtering palindromes")
        sumstats.filter_palindromic(inplace=True, mode="out")
    if settings.remove_hla:
        logger.debug("filtering hla region")
        sumstats.exclude_hla(inplace=True)
    if settings.keep_only_hapmap:
        sumstats.filter_hapmap3(inplace=True, build=build)
    logger.debug("dropping duplicate rsids")
    len_before = len(sumstats.data)
    sumstats.data = sumstats.data.drop_duplicates(subset=[GWASLAB_RSID_COL], keep=False)
    len_after = len(sumstats.data)
    logger.debug(f"dropped {len_before - len_after} variants with identical rsids")


@frozen
class SumstatsSource:
    """
    A source of GWASlab sumstats to use for computing genetic correlation
    """

    task: Task
    alias: str
    pipe: DataProcessingPipe = IdentityPipe()

    @property
    def asset_id(self) -> AssetId:
        return self.task.asset_id


@frozen
class GeneticCorrelationByCTLDSCTask(Task):
    """
    Estimate genetic correlation by cross-trait linkage disequilibrium score regression.
    Note that when CT-LDSC runs, it will also print observed scale heritabilities.
    For binary traits, observed scale heritability is not directly meaningful, and must be converted to liability scale.

    See: Bulik-Sullivan, Brendan, et al. "An atlas of genetic correlations across human diseases and traits." Nature genetics 47.11 (2015): 1236-1241.
    """

    source_sumstats_tasks: Sequence[SumstatsSource]
    ld_ref_task: Task
    ld_file_filename_pattern: str
    _meta: Meta
    build: GenomeBuild
    filter_settings: FilterSettings = FilterSettings()

    def __attrs_post_init__(self):
        assert len(self.source_sumstats_tasks) > 1

    @property
    def _ld_ref_id(self) -> AssetId:
        return self.ld_ref_task.asset_id

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        result = [self.ld_ref_task]
        for item in self.source_sumstats_tasks:
            result.append(item.task)
        return result

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        ref_id = self._ld_ref_id
        ref_asset = fetch(ref_id)
        assert isinstance(ref_asset, DirectoryAsset)
        results = []
        for i in range(len(self.source_sumstats_tasks) - 1):
            i_sumstats, i_name = load_and_preprocess_sumstats(
                source=self.source_sumstats_tasks[i],
                fetch=fetch,
                settings=self.filter_settings,
                build=self.build,
            )
            for j in range(i + 1, len(self.source_sumstats_tasks)):
                j_sumstats, j_name = load_and_preprocess_sumstats(
                    source=self.source_sumstats_tasks[j],
                    fetch=fetch,
                    settings=self.filter_settings,
                    build=self.build,
                )
                compatible_rsids = get_compatible_snps_polars(
                    i_sumstats.data, j_sumstats.data
                )
                j_sumstats.data = j_sumstats.data.loc[
                    j_sumstats.data[GWASLAB_RSID_COL].isin(
                        compatible_rsids[GWASLAB_RSID_COL]
                    )
                ]
                i_sumstats.estimate_rg_by_ldsc(
                    other_traits=[j_sumstats],
                    rg=f"{i_name},{j_name}",
                    ref_ld_chr=str(ref_asset.path) + self.ld_file_filename_pattern,
                    w_ld_chr=str(ref_asset.path) + self.ld_file_filename_pattern,
                    build=self.build,
                )
            results.append(i_sumstats.ldsc_rg)

        out_df: pd.DataFrame = pd.concat(results, ignore_index=True)
        logger.debug(
            f"rg_ldsc_results:\n{out_df}",
        )
        out_path = scratch_dir / "ldsc_rg.csv"
        out_df.to_csv(out_path, index=False)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        sources: Sequence[SumstatsSource],
        ld_ref_task: Task,
        build: GenomeBuild,
        ld_file_filename_pattern: str = "/LDscore.@",
    ):
        assert len(sources) > 1
        sumstats_meta = sources[0].task.meta
        assert isinstance(sumstats_meta, GWASLabSumStatsMeta)
        meta = ResultTableMeta(
            id=asset_id,
            trait="multi_trait",
            project="genetic_correlation",
            sub_dir=PurePath("ct_ldsc"),
            extension=".csv",
            read_spec=DataFrameReadSpec(DataFrameTextFormat(",")),
        )
        return cls(
            meta=meta,
            ld_ref_task=ld_ref_task,
            ld_file_filename_pattern=ld_file_filename_pattern,
            source_sumstats_tasks=sources,
            build=build,
        )


def load_and_preprocess_sumstats(
    source: SumstatsSource, fetch: Fetch, settings: FilterSettings, build: GenomeBuild
) -> tuple[gwaslab.Sumstats, str]:
    name = source.alias
    sumstats_asset = fetch(source.asset_id)
    logger.debug(f"reading sumstats for {name}")
    sumstats = read_sumstats(sumstats_asset)
    assert GWASLAB_RSID_COL in sumstats.data.columns
    sumstats.data = source.pipe.process_pandas(sumstats.data)
    filter_sumstats(sumstats, settings, build=build)
    return sumstats, name


def get_compatible_snps_polars(df_i: pd.DataFrame, df_j: pd.DataFrame) -> pd.DataFrame:
    """
    Get a DataFrame of rsIDs corresponding to SNPs for which both dataframes have compatible alleles.
    This is necessary because certain dbSNP rsIDs actually correspond to multiple allele pairs.

    pandas joins on string columns can very slow, so this operation uses polars.
    """
    pl_i = pl.from_pandas(
        df_i[
            [GWASLAB_RSID_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL]
        ]
    ).with_columns(
        pl.col(GWASLAB_RSID_COL).cast(pl.String()).alias("i_rsid"),
        pl.col(GWASLAB_EFFECT_ALLELE_COL).cast(pl.String()).alias("i_ea"),
        pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).cast(pl.String()).alias("i_nea"),
    )
    pl_j = pl.from_pandas(
        df_j[
            [GWASLAB_RSID_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL]
        ]
    ).with_columns(
        pl.col(GWASLAB_RSID_COL).cast(pl.String()).alias("j_rsid"),
        pl.col(GWASLAB_EFFECT_ALLELE_COL).cast(pl.String()).alias("j_ea"),
        pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).cast(pl.String()).alias("j_nea"),
    )
    pl_j = pl_j.with_columns(
        complement_reverse_expr("j_ea").alias("j_ea_other_strand"),
        complement_reverse_expr("j_nea").alias("j_nea_other_strand"),
    )
    joined = pl_i.join(pl_j, left_on="i_rsid", right_on="j_rsid")
    joined = joined.with_columns(
        (
            match_reference_expr(
                ea_col="i_ea",
                nea_col="i_nea",
                ref_ea_col="j_ea",
                ref_nea_col="j_nea",
            )
            | match_flipped_reference_expr(
                ea_col="i_ea",
                nea_col="i_nea",
                ref_ea_col="j_ea",
                ref_nea_col="j_nea",
            )
            | match_reference_expr(
                ea_col="i_ea",
                nea_col="i_nea",
                ref_ea_col="j_ea_other_strand",
                ref_nea_col="j_nea_other_strand",
            )
            | match_flipped_reference_expr(
                ea_col="i_ea",
                nea_col="i_nea",
                ref_ea_col="j_ea_other_strand",
                ref_nea_col="j_nea_other_strand",
            )
        ).alias("match")
    )
    logger.debug(
        f"Out of {len(joined)} matching rsids, {joined['match'].sum()} had matching alleles"
    )
    result = joined.filter(pl.col("match")).select(GWASLAB_RSID_COL).to_pandas()
    return result
