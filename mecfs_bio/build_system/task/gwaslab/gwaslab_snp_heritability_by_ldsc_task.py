"""
This task uses GWASLAB's wrapper around the original LDSC code to estimate compute heritability
by ld-score regression.

While LDSC is a useful and powerful technique, its results should be taken with a grain of salt.

Here is a discussion (which covers both single and cross-trait analysis) from chapter 23 of the Handbook of Statistical Genomics:
Balding, David J., Ida Moltke, and John Marioni, eds. Handbook of statistical genomics. John Wiley & Sons, 2019.
Chapter authors are S. Burgess, C.N. Foley and V. Zuber

"
*A criticism of LD score regression is that every analysis for each pair of traits uses the same LD scores as the dependent variable in the regression model (and as LD scores have been precomputed by its proponents, literally the same LD scores are used in the majority of applied analyses). This means that any influential points in the regression will affect not only one LD score regression analysis, but all such analyses. LD scores are also likely to be a ‘weak instrument’ in the language of Mendelian randomization, as they will only explain a small proportion of variance in the dependent variable. Additionally, due to the scale of the data, it is not possible to provide a visual representation of an LD score regression analysis. Standard regression diagnostics are rarely, if ever, performed. Finally, results from LD score regression are not always consistent with known causal relationships; for example, the method did not find evidence for a genetic correlation between LDL cholesterol and CHD risk that survived a multiple testing correction (Bulik-Sullivan et al., 2015). The method has utility in mapping the genetic distance between related phenotypes, such as determining how closely related different psychiatric disorders are in terms of their genetic predictors (Cross-Disorder Group of the Psychiatric Genomics Consortium, 2013). However, the reliance of the method on numerous linearity and independence assumptions, incorrect weighting in the linear regression model (correct weights would require computation of the Cholesky decomposition of a matrix with dimension equal to the number of genetic variants in the model – misspecified weights are recommended for use in practice), and lack of validation against known causal relationships mean that results from the method should not be treated too seriously as an assessment of causality.*
"

"""

from pathlib import Path, PurePath

import numpy as np
import pandas as pd
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
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    PhenotypeInfo,
    QuantPhenotype,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.ldsc_constants import LDSC_Z_COL
from mecfs_bio.util.type_related.unwrap import unwrap

logger = structlog.get_logger()


@frozen
class SNPHeritabilityByLDSCTask(Task):
    meta: Meta
    source_sumstats_task: Task
    ld_ref_task: Task
    ld_file_filename_pattern: str
    phenotype_info: PhenotypeInfo
    pipe: DataProcessingPipe
    build: GenomeBuild
    set_N: int | None

    @property
    def _source_sumstats_id(self) -> AssetId:
        return self.source_sumstats_task.asset_id

    def _ld_ref_id(self) -> AssetId:
        return self.ld_ref_task.asset_id

    @property
    def deps(self) -> list["Task"]:
        return [self.source_sumstats_task, self.ld_ref_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        sumstats_asset = fetch(self._source_sumstats_id)
        sumstats = read_sumstats(sumstats_asset)
        sumstats.data = self.pipe.process_pandas(sumstats.data)
        sumstats.data = _drop_variants_with_degenerate_z(sumstats.data)
        ref_id = self._ld_ref_id()
        ref_asset = fetch(ref_id)
        assert isinstance(ref_asset, DirectoryAsset)
        if self.set_N is not None:
            sumstats.data[GWASLAB_SAMPLE_SIZE_COLUMN] = self.set_N
        elif GWASLAB_SAMPLE_SIZE_COLUMN not in sumstats.data.columns:
            assert self.phenotype_info.total_sample_size is not None
            sumstats.data[GWASLAB_SAMPLE_SIZE_COLUMN] = (
                self.phenotype_info.total_sample_size
            )

        sumstats.infer_build()
        assert sumstats.meta["gwaslab"]["genome_build"] == self.build
        sumstats.estimate_h2_by_ldsc(
            ref_ld_chr=str(ref_asset.path) + self.ld_file_filename_pattern,
            w_ld_chr=str(ref_asset.path) + self.ld_file_filename_pattern,
            **_get_prev_params(self.phenotype_info),
        )
        ldsc_h2 = unwrap(sumstats.ldsc_h2)
        assert isinstance(ldsc_h2, pd.DataFrame)
        out_df = ldsc_h2
        logger.debug(
            f"ldsc_h2_results:\n{out_df}",
        )
        out_path = scratch_dir / "ldsc_h2.csv"
        out_df.to_csv(out_path, index=False)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        source_sumstats_task: Task,
        ld_ref_task: Task,
        phenotype_info: PhenotypeInfo,
        pipe: DataProcessingPipe,
        build: GenomeBuild,
        set_sample_size: int | None = None,
        ld_file_filename_pattern: str = "/LDscore.@",
    ):
        sumstats_meta = source_sumstats_task.meta
        assert isinstance(sumstats_meta, GWASLabSumStatsMeta)
        meta = ResultTableMeta(
            id=asset_id,
            trait=sumstats_meta.trait,
            project=sumstats_meta.project,
            sub_dir=PurePath("analysis"),
            extension=".csv",
            read_spec=DataFrameReadSpec(DataFrameTextFormat(",")),
        )
        return cls(
            meta=meta,
            source_sumstats_task=source_sumstats_task,
            ld_ref_task=ld_ref_task,
            ld_file_filename_pattern=ld_file_filename_pattern,
            set_N=set_sample_size,
            phenotype_info=phenotype_info,
            pipe=pipe,
            build=build,
        )


def _drop_variants_with_degenerate_z(data: pd.DataFrame) -> pd.DataFrame:
    """Drop variants whose LDSC Z score would be non-finite.

    gwaslab builds the LDSC Z score as BETA / SE (or uses an existing Z column). A
    variant with SE == 0 therefore yields an infinite Z (or NaN, when BETA is also 0,
    as happens when a source reports an odds ratio rounded to 1.00). deCODE summary
    statistics contain many such variants: odds ratios rounded to 1.00 give BETA == 0
    and SE == 0, and underflowed p-values give SE == 0 with a non-zero BETA. A
    non-finite Z makes the single-trait two-step estimator's IRWLS reweighting produce
    a non-finite design matrix, which aborts the underlying SVD.

    LDSC's own munge step drops these variants; single-trait estimate_h2_by_ldsc does
    not (unlike the stratified path, which caps chi-square unconditionally), so we drop
    them here rather than at the call site, since the requirement is intrinsic to this
    regression. Variants are matched to gwaslab's Z-resolution order: prefer an existing
    Z column, otherwise derive the finiteness requirement from BETA and SE.
    """
    if LDSC_Z_COL in data.columns:
        keep = np.isfinite(data[LDSC_Z_COL])
    elif GWASLAB_BETA_COL in data.columns and GWASLAB_SE_COL in data.columns:
        keep = (
            np.isfinite(data[GWASLAB_BETA_COL])
            & np.isfinite(data[GWASLAB_SE_COL])
            & (data[GWASLAB_SE_COL] > 0)
        )
    else:
        return data
    n_dropped = int((~keep).sum())
    if n_dropped:
        logger.info(
            "Dropped variants with degenerate (non-finite) LDSC Z score",
            n_dropped=n_dropped,
            n_remaining=int(keep.sum()),
        )
    return data.loc[keep].copy()


def _get_prev_params(phenotype_info: PhenotypeInfo) -> dict:
    if isinstance(phenotype_info, QuantPhenotype):
        return {
            "samp_prev": float("nan"),
            "pop_prev": float("nan"),
        }
    return {
        "samp_prev": phenotype_info.sample_prevalence,
        "pop_prev": phenotype_info.estimated_population_prevalence,
    }
