"""
Cell-type-specific stratified LD score regression ("Heritability enrichment of specifically
expressed genes identifies disease-relevant tissues and cell types", Finucane et al. 2018) over
every model of a prepared annotation set at once.

weighting="exact" reproduces ldsc's cell_type_specific (as run by CellAnalysisByLDSCTask) to
floating-point tolerance; "shared_controls" is a faster approximation (see batched_sldsc_cts).
The annotation set may come from any source with an adapter to the prepared S-LDSC contract.
"""

import gc
from pathlib import Path, PurePath

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.sldsc.batched_sldsc_cts import (
    BatchedSLDSCResult,
    batched_sldsc_cts,
)
from mecfs_bio.build_system.task.sldsc.prepared_sldsc_reference import (
    read_prepared_annotation_set,
    read_prepared_baseline,
)
from mecfs_bio.build_system.task.sldsc.sldsc_regression_inputs import (
    AlignedSLDSCRegression,
    build_regression_inputs,
    load_sldsc_sumstats,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.sldsc_constants import (
    SLDSC_CTS_COEF_COL,
    SLDSC_CTS_COEF_SE_COL,
    SLDSC_CTS_NAME_COL,
    SLDSC_CTS_P_COL,
    SLDSC_DEFAULT_MIN_HAPMAP3_SNPS,
    SLDSCWeighting,
)

_RESULT_FILENAME = "ldsc_h2_cts.csv"


@frozen(slots=True)
class VectorizedCellAnalysisBySLDSCTask(Task):
    """
    source_sumstats_task: gwaslab sumstats (pickle).
    baseline_task: prepared S-LDSC baseline (PrepareSLDSCBaselineTask).
    annotation_set_task: prepared S-LDSC annotation set (e.g. PrepareLdctsAnnotationSetTask).
    weighting: exact (default) or shared_controls.
    pre_pipe: applied to the sumstats before analysis.
    min_hapmap3_snps: minimum number of sumstats SNPs in the baseline's HapMap3 regression-SNP
        set; fewer raises, since it usually means an upstream processing problem.
    """

    meta: Meta
    source_sumstats_task: Task
    baseline_task: Task
    annotation_set_task: Task
    weighting: SLDSCWeighting = "exact"
    pre_pipe: DataProcessingPipe = IdentityPipe()
    min_hapmap3_snps: int = SLDSC_DEFAULT_MIN_HAPMAP3_SNPS

    @property
    def deps(self) -> list[Task]:
        return [self.source_sumstats_task, self.baseline_task, self.annotation_set_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        aligned = _load_aligned_regression(
            sumstats_asset=fetch(self.source_sumstats_task.asset_id),
            pre_pipe=self.pre_pipe,
            baseline_dir=_directory(fetch, self.baseline_task),
            annotation_set_dir=_directory(fetch, self.annotation_set_task),
            min_hapmap3_snps=self.min_hapmap3_snps,
        )
        result = batched_sldsc_cts(
            inputs=aligned.inputs, models=aligned.models, weighting=self.weighting
        )
        out_path = scratch_dir / _RESULT_FILENAME
        _result_frame(result).write_csv(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        source_sumstats_task: Task,
        baseline_task: Task,
        annotation_set_task: Task,
        weighting: SLDSCWeighting = "exact",
        pre_pipe: DataProcessingPipe = IdentityPipe(),
        min_hapmap3_snps: int = SLDSC_DEFAULT_MIN_HAPMAP3_SNPS,
    ) -> "VectorizedCellAnalysisBySLDSCTask":
        sumstats_meta = source_sumstats_task.meta
        assert isinstance(sumstats_meta, GWASLabSumStatsMeta)
        return cls(
            meta=ResultTableMeta(
                id=asset_id,
                trait=sumstats_meta.trait,
                project=sumstats_meta.project,
                sub_dir=PurePath("analysis"),
                extension=".csv",
                read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
            ),
            source_sumstats_task=source_sumstats_task,
            baseline_task=baseline_task,
            annotation_set_task=annotation_set_task,
            weighting=weighting,
            pre_pipe=pre_pipe,
            min_hapmap3_snps=min_hapmap3_snps,
        )


def _load_aligned_regression(
    sumstats_asset: Asset,
    pre_pipe: DataProcessingPipe,
    baseline_dir: Path,
    annotation_set_dir: Path,
    min_hapmap3_snps: int,
) -> AlignedSLDSCRegression:
    """
    Kernel inputs aligned to the regression SNPs. The full prepared references (the annotation
    LD matrix alone is ~2.3 GB for Multi_tissue_chromatin) are local to this function, so they
    are freed before the kernel runs rather than held alongside its row-aligned copies.
    """
    aligned = build_regression_inputs(
        sumstats=load_sldsc_sumstats(asset=sumstats_asset, pre_pipe=pre_pipe),
        baseline=read_prepared_baseline(baseline_dir),
        annotation_set=read_prepared_annotation_set(annotation_set_dir),
        min_hapmap3_snps=min_hapmap3_snps,
    )
    # The loaded gwaslab Sumstats object sits in reference cycles, so its full dataframe is
    # only freed by the cyclic garbage collector.
    gc.collect()
    return aligned


def _directory(fetch: Fetch, task: Task) -> Path:
    asset = fetch(task.asset_id)
    assert isinstance(asset, DirectoryAsset)
    return asset.path


def _result_frame(result: BatchedSLDSCResult) -> pl.DataFrame:
    return pl.DataFrame(
        {
            SLDSC_CTS_NAME_COL: list(result.labels),
            SLDSC_CTS_COEF_COL: result.coef,
            SLDSC_CTS_COEF_SE_COL: result.coef_se,
            SLDSC_CTS_P_COL: result.p,
        }
    ).sort(SLDSC_CTS_P_COL, maintain_order=True)
