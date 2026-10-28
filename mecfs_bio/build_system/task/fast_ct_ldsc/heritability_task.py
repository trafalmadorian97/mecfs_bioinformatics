"""
Per-trait LDSC heritability and diagnostics for the traits of a fast CT-LDSC study.

Each trait's heritability is the same own-kept-set fit the all-pairs task uses
(estimate_trait_context), so the two tables agree. Columns: the alias p; observed-scale h2_obs and
its delete-block jackknife SE; for binary traits the liability-scale h2_liab and SE (the
observed-to-liability conversion is linear in h2_obs, so both are scaled by
liability_h2(1.0, info)), null for quantitative traits; the free LDSC Intercept and its SE;
lambda_gc (median kept chi-square / 0.4549); mean_chi2; n_snps (kept SNPs); n_bar (mean kept N);
and the alignment diagnostics coverage, allele_mismatch_fraction and position_concordance.
"""

from __future__ import annotations

from pathlib import Path
from typing import Self, Sequence

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.fast_ct_ldsc.aligned_trait_task import (
    FastCtLdscAlignedTraitTask,
)
from mecfs_bio.build_system.task.fast_ct_ldsc.trait_matrix import (
    TraitMatrix,
    assert_valid_study,
    estimate_trait_contexts,
    load_trait_matrix,
    multi_trait_table_meta,
    study_deps,
)
from mecfs_bio.build_system.task.fast_ct_ldsc.variant_index_task import (
    FastCtLdscVariantIndexTask,
)
from mecfs_bio.build_system.task.gwaslab.ldsc_diagnostic_plot_task import liability_h2
from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_h2 import DEFAULT_N_BLOCKS
from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_rg import TraitLdscContext
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.fast_ct_ldsc_constants import (
    FAST_CT_LDSC_ALLELE_MISMATCH_FRACTION_COL,
    FAST_CT_LDSC_COVERAGE_COL,
    FAST_CT_LDSC_H2_LIAB_COL,
    FAST_CT_LDSC_H2_LIAB_SE_COL,
    FAST_CT_LDSC_H2_OBS_COL,
    FAST_CT_LDSC_H2_OBS_SE_COL,
    FAST_CT_LDSC_INTERCEPT_COL,
    FAST_CT_LDSC_INTERCEPT_SE_COL,
    FAST_CT_LDSC_LAMBDA_GC_COL,
    FAST_CT_LDSC_MEAN_CHI2_COL,
    FAST_CT_LDSC_N_BAR_COL,
    FAST_CT_LDSC_N_SNPS_COL,
    FAST_CT_LDSC_POSITION_CONCORDANCE_COL,
    FAST_CT_LDSC_TRAIT_COL,
)


@frozen(slots=True)
class FastCtLdscHeritabilityTask(Task):
    """One row of heritability and alignment diagnostics per aligned trait."""

    meta: Meta
    index_task: FastCtLdscVariantIndexTask
    aligned_tasks: tuple[FastCtLdscAlignedTraitTask, ...]
    n_blocks: int = DEFAULT_N_BLOCKS

    def __attrs_post_init__(self) -> None:
        assert_valid_study(self.aligned_tasks)

    @property
    def deps(self) -> list[Task]:
        return study_deps(self.index_task, self.aligned_tasks)

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        matrix = load_trait_matrix(self.index_task, self.aligned_tasks, fetch)
        contexts = estimate_trait_contexts(matrix, n_blocks=self.n_blocks)
        out_path = scratch_dir / "fast_ct_ldsc_heritability.parquet"
        heritability_table(matrix, contexts).write_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        *,
        index_task: FastCtLdscVariantIndexTask,
        aligned_tasks: Sequence[FastCtLdscAlignedTraitTask],
        n_blocks: int = DEFAULT_N_BLOCKS,
    ) -> Self:
        return cls(
            meta=multi_trait_table_meta(asset_id),
            index_task=index_task,
            aligned_tasks=tuple(aligned_tasks),
            n_blocks=n_blocks,
        )


def heritability_table(
    matrix: TraitMatrix, contexts: Sequence[TraitLdscContext]
) -> pl.DataFrame:
    rows = []
    for alias, ctx, info, diag in zip(
        matrix.aliases, contexts, matrix.sample_infos, matrix.diagnostics
    ):
        factor = liability_h2(1.0, info)
        rows.append(
            {
                FAST_CT_LDSC_TRAIT_COL: alias,
                FAST_CT_LDSC_H2_OBS_COL: ctx.h2,
                FAST_CT_LDSC_H2_OBS_SE_COL: ctx.h2_se,
                FAST_CT_LDSC_H2_LIAB_COL: None if factor is None else factor * ctx.h2,
                FAST_CT_LDSC_H2_LIAB_SE_COL: None
                if factor is None
                else factor * ctx.h2_se,
                FAST_CT_LDSC_INTERCEPT_COL: ctx.intercept,
                FAST_CT_LDSC_INTERCEPT_SE_COL: ctx.intercept_se,
                FAST_CT_LDSC_LAMBDA_GC_COL: ctx.lambda_gc,
                FAST_CT_LDSC_MEAN_CHI2_COL: ctx.mean_chi2,
                FAST_CT_LDSC_N_SNPS_COL: ctx.n_snps,
                FAST_CT_LDSC_N_BAR_COL: ctx.n_bar,
                FAST_CT_LDSC_COVERAGE_COL: diag.coverage,
                FAST_CT_LDSC_ALLELE_MISMATCH_FRACTION_COL: diag.allele_mismatch_fraction,
                FAST_CT_LDSC_POSITION_CONCORDANCE_COL: diag.position_concordance,
            }
        )
    schema_overrides = {
        FAST_CT_LDSC_H2_LIAB_COL: pl.Float64,
        FAST_CT_LDSC_H2_LIAB_SE_COL: pl.Float64,
        FAST_CT_LDSC_POSITION_CONCORDANCE_COL: pl.Float64,
    }
    return pl.DataFrame(rows, schema_overrides=schema_overrides)
