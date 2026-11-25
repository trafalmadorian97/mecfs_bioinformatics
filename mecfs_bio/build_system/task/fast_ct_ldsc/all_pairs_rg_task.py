"""
All-pairs cross-trait LDSC genetic correlation for a set of aligned traits.

Every trait is aligned to the same context SNP set, so the study is a stack of z-scores
Z (S, T) with sample sizes N (S, T) (see trait_matrix). Each trait's heritability is estimated
once on its own kept set (estimate_trait_context). Then, for each trait i, the covariances with
all later traits j > i are computed in column batches of batch_size by the batched rg kernel
(batched_rg), which also re-estimates each h2_j on its own kept set, so the diagonal h2 is the
same in every pair. batch_size bounds the kernel's (S, batch_size) temporaries and hence peak
memory.

Output: one row per unordered pair, p1 = trait i and p2 = trait j (traits ordered by lower-cased
alias), with rg and its delete-block jackknife SE / z / p, the observed-scale genetic covariance
gcov and the cross-trait intercept gcov_int (each with its jackknife SE), both observed-scale
heritabilities, h2_int (the intercept of p2's own heritability regression, which LAVA's
sample-overlap reader uses for the diagonal), and n_snps (SNPs entering the covariance). rg is
NaN where either heritability is not positive.
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
from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_h2 import DEFAULT_N_BLOCKS
from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_rg import (
    BatchedRgResult,
    TraitLdscContext,
    batched_rg,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.fast_ct_ldsc_constants import (
    FAST_CT_LDSC_GCOV_COL,
    FAST_CT_LDSC_GCOV_INT_COL,
    FAST_CT_LDSC_GCOV_INT_SE_COL,
    FAST_CT_LDSC_GCOV_SE_COL,
    FAST_CT_LDSC_H2_INT_COL,
    FAST_CT_LDSC_H2_OBS_P1_COL,
    FAST_CT_LDSC_H2_OBS_P2_COL,
    FAST_CT_LDSC_N_SNPS_COL,
    FAST_CT_LDSC_P1_COL,
    FAST_CT_LDSC_P2_COL,
    FAST_CT_LDSC_RG_COL,
    FAST_CT_LDSC_RG_P_COL,
    FAST_CT_LDSC_RG_SE_COL,
    FAST_CT_LDSC_RG_Z_COL,
)

DEFAULT_RG_BATCH_SIZE = 8


@frozen(slots=True)
class FastCtLdscAllPairsRgTask(Task):
    """Genetic correlation for every unordered pair of the aligned traits."""

    meta: Meta
    index_task: FastCtLdscVariantIndexTask
    aligned_tasks: tuple[FastCtLdscAlignedTraitTask, ...]
    n_blocks: int = DEFAULT_N_BLOCKS
    batch_size: int = DEFAULT_RG_BATCH_SIZE

    def __attrs_post_init__(self) -> None:
        assert_valid_study(self.aligned_tasks)
        assert self.batch_size >= 1

    @property
    def deps(self) -> list[Task]:
        return study_deps(self.index_task, self.aligned_tasks)

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        matrix = load_trait_matrix(self.index_task, self.aligned_tasks, fetch)
        contexts = estimate_trait_contexts(matrix, n_blocks=self.n_blocks)
        table = all_pairs_table(
            matrix, contexts, n_blocks=self.n_blocks, batch_size=self.batch_size
        )
        out_path = scratch_dir / "fast_ct_ldsc_pairs.parquet"
        table.write_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        *,
        index_task: FastCtLdscVariantIndexTask,
        aligned_tasks: Sequence[FastCtLdscAlignedTraitTask],
        n_blocks: int = DEFAULT_N_BLOCKS,
        batch_size: int = DEFAULT_RG_BATCH_SIZE,
    ) -> Self:
        return cls(
            meta=multi_trait_table_meta(asset_id),
            index_task=index_task,
            aligned_tasks=tuple(aligned_tasks),
            n_blocks=n_blocks,
            batch_size=batch_size,
        )


def all_pairs_table(
    matrix: TraitMatrix,
    contexts: Sequence[TraitLdscContext],
    *,
    n_blocks: int,
    batch_size: int,
) -> pl.DataFrame:
    """One row per unordered pair (i < j), batching the later traits j for each trait i."""
    rows: list[dict[str, str | float | int]] = []
    n_traits = len(matrix.aliases)
    for i in range(n_traits):
        for start in range(i + 1, n_traits, batch_size):
            cols = list(range(start, min(start + batch_size, n_traits)))
            result = batched_rg(
                contexts[i],
                z_protein=matrix.z[:, cols],
                ld=matrix.context.ld,
                n_protein=matrix.n[:, cols],
                m=matrix.context.m,
                n_blocks=n_blocks,
            )
            rows.extend(_pair_rows(matrix, contexts, i, cols, result))
    return pl.DataFrame(rows)


def _pair_rows(
    matrix: TraitMatrix,
    contexts: Sequence[TraitLdscContext],
    i: int,
    cols: list[int],
    result: BatchedRgResult,
) -> list[dict[str, str | float | int]]:
    return [
        {
            FAST_CT_LDSC_P1_COL: matrix.aliases[i],
            FAST_CT_LDSC_P2_COL: matrix.aliases[j],
            FAST_CT_LDSC_RG_COL: float(result.rg[k]),
            FAST_CT_LDSC_RG_SE_COL: float(result.rg_se[k]),
            FAST_CT_LDSC_RG_Z_COL: float(result.rg_z[k]),
            FAST_CT_LDSC_RG_P_COL: float(result.rg_p[k]),
            FAST_CT_LDSC_GCOV_COL: float(result.gcov[k]),
            FAST_CT_LDSC_GCOV_SE_COL: float(result.gcov_se[k]),
            FAST_CT_LDSC_GCOV_INT_COL: float(result.gcov_intercept[k]),
            FAST_CT_LDSC_GCOV_INT_SE_COL: float(result.gcov_intercept_se[k]),
            FAST_CT_LDSC_H2_OBS_P1_COL: contexts[i].h2,
            FAST_CT_LDSC_H2_OBS_P2_COL: contexts[j].h2,
            FAST_CT_LDSC_H2_INT_COL: contexts[j].intercept,
            FAST_CT_LDSC_N_SNPS_COL: int(result.n_snps[k]),
        }
        for k, j in enumerate(cols)
    ]
