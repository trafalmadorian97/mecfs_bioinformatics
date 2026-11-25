"""
Stack the aligned traits of a fast CT-LDSC study into (S, T) matrices over the shared context.

Glossary: S is the number of context (regression) SNPs and T the number of traits. Column t of z
and n holds trait t's aligned z-score and per-SNP sample size (NaN where the trait lacks the
variant), in context row order.
"""

from __future__ import annotations

from pathlib import PurePath
from typing import Sequence

import numpy as np
from attrs import frozen

from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.fast_ct_ldsc.aligned_trait_task import (
    FastCtLdscAlignedTraitTask,
    read_aligned_trait,
)
from mecfs_bio.build_system.task.fast_ct_ldsc.alignment_guards import (
    AlignmentDiagnostics,
)
from mecfs_bio.build_system.task.fast_ct_ldsc.variant_index_task import (
    FastCtLdscVariantIndexTask,
    load_indexed_ldsc_context,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    MULTI_TRAIT,
    PhenotypeInfo,
)
from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_rg import (
    TraitLdscContext,
    estimate_trait_context,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_ldsc_context import PppLdscContext
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL


@frozen(slots=True)
class TraitMatrix:
    """The study's traits over the shared context: z and n are (S, T); aliases, diagnostics and
    sample_infos are length T, in column order."""

    context: PppLdscContext
    aliases: tuple[str, ...]
    z: np.ndarray
    n: np.ndarray
    diagnostics: tuple[AlignmentDiagnostics, ...]
    sample_infos: tuple[PhenotypeInfo, ...]

    def __attrs_post_init__(self) -> None:
        shape = (self.context.n_snps, len(self.aliases))
        for name, arr in (("z", self.z), ("n", self.n)):
            assert arr.shape == shape, f"{name} has shape {arr.shape}, expected {shape}"
            assert arr.dtype.kind == "f", f"{name} must be float, got {arr.dtype}"
        assert len(self.diagnostics) == len(self.aliases)
        assert len(self.sample_infos) == len(self.aliases)


def load_trait_matrix(
    index_task: FastCtLdscVariantIndexTask,
    aligned_tasks: Sequence[FastCtLdscAlignedTraitTask],
    fetch: Fetch,
) -> TraitMatrix:
    """Read every aligned trait and stack it, asserting each was aligned to this context."""
    indexed = load_indexed_ldsc_context(index_task, fetch)
    context_rsids = indexed.variants[GWASLAB_RSID_COL]
    arrays = [read_aligned_trait(task, fetch) for task in aligned_tasks]
    for task, aligned in zip(aligned_tasks, arrays):
        assert aligned.rsid.equals(context_rsids), (
            f"{task.alias}: aligned rsIDs differ from the context's; the aligned trait was "
            "built against a different variant index"
        )
    return TraitMatrix(
        context=indexed.context,
        aliases=tuple(task.alias for task in aligned_tasks),
        z=np.column_stack([a.z for a in arrays]),
        n=np.column_stack([a.n for a in arrays]),
        diagnostics=tuple(a.diagnostics for a in arrays),
        sample_infos=tuple(task.source.sample_info for task in aligned_tasks),
    )


def estimate_trait_contexts(
    matrix: TraitMatrix, *, n_blocks: int
) -> list[TraitLdscContext]:
    """Each trait's own-set LDSC heritability fit (estimate_trait_context), in column order."""
    return [
        estimate_trait_context(
            z_trait=matrix.z[:, t],
            n_trait=matrix.n[:, t],
            ld=matrix.context.ld,
            m=matrix.context.m,
            n_blocks=n_blocks,
        )
        for t in range(len(matrix.aliases))
    ]


def assert_valid_study(aligned_tasks: Sequence[FastCtLdscAlignedTraitTask]) -> None:
    """A study needs at least two traits with unique aliases, ordered by lower-cased alias."""
    aliases = [task.alias for task in aligned_tasks]
    assert len(aliases) >= 2, "an all-pairs study needs at least two traits"
    assert len(set(aliases)) == len(aliases), f"duplicate trait aliases in {aliases}"
    assert aliases == sorted(aliases, key=str.lower), (
        f"aligned tasks must be ordered by lower-cased alias, got {aliases}"
    )


def study_deps(
    index_task: FastCtLdscVariantIndexTask,
    aligned_tasks: Sequence[FastCtLdscAlignedTraitTask],
) -> list[Task]:
    return [index_task, index_task.ld_scores_task, *aligned_tasks]


def multi_trait_table_meta(asset_id: str) -> ResultTableMeta:
    return ResultTableMeta(
        id=asset_id,
        trait=MULTI_TRAIT,
        project="genetic_correlation",
        sub_dir=PurePath("fast_ct_ldsc"),
        extension=".parquet",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    )
