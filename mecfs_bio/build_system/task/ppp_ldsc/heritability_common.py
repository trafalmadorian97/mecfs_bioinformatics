"""
Pieces shared by the per-protein LDSC heritability tasks (UKB-PPP, deCODE SomaScan): the run
config, the shared regression context built from a variant index and reference LD scores, the
per-protein cis-exclusion masks, and whole-table reads of dependencies.
"""

from collections.abc import Mapping, Sequence

import numpy as np
import polars as pl
from attrs import frozen

from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_h2 import DEFAULT_N_BLOCKS
from mecfs_bio.build_system.task.ppp_ldsc.ppp_ldsc_context import (
    PppLdscContext,
    build_cis_mask,
    build_ppp_ldsc_context,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)
from mecfs_bio.constants.ppp_database_constants import (
    PPP_INDEX_IS_STRAND_AMBIGUOUS_COL,
)
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

_INDEX_CONTEXT_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    PPP_INDEX_IS_STRAND_AMBIGUOUS_COL,
]


@frozen(slots=True)
class PppHeritabilityConfig:
    drop_strand_ambiguous: bool = True
    exclude_mhc: bool = True
    cis_window_bp: int = 1_000_000
    n_blocks: int = DEFAULT_N_BLOCKS
    # Proteins processed together per batched regression. Peak memory ~ n_snps * batch *
    # a few float64 arrays
    batch_size: int = 50


def build_context(
    index_task: Task,
    ld_scores_task: Task,
    config: PppHeritabilityConfig,
    fetch: Fetch,
) -> PppLdscContext:
    """The shared regression context of the variant index and the reference LD scores."""
    index_df = (
        scan_dataframe_asset(
            fetch(index_task.asset_id),
            meta=index_task.meta,
            parquet_backend="polars",
        )
        .to_native()
        .select(_INDEX_CONTEXT_COLUMNS)
        .collect()
    )
    return build_ppp_ldsc_context(
        index_df=index_df,
        ld_df=read_table(ld_scores_task, fetch),
        drop_strand_ambiguous=config.drop_strand_ambiguous,
        exclude_mhc=config.exclude_mhc,
    )


def build_cis_exclude[K](
    protein_ids: Sequence[K],
    context: PppLdscContext,
    gene_coords: Mapping[K, GenomicInterval],
    cis_window_bp: int,
) -> np.ndarray:
    """(n_snps, n_proteins) bool: True where a variant is cis to the protein's gene; all
    False for a protein without coordinates."""
    exclude = np.zeros((context.n_snps, len(protein_ids)), dtype=bool)
    for j, protein_id in enumerate(protein_ids):
        coords = gene_coords.get(protein_id)
        if coords is not None:
            exclude[:, j] = build_cis_mask(
                context,
                coords.chrom,
                coords.start,
                coords.end,
                cis_window_bp,
            )
    return exclude


def read_table(task: Task, fetch: Fetch) -> pl.DataFrame:
    """A dependency's whole dataframe asset, in memory."""
    return (
        scan_dataframe_asset(
            fetch(task.asset_id),
            meta=task.meta,
            parquet_backend="polars",
        )
        .to_native()
        .collect()
    )
