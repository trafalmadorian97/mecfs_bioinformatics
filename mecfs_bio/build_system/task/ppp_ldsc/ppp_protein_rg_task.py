"""
Cross-trait genetic correlation of one trait against every UKB-PPP protein via LD-score regression.

Given a trait's hg38 summary statistics, this computes the genetic correlation (rg) of the trait
with each protein for a single variant set -- all variants, or cis-excluded (variants within a
window of the protein's gene, from Sun et al. 2023 ST3 hg38 gene coordinates) -- selected by the
config. Run the task twice to get both variant sets; each run's asset id encodes the variant set.

All proteins share one variant-index row order, so the index<->LD-score alignment, LD scores, M and
jackknife blocks are built once (PppLdscContext), the trait is aligned onto that set once
(align_trait_to_context) and its heritability estimated once (estimate_trait_context), and the
per-protein work is a cheap batched weighted regression pair -- h2_protein and the trait-protein
covariance (batched_rg). Following the GenomicSEM convention, the trait heritability is computed on
its own kept set and reused for every protein.

Output: one row per protein with rg, its jackknife SE / z / p, the genetic covariance and its
cross-trait intercept, and the two heritabilities. No significance filtering is applied.
"""

from __future__ import annotations

from pathlib import Path, PurePath

import numpy as np
import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.task.ppp_database.build_slim_protein_parquet_task import (
    BuildSlimProteinParquetTask,
)
from mecfs_bio.build_system.task.ppp_database.protein_sample_size_task import (
    PppProteinSampleSizeTask,
)
from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_rg import (
    BatchedRgResult,
    TraitLdscContext,
    batched_rg,
    estimate_trait_context,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_ldsc_context import (
    PppLdscContext,
    build_cis_mask,
    build_ppp_ldsc_context,
    read_ld_scores,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_heritability_task import (
    ProteinSampleInfo,
    read_st3_gene_coords,
)
from mecfs_bio.build_system.task.ppp_ldsc.trait_alignment import align_trait_to_context
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
from mecfs_bio.constants.ppp_index_constants import PPP_INDEX_IS_STRAND_AMBIGUOUS_COL
from mecfs_bio.constants.ppp_ldsc_constants import (
    PPP_N_GENE_COL,
    PPP_N_OID_COL,
    PPP_N_SAMPLE_SIZE_COL,
    PPP_N_SYNID_COL,
    PPP_RG_GCOV_COL,
    PPP_RG_GCOV_INTERCEPT_COL,
    PPP_RG_GENE_COL,
    PPP_RG_H2_PROTEIN_COL,
    PPP_RG_H2_TRAIT_COL,
    PPP_RG_N_BAR_PROTEIN_COL,
    PPP_RG_N_BAR_TRAIT_COL,
    PPP_RG_N_SNPS_COL,
    PPP_RG_OID_COL,
    PPP_RG_RG_COL,
    PPP_RG_RG_P_COL,
    PPP_RG_RG_SE_COL,
    PPP_RG_RG_Z_COL,
    PPP_RG_VARIANT_SET_COL,
    PPP_VARIANT_SET_ALL,
    PPP_VARIANT_SET_CIS_EXCLUDED,
    PppVariantSet,
)
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

logger = structlog.get_logger()

# Index columns read for the context: the context builder needs CHR/POS/rsID/strand-ambiguity,
# and trait alignment additionally needs the index EA/NEA to orient the trait's z-score.
_INDEX_CONTEXT_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    PPP_INDEX_IS_STRAND_AMBIGUOUS_COL,
]
_CONTEXT_VARIANT_COLUMNS = [
    GWASLAB_RSID_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
]


@frozen
class PppRgConfig:
    variant_set: PppVariantSet = PPP_VARIANT_SET_ALL
    # Constant trait N used only when the trait sumstats lack a per-SNP N column.
    trait_total_sample_size: int | None = None
    drop_strand_ambiguous: bool = True
    exclude_mhc: bool = True
    cis_window_bp: int = 1_000_000
    n_blocks: int = 200
    # Proteins processed together per batched regression.
    batch_size: int = 50
    # Fail if fewer than this many context SNPs match the trait (harmonization guard).
    min_trait_snps: int = 200_000


@frozen
class PppProteinCrossTraitRgTask(GeneratingTask):
    """Compute the trait-vs-protein genetic correlation for every UKB-PPP protein.

    trait_task: a task producing the trait's hg38 sumstats dataframe (gwaslab columns rsID, EA,
        NEA, BETA, SE, and optionally N).
    protein_tasks: the slim per-protein beta/se tasks (aligned to index_task).
    index_task: the shared ConstructPppVariantIndexTask (variant identity / rsID / alleles).
    ld_ref_task: an extracted reference LD-score directory (LDscore.<chr>.l2.* files).
    sample_size_task: the per-protein N table (PppProteinSampleSizeTask).
    st3_task: Sun et al. 2023 supplementary xlsx (ST3 hg38 gene coordinates).
    """

    meta: Meta
    trait_task: Task
    protein_tasks: tuple[BuildSlimProteinParquetTask, ...]
    index_task: Task
    ld_ref_task: Task
    sample_size_task: PppProteinSampleSizeTask
    st3_task: Task
    config: PppRgConfig

    @property
    def deps(self) -> list[Task]:
        return [
            self.trait_task,
            *self.protein_tasks,
            self.index_task,
            self.ld_ref_task,
            self.sample_size_task,
            self.st3_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        index_df = _read_index(self.index_task, fetch)
        context = build_ppp_ldsc_context(
            index_df,
            *_read_ld(self.ld_ref_task, fetch),
            drop_strand_ambiguous=self.config.drop_strand_ambiguous,
            exclude_mhc=self.config.exclude_mhc,
        )
        context_variants = index_df.select(_CONTEXT_VARIANT_COLUMNS).gather(
            pl.Series(context.row_pos)
        )
        trait_ctx = _build_trait_context(
            self.trait_task, context, context_variants, self.config, fetch
        )
        sample_sizes = _read_sample_sizes(self.sample_size_task, fetch)
        gene_coords = read_st3_gene_coords(_st3_path(self.st3_task, fetch))

        tasks = _proteins_to_process(
            self.protein_tasks, sample_sizes, gene_coords, self.config.variant_set
        )
        logger.info(
            "built ppp cross-trait rg context",
            n_snps=context.n_snps,
            m=context.m,
            h2_trait=trait_ctx.h2,
            variant_set=self.config.variant_set,
            n_proteins=len(tasks),
        )

        rows: list[dict] = []
        for start in range(0, len(tasks), self.config.batch_size):
            batch = tasks[start : start + self.config.batch_size]
            rows.extend(
                _process_batch(
                    batch,
                    context,
                    trait_ctx,
                    sample_sizes,
                    gene_coords,
                    self.config,
                    fetch,
                )
            )

        table = pl.DataFrame(rows)
        out_path = scratch_dir / f"{self.meta.asset_id}.parquet"
        table.write_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        trait_task: Task,
        protein_tasks: tuple[BuildSlimProteinParquetTask, ...],
        index_task: Task,
        ld_ref_task: Task,
        sample_size_task: PppProteinSampleSizeTask,
        st3_task: Task,
        config: PppRgConfig = PppRgConfig(),
    ) -> PppProteinCrossTraitRgTask:
        # The trait identity comes from the trait task's meta (never passed in), per convention.
        trait_meta = trait_task.meta
        assert isinstance(
            trait_meta, (FilteredGWASDataMeta, GWASSummaryDataFileMeta)
        ), f"trait_task must produce a GWAS dataframe meta, got {type(trait_meta)}"
        meta = ResultTableMeta(
            id=f"{asset_id}_{config.variant_set}",
            trait=trait_meta.trait,
            project="ppp_genetic_correlation",
            sub_dir=PurePath("analysis"),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
        return cls(
            meta=meta,
            trait_task=trait_task,
            protein_tasks=protein_tasks,
            index_task=index_task,
            ld_ref_task=ld_ref_task,
            sample_size_task=sample_size_task,
            st3_task=st3_task,
            config=config,
        )


def _read_index(index_task: Task, fetch: Fetch) -> pl.DataFrame:
    return (
        scan_dataframe_asset(
            fetch(index_task.asset_id),
            meta=index_task.meta,
            parquet_backend="polars",
        )
        .to_native()
        .select(_INDEX_CONTEXT_COLUMNS)
        .collect()
    )


def _read_ld(ld_ref_task: Task, fetch: Fetch) -> tuple[pl.DataFrame, float]:
    ld_asset = fetch(ld_ref_task.asset_id)
    assert isinstance(ld_asset, DirectoryAsset)
    return read_ld_scores(ld_asset.path)


def _build_trait_context(
    trait_task: Task,
    context: PppLdscContext,
    context_variants: pl.DataFrame,
    config: PppRgConfig,
    fetch: Fetch,
) -> TraitLdscContext:
    trait_df = (
        scan_dataframe_asset(
            fetch(trait_task.asset_id),
            meta=trait_task.meta,
            parquet_backend="polars",
        )
        .to_native()
        .collect()
    )
    aligned = align_trait_to_context(
        trait_df,
        context_variants,
        trait_total_sample_size=config.trait_total_sample_size,
        min_trait_snps=config.min_trait_snps,
    )
    return estimate_trait_context(
        aligned.z, aligned.n, context.ld, context.m, n_blocks=config.n_blocks
    )


def _read_sample_sizes(
    sample_size_task: PppProteinSampleSizeTask, fetch: Fetch
) -> dict[str, ProteinSampleInfo]:
    """Synapse id -> ProteinSampleInfo (oid, gene, n) from the sample-size table."""
    table = (
        scan_dataframe_asset(
            fetch(sample_size_task.asset_id),
            meta=sample_size_task.meta,
            parquet_backend="polars",
        )
        .to_native()
        .collect()
    )
    return {
        row[PPP_N_SYNID_COL]: ProteinSampleInfo(
            oid=row[PPP_N_OID_COL],
            gene=row[PPP_N_GENE_COL],
            n=int(row[PPP_N_SAMPLE_SIZE_COL]),
        )
        for row in table.iter_rows(named=True)
    }


def _st3_path(st3_task: Task, fetch: Fetch) -> Path:
    asset = fetch(st3_task.asset_id)
    assert isinstance(asset, FileAsset)
    return asset.path


def _proteins_to_process(
    protein_tasks: tuple[BuildSlimProteinParquetTask, ...],
    sample_sizes: dict[str, ProteinSampleInfo],
    gene_coords: dict[str, GenomicInterval],
    variant_set: PppVariantSet,
) -> list[BuildSlimProteinParquetTask]:
    """In cis-excluded mode, skip proteins with no ST3 gene coordinates (no cis window to
    exclude); in all-variants mode, keep every protein."""
    if variant_set == PPP_VARIANT_SET_ALL:
        return list(protein_tasks)
    return [
        task for task in protein_tasks if sample_sizes[task.synid].oid in gene_coords
    ]


def _process_batch(
    batch: list[BuildSlimProteinParquetTask],
    context: PppLdscContext,
    trait_ctx: TraitLdscContext,
    sample_sizes: dict[str, ProteinSampleInfo],
    gene_coords: dict[str, GenomicInterval],
    config: PppRgConfig,
    fetch: Fetch,
) -> list[dict]:
    oids, genes, n_vec, z_cols = [], [], [], []
    for protein_task in batch:
        info = sample_sizes[protein_task.synid]
        oids.append(info.oid)
        genes.append(info.gene)
        n_vec.append(float(info.n))
        z_cols.append(_read_protein_z(protein_task, context, fetch))
    z_protein = np.column_stack(z_cols)
    n_arr = np.asarray(n_vec)

    exclude = (
        _build_cis_exclude(oids, context, gene_coords, config.cis_window_bp)
        if config.variant_set == PPP_VARIANT_SET_CIS_EXCLUDED
        else None
    )
    res = batched_rg(
        trait_ctx,
        z_protein,
        context.ld,
        n_arr,
        context.m,
        n_blocks=config.n_blocks,
        exclude=exclude,
    )
    return [
        _make_row(oids[j], genes[j], config.variant_set, res, j)
        for j in range(len(oids))
    ]


def _read_protein_z(
    protein_task: BuildSlimProteinParquetTask, context: PppLdscContext, fetch: Fetch
) -> np.ndarray:
    """Signed protein z-score (BETA/SE, index-effect-allele oriented) at the context SNPs."""
    frame = (
        scan_dataframe_asset(
            fetch(protein_task.asset_id),
            meta=protein_task.meta,
            parquet_backend="polars",
        )
        .to_native()
        .select(GWASLAB_BETA_COL, GWASLAB_SE_COL)
        .collect()
    )
    beta = frame[GWASLAB_BETA_COL].to_numpy().astype(float)[context.row_pos]
    se = frame[GWASLAB_SE_COL].to_numpy().astype(float)[context.row_pos]
    with np.errstate(invalid="ignore", divide="ignore"):
        return beta / se


def _build_cis_exclude(
    oids: list[str],
    context: PppLdscContext,
    gene_coords: dict[str, GenomicInterval],
    cis_window_bp: int,
) -> np.ndarray:
    exclude = np.zeros((context.n_snps, len(oids)), dtype=bool)
    for j, oid in enumerate(oids):
        coords = gene_coords.get(oid)
        if coords is not None:
            exclude[:, j] = build_cis_mask(
                context, coords.chrom, coords.start, coords.end, cis_window_bp
            )
    return exclude


def _make_row(
    oid: str,
    gene: str,
    variant_set: str,
    res: BatchedRgResult,
    j: int,
) -> dict:
    return {
        PPP_RG_OID_COL: oid,
        PPP_RG_GENE_COL: gene,
        PPP_RG_VARIANT_SET_COL: variant_set,
        PPP_RG_RG_COL: float(res.rg[j]),
        PPP_RG_RG_SE_COL: float(res.rg_se[j]),
        PPP_RG_RG_Z_COL: float(res.rg_z[j]),
        PPP_RG_RG_P_COL: float(res.rg_p[j]),
        PPP_RG_GCOV_COL: float(res.gcov[j]),
        PPP_RG_GCOV_INTERCEPT_COL: float(res.gcov_intercept[j]),
        PPP_RG_H2_TRAIT_COL: float(res.h2_trait),
        PPP_RG_H2_PROTEIN_COL: float(res.h2_protein[j]),
        PPP_RG_N_SNPS_COL: int(res.n_snps[j]),
        PPP_RG_N_BAR_TRAIT_COL: float(res.n_bar_trait),
        PPP_RG_N_BAR_PROTEIN_COL: float(res.n_bar_protein[j]),
    }
