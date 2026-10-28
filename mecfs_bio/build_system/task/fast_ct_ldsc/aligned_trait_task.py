"""
Align one trait's z-scores to the fast CT-LDSC regression SNPs, once, behind upstream-bug guards.

The trait is read either from a pickled gwaslab Sumstats object or from a dataframe asset, passed
through the source's pipe, and aligned by rsID to the context derived from the shared variant
index (align_trait_to_context): z = BETA/SE oriented to the index effect allele, NaN where the
trait lacks the variant or its alleles do not match. Before anything is written, the alignment
diagnostics are checked against the guards (see alignment_guards); a source's own guards, when
set, replace the generator's. The task also fails when a variant has a finite z but no sample
size, which would otherwise poison the heritability and N_bar.

The trait's declared genome build, which decides whether the position guard applies, is the
source meta's genome_build when set; otherwise, for a gwaslab pickle, the build the pickle
records (when it is 19 or 38); otherwise unknown. A pickle whose meta and recorded builds
disagree is itself an upstream bug and fails.

Output: a directory holding aligned.parquet (rsID, z, n over the context SNPs, in context row
order) and diagnostics.json (the alignment diagnostics).
"""

from __future__ import annotations

import json
from pathlib import Path, PurePath
from typing import Self

import narwhals
import numpy as np
import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.read_spec.read_sumstats import read_sumstats
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.fast_ct_ldsc.alignment_guards import (
    AlignmentDiagnostics,
    AlignmentGuards,
    alignment_diagnostics,
    check_alignment_guards,
)
from mecfs_bio.build_system.task.fast_ct_ldsc.variant_index_task import (
    FastCtLdscVariantIndexTask,
    load_indexed_ldsc_context,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    SumstatsSource,
)
from mecfs_bio.build_system.task.ppp_ldsc.trait_alignment import (
    TraitAligned,
    align_trait_to_context,
)
from mecfs_bio.build_system.task.task_util import produces_dataframe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.fast_ct_ldsc_constants import (
    FAST_CT_LDSC_ALIGNED_FILENAME,
    FAST_CT_LDSC_DIAGNOSTICS_FILENAME,
    FAST_CT_LDSC_N_COL,
    FAST_CT_LDSC_Z_COL,
)
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

_REQUIRED_TRAIT_COLUMNS = (
    GWASLAB_RSID_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_BETA_COL,
    GWASLAB_SE_COL,
)
_OPTIONAL_TRAIT_COLUMNS = (
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
)


@frozen(slots=True)
class LoadedTrait:
    """A trait's sumstats as a lazy frame, plus its declared genome build (None if unknown)."""

    frame: narwhals.LazyFrame
    declared_build: GenomeBuild | None


@frozen(slots=True)
class AlignedTraitArrays:
    """A written aligned trait read back: rsID, z and n over the context SNPs, in context row
    order, and the alignment diagnostics."""

    rsid: pl.Series
    z: np.ndarray
    n: np.ndarray
    diagnostics: AlignmentDiagnostics

    def __attrs_post_init__(self) -> None:
        s = len(self.rsid)
        for name, arr in (("z", self.z), ("n", self.n)):
            assert arr.shape == (s,), f"{name} has shape {arr.shape}, expected ({s},)"
            assert arr.dtype.kind == "f", f"{name} must be float, got {arr.dtype}"


@frozen(slots=True)
class FastCtLdscAlignedTraitTask(Task):
    """Align one trait to the fast CT-LDSC context and enforce the alignment guards."""

    meta: Meta
    source: SumstatsSource
    index_task: FastCtLdscVariantIndexTask
    guards: AlignmentGuards

    @property
    def alias(self) -> str:
        return self.source.alias

    @property
    def deps(self) -> list[Task]:
        return [self.source.task, self.index_task, self.index_task.ld_scores_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        indexed = load_indexed_ldsc_context(self.index_task, fetch)
        loaded = _load_trait(self.source, fetch)
        aligned = _align(self.source, loaded.frame, indexed.variants)
        diagnostics = alignment_diagnostics(aligned, loaded.declared_build)
        check_alignment_guards(
            diagnostics, self.source.alignment_guards or self.guards, self.alias
        )
        _assert_n_finite_where_z_finite(aligned, self.alias)
        return _write_aligned(
            scratch_dir, indexed.variants[GWASLAB_RSID_COL], aligned, diagnostics
        )

    @classmethod
    def create(
        cls,
        asset_id: str,
        *,
        source: SumstatsSource,
        index_task: FastCtLdscVariantIndexTask,
        guards: AlignmentGuards,
    ) -> Self:
        src_meta = source.task.meta
        if not isinstance(src_meta, GWASLabSumStatsMeta):
            assert produces_dataframe(source.task), (
                f"{source.alias}: source must be a gwaslab pickle or a dataframe asset"
            )
            assert isinstance(
                src_meta, FilteredGWASDataMeta | GWASSummaryDataFileMeta
            ), (
                f"{source.alias}: dataframe source meta must be FilteredGWASDataMeta or "
                f"GWASSummaryDataFileMeta, got {type(src_meta).__name__}"
            )
        meta = ResultDirectoryMeta(
            id=asset_id,
            trait=src_meta.trait,
            project=src_meta.project,
            sub_dir=PurePath("analysis/fast_ct_ldsc_aligned"),
        )
        return cls(meta=meta, source=source, index_task=index_task, guards=guards)


def _load_trait(source: SumstatsSource, fetch: Fetch) -> LoadedTrait:
    meta = source.task.meta
    asset = fetch(source.task.asset_id)
    if isinstance(meta, GWASLabSumStatsMeta):
        sumstats = read_sumstats(asset)
        frame = pl.from_pandas(sumstats.data).with_columns(
            pl.col(pl.Categorical, pl.Enum).cast(pl.String)
        )
        return LoadedTrait(
            frame=narwhals.from_native(frame.lazy()),
            declared_build=_pickle_declared_build(
                meta_build=meta.genome_build,
                pickle_build=str(sumstats.build),
                alias=source.alias,
            ),
        )
    assert isinstance(meta, FilteredGWASDataMeta | GWASSummaryDataFileMeta)
    return LoadedTrait(
        frame=scan_dataframe_asset(asset, meta=meta, parquet_backend="polars"),
        declared_build=meta.genome_build,
    )


def _pickle_declared_build(
    *, meta_build: GenomeBuild | None, pickle_build: str, alias: str
) -> GenomeBuild | None:
    recorded: GenomeBuild | None = None
    if pickle_build == "19":
        recorded = "19"
    elif pickle_build == "38":
        recorded = "38"
    if meta_build is None:
        return recorded
    assert recorded is None or recorded == meta_build, (
        f"{alias}: meta declares genome build {meta_build} but the gwaslab pickle records "
        f"build {pickle_build}"
    )
    return meta_build


def _align(
    source: SumstatsSource, frame: narwhals.LazyFrame, variants: pl.DataFrame
) -> TraitAligned:
    processed = source.pipe.process(frame)
    present = set(processed.collect_schema().names())
    missing = [c for c in _REQUIRED_TRAIT_COLUMNS if c not in present]
    assert not missing, f"{source.alias}: trait lacks required columns {missing}"
    columns = [*_REQUIRED_TRAIT_COLUMNS] + [
        c for c in _OPTIONAL_TRAIT_COLUMNS if c in present
    ]
    return align_trait_to_context(
        processed.select(columns),
        variants,
        trait_total_sample_size=source.sample_info.total_sample_size,
        min_trait_snps=0,
    )


def _assert_n_finite_where_z_finite(aligned: TraitAligned, alias: str) -> None:
    bad = int((np.isfinite(aligned.z) & ~np.isfinite(aligned.n)).sum())
    assert bad == 0, (
        f"{alias}: {bad} aligned variants have a finite z but a missing or non-finite N"
    )


def _write_aligned(
    scratch_dir: Path,
    rsid: pl.Series,
    aligned: TraitAligned,
    diagnostics: AlignmentDiagnostics,
) -> DirectoryAsset:
    out_dir = scratch_dir / "aligned_trait"
    out_dir.mkdir()
    pl.DataFrame(
        {
            GWASLAB_RSID_COL: rsid,
            FAST_CT_LDSC_Z_COL: aligned.z,
            FAST_CT_LDSC_N_COL: aligned.n,
        }
    ).write_parquet(out_dir / FAST_CT_LDSC_ALIGNED_FILENAME)
    (out_dir / FAST_CT_LDSC_DIAGNOSTICS_FILENAME).write_text(
        json.dumps(diagnostics.to_json_dict(), indent=2)
    )
    return DirectoryAsset(out_dir)


def read_aligned_trait(
    task: FastCtLdscAlignedTraitTask, fetch: Fetch
) -> AlignedTraitArrays:
    asset = fetch(task.asset_id)
    assert isinstance(asset, DirectoryAsset)
    frame = pl.read_parquet(asset.path / FAST_CT_LDSC_ALIGNED_FILENAME)
    diagnostics = AlignmentDiagnostics.from_json_dict(
        json.loads((asset.path / FAST_CT_LDSC_DIAGNOSTICS_FILENAME).read_text())
    )
    return AlignedTraitArrays(
        rsid=frame[GWASLAB_RSID_COL],
        z=frame[FAST_CT_LDSC_Z_COL].to_numpy().astype(float),
        n=frame[FAST_CT_LDSC_N_COL].to_numpy().astype(float),
        diagnostics=diagnostics,
    )
