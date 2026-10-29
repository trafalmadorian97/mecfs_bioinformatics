"""
Colocalize one trait GWAS with every proteomics GWAS that has a hit at one of its loci.

The algorithm, for one trait:
1. Find the trait's genome-wide significant loci (windows around lead variants); loci in the
   extended MHC are reported but excluded.
2. Pair each locus with every protein whose hit segment (from the protein hit table) overlaps
   it, and select pairs by score within the pair budget (see candidate_selection).
3. For each selected protein, fetch its dense summary statistics for all its selected windows
   in one download, orient each window to GRCh38 with the untrusted genome-reference path,
   align it exactly to the trait window, and run the colocalization method.
4. Write the raw results table and the aligned window of every tested pair.

The build graph is static and which proteins to fetch depends on the hit table's contents, so
all dense fetching happens inside this one Task, bounded by the pair budget.

The trait must be genome-reference harmonized to GRCh38 (NEA = reference allele). Proteins are
modelled as quantitative traits with sdY = 1 (UKB-PPP protein levels are inverse-rank
normalized). A pair with fewer than min_aligned_variants usable aligned variants gets status
failed with a reason instead of posteriors.

Output directory:
- raw_results.parquet: one row per candidate pair (tested, failed or not selected) plus one row
  per excluded locus, integers as Int32.
- windows/{pair}.parquet: the aligned trait x protein window of every tested pair.
"""

from collections.abc import Mapping, Sequence
from functools import partial
from pathlib import Path
from typing import get_args

import numpy as np
import polars as pl
import structlog
from attrs import frozen
from scipy.stats import norm

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    load_fasta,
    load_panel_path,
    resolve_harmonized_build,
    resolve_panel_af_col,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    DropReason,
)
from mecfs_bio.build_system.task.ppp_ldsc.gene_coords import read_gene_coords
from mecfs_bio.build_system.task.protein_coloc.candidate_selection import (
    SelectionConfig,
    build_candidates,
    select_candidates,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import (
    ColocAbf,
    ColocInput,
    ColocMethod,
    ColocPosteriors,
    QuantitativeTraitModel,
    TraitModel,
)
from mecfs_bio.build_system.task.protein_coloc.orient_and_align import (
    OrientedWindow,
    align_windows,
    orient_window,
)
from mecfs_bio.build_system.task.protein_coloc.ppp_window_fetcher import (
    DenseWindowFetcher,
    LocusWindow,
)
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    PppProteinHitTableTask,
    two_sided_mlog10p,
)
from mecfs_bio.build_system.task.protein_coloc.trait_loci import find_trait_loci
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.allele_frequency_panel_constants import PanelAncestry
from mecfs_bio.constants.genomic_coordinate_constants import extended_mhc_interval
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_POS_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    BETA_PROTEIN_COL,
    BETA_TRAIT_COL,
    EAF_PROTEIN_COL,
    EAF_TRAIT_COL,
    FAILURE_REASON_COL,
    LOCUS_CHR_COL,
    LOCUS_END_COL,
    LOCUS_ID_COL,
    LOCUS_START_COL,
    METHOD_COL,
    N_ALIGNED_COL,
    N_DROPPED_PREFIX,
    N_DUPLICATE_KEYS_DROPPED_COL,
    N_PROTEIN_COL,
    N_PROTEIN_VARIANTS_COL,
    N_TRAIT_COL,
    N_TRAIT_VARIANTS_COL,
    PP_H0_COL,
    PP_H1_COL,
    PP_H2_COL,
    PP_H3_COL,
    PP_H4_COL,
    PRIOR_P1_COL,
    PRIOR_P2_COL,
    PRIOR_P12_COL,
    PROTEIN_ID_COL,
    PROTEIN_LEAD_MLOG10P_COL,
    PROTEIN_LEAD_POS_COL,
    RAW_RESULTS_FILENAME,
    SE_PROTEIN_COL,
    SE_TRAIT_COL,
    STATUS_COL,
    STATUS_EXCLUDED_MHC,
    STATUS_FAILED,
    STATUS_SELECTED,
    STATUS_TESTED,
    WINDOW_SCHEMA,
    WINDOWS_DIRNAME,
    pair_file_stem,
)
from mecfs_bio.util.retry.call_with_retries import call_with_retries

logger = structlog.get_logger()

# UKB-PPP protein levels are inverse-rank normalized, so their phenotypic SD is 1.
PROTEIN_TRAIT_MODEL = QuantitativeTraitModel(sd_y=1.0)
# Dense fetches hit remote services; network and filesystem errors are worth retrying.
_FETCH_TRANSIENT_ERRORS: tuple[type[Exception], ...] = (OSError,)
_RESULT_STATUS_COL = "_result_status"
DROP_REASONS: tuple[str, ...] = get_args(DropReason)


def drop_count_col(reason: str) -> str:
    return f"{N_DROPPED_PREFIX}{reason.upper()}"


@frozen(slots=True)
class ColocRunConfig:
    """
    lead_p_threshold: a trait variant can lead a locus when its p is below this.
    half_window_bp: a locus window is lead +- half_window_bp (then overlapping windows merge).
    exclude_mhc: report, but do not test, loci overlapping the GRCh38 extended MHC.
    selection: the pair budget and scoring.
    min_aligned_variants: a pair with fewer usable aligned variants fails instead of being tested.
    harmonization: options of the untrusted orientation of protein windows.
    """

    lead_p_threshold: float = 5e-8
    half_window_bp: int = 500_000
    exclude_mhc: bool = True
    selection: SelectionConfig = SelectionConfig()
    min_aligned_variants: int = 50
    harmonization: GenomeReferenceHarmonizationOptions = (
        GenomeReferenceHarmonizationOptions()
    )

    def __attrs_post_init__(self) -> None:
        assert 0 < self.lead_p_threshold < 1
        assert self.half_window_bp > 0
        assert self.min_aligned_variants >= 1


@frozen(slots=True)
class _OrientationReference:
    fasta: IndexedFasta
    panel_path: Path
    panel_af_col: str
    options: GenomeReferenceHarmonizationOptions


@frozen(slots=True)
class _PairOutcome:
    """The result of testing one selected (locus, protein) pair."""

    locus_id: str
    protein_id: str
    oriented: OrientedWindow
    n_trait: int
    n_aligned: int
    n_duplicate_keys_dropped: int
    posteriors: ColocPosteriors | None
    failure_reason: str | None
    protein_lead_pos: int | None
    protein_lead_mlog10p: float | None


def _scan_trait(trait_task: Task, fetch: Fetch) -> pl.LazyFrame:
    """The trait in WINDOW_SCHEMA; EAF and N are null when the trait has none."""
    lazy = scan_dataframe_asset(
        fetch(trait_task.asset_id), meta=trait_task.meta, parquet_backend="polars"
    ).to_native()
    assert isinstance(lazy, pl.LazyFrame)
    names = lazy.collect_schema().names()
    return lazy.select(
        [
            pl.col(name) if name in names else pl.lit(None).alias(name)
            for name in WINDOW_SCHEMA.names()
        ]
    ).cast(WINDOW_SCHEMA)


def _assert_trait_supports_model(trait: pl.LazyFrame, model: TraitModel) -> None:
    if isinstance(model, QuantitativeTraitModel) and model.sd_y is None:
        n_usable = (
            trait.select(
                (
                    pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL).is_not_null()
                    & pl.col(GWASLAB_SAMPLE_SIZE_COLUMN).is_not_null()
                ).sum()
            )
            .collect()
            .item()
        )
        assert n_usable > 0, (
            "a quantitative trait without sd_y needs EAF and N to estimate sdY"
        )


def _significant_trait_rows(trait: pl.LazyFrame, p_threshold: float) -> pl.DataFrame:
    """A superset of the variants below p_threshold, filtered lazily to bound memory."""
    z_threshold = float(norm.isf(p_threshold / 2)) * (1 - 1e-6)
    return trait.filter(
        (pl.col(GWASLAB_BETA_COL) / pl.col(GWASLAB_SE_COL)).abs() >= z_threshold
    ).collect()


def _read_table(task: Task, fetch: Fetch) -> pl.DataFrame:
    lazy = scan_dataframe_asset(
        fetch(task.asset_id), meta=task.meta, parquet_backend="polars"
    ).to_native()
    assert isinstance(lazy, pl.LazyFrame)
    return lazy.collect()


def _locus_window(loci: pl.DataFrame, locus_id: str) -> LocusWindow:
    row = loci.filter(pl.col(LOCUS_ID_COL) == locus_id).row(0, named=True)
    return LocusWindow(
        locus_id=locus_id,
        chrom=row[LOCUS_CHR_COL],
        start=row[LOCUS_START_COL],
        end=row[LOCUS_END_COL],
    )


def _trait_window(trait: pl.LazyFrame, window: LocusWindow) -> pl.DataFrame:
    return trait.filter(
        pl.col(GWASLAB_CHROM_COL) == window.chrom,
        pl.col(GWASLAB_POS_COL).is_between(window.start, window.end),
    ).collect()


def _usable_variants(aligned: pl.DataFrame) -> pl.DataFrame:
    """Aligned variants with finite effects and positive standard errors on both sides."""
    return aligned.filter(
        pl.col(BETA_TRAIT_COL).is_finite(),
        pl.col(BETA_PROTEIN_COL).is_finite(),
        pl.col(SE_TRAIT_COL).is_finite() & (pl.col(SE_TRAIT_COL) > 0),
        pl.col(SE_PROTEIN_COL).is_finite() & (pl.col(SE_PROTEIN_COL) > 0),
    )


def _optional_array(frame: pl.DataFrame, column: str) -> np.ndarray | None:
    values = frame[column].cast(pl.Float64).to_numpy()
    return values if np.isfinite(values).any() else None


def _coloc_input(window: pl.DataFrame, trait_model: TraitModel) -> ColocInput:
    """EAF stands in for MAF: coloc only uses f * (1 - f), which is symmetric."""
    return ColocInput(
        beta1=window[BETA_TRAIT_COL].to_numpy().astype(np.float64),
        se1=window[SE_TRAIT_COL].to_numpy().astype(np.float64),
        beta2=window[BETA_PROTEIN_COL].to_numpy().astype(np.float64),
        se2=window[SE_PROTEIN_COL].to_numpy().astype(np.float64),
        maf1=_optional_array(window, EAF_TRAIT_COL),
        n1=_optional_array(window, N_TRAIT_COL),
        maf2=_optional_array(window, EAF_PROTEIN_COL),
        n2=_optional_array(window, N_PROTEIN_COL),
        model1=trait_model,
        model2=PROTEIN_TRAIT_MODEL,
    )


def _orient(
    rows: pl.DataFrame, window: LocusWindow, reference: _OrientationReference
) -> OrientedWindow:
    if rows.height == 0:
        return OrientedWindow(frame=WINDOW_SCHEMA.to_frame(), n_input=0, drop_counts={})
    return orient_window(
        rows,
        chrom=window.chrom,
        fasta=reference.fasta,
        panel_path=reference.panel_path,
        panel_af_col=reference.panel_af_col,
        options=reference.options,
    )


@frozen(slots=True)
class _ProteinLead:
    pos: int | None
    mlog10p: float | None


def _protein_lead(window: pl.DataFrame) -> _ProteinLead:
    """The most significant protein variant of a window."""
    if window.height == 0:
        return _ProteinLead(pos=None, mlog10p=None)
    with np.errstate(invalid="ignore", divide="ignore"):
        mlog10p = two_sided_mlog10p(
            (window[BETA_PROTEIN_COL] / window[SE_PROTEIN_COL]).to_numpy()
        )
    best = int(np.nanargmax(mlog10p))
    return _ProteinLead(
        pos=int(window[GWASLAB_POS_COL][best]), mlog10p=float(mlog10p[best])
    )


def _test_pair(
    window: LocusWindow,
    protein_id: str,
    protein_rows: pl.DataFrame,
    trait_rows: pl.DataFrame,
    trait_model: TraitModel,
    method: ColocMethod,
    reference: _OrientationReference,
    min_aligned_variants: int,
    windows_dir: Path,
) -> _PairOutcome:
    oriented = _orient(protein_rows, window, reference)
    aligned = align_windows(trait=trait_rows, protein=oriented.frame)
    usable = _usable_variants(aligned.frame)
    lead = _protein_lead(usable)
    enough = usable.height >= min_aligned_variants
    if enough:
        usable.write_parquet(
            windows_dir / f"{pair_file_stem(window.locus_id, protein_id)}.parquet"
        )
    return _PairOutcome(
        locus_id=window.locus_id,
        protein_id=protein_id,
        oriented=oriented,
        n_trait=aligned.n_trait,
        n_aligned=usable.height,
        n_duplicate_keys_dropped=aligned.n_duplicate_keys_dropped,
        posteriors=method.run(_coloc_input(usable, trait_model)) if enough else None,
        failure_reason=(
            None
            if enough
            else f"{usable.height} usable aligned variants, fewer than {min_aligned_variants}"
        ),
        protein_lead_pos=lead.pos,
        protein_lead_mlog10p=lead.mlog10p,
    )


def _test_protein(
    protein_id: str,
    windows: Sequence[LocusWindow],
    trait: pl.LazyFrame,
    task: "TraitProteinColocTask",
    reference: _OrientationReference,
    scratch_dir: Path,
    windows_dir: Path,
    wf: WF,
) -> list[_PairOutcome]:
    """One dense fetch for all of a protein's selected windows, then test each pair."""
    frames = call_with_retries(
        partial(task.fetcher.fetch, protein_id, windows, scratch_dir, wf),
        retry_on=_FETCH_TRANSIENT_ERRORS,
    )
    logger.info(
        "fetched dense protein windows", protein_id=protein_id, n_windows=len(windows)
    )
    return [
        _test_pair(
            window=window,
            protein_id=protein_id,
            protein_rows=frames.get(window.locus_id, WINDOW_SCHEMA.to_frame()),
            trait_rows=_trait_window(trait, window),
            trait_model=task.trait_model,
            method=task.method,
            reference=reference,
            min_aligned_variants=task.config.min_aligned_variants,
            windows_dir=windows_dir,
        )
        for window in windows
    ]


def _windows_by_protein(
    ranked: pl.DataFrame, loci: pl.DataFrame
) -> dict[str, list[LocusWindow]]:
    selected = ranked.filter(pl.col(STATUS_COL) == STATUS_SELECTED)
    by_protein: dict[str, list[LocusWindow]] = {}
    for protein_id, locus_id in selected.select(PROTEIN_ID_COL, LOCUS_ID_COL).rows():
        by_protein.setdefault(protein_id, []).append(_locus_window(loci, locus_id))
    return by_protein


def _test_selected_pairs(
    ranked: pl.DataFrame,
    loci: pl.DataFrame,
    trait: pl.LazyFrame,
    task: "TraitProteinColocTask",
    reference: _OrientationReference,
    scratch_dir: Path,
    windows_dir: Path,
    wf: WF,
) -> list[_PairOutcome]:
    """Test every selected pair, fetching each selected protein once."""
    return [
        outcome
        for protein_id, windows in _windows_by_protein(ranked, loci).items()
        for outcome in _test_protein(
            protein_id, windows, trait, task, reference, scratch_dir, windows_dir, wf
        )
    ]


def _outcome_row(outcome: _PairOutcome, method: ColocMethod) -> dict[str, object]:
    post = outcome.posteriors
    row: dict[str, object] = {
        LOCUS_ID_COL: outcome.locus_id,
        PROTEIN_ID_COL: outcome.protein_id,
        _RESULT_STATUS_COL: STATUS_TESTED if post is not None else STATUS_FAILED,
        FAILURE_REASON_COL: outcome.failure_reason,
        N_TRAIT_VARIANTS_COL: outcome.n_trait,
        N_PROTEIN_VARIANTS_COL: outcome.oriented.n_input,
        N_ALIGNED_COL: outcome.n_aligned,
        N_DUPLICATE_KEYS_DROPPED_COL: outcome.n_duplicate_keys_dropped,
        PROTEIN_LEAD_POS_COL: outcome.protein_lead_pos,
        PROTEIN_LEAD_MLOG10P_COL: outcome.protein_lead_mlog10p,
        PP_H0_COL: post.pp_h0 if post is not None else None,
        PP_H1_COL: post.pp_h1 if post is not None else None,
        PP_H2_COL: post.pp_h2 if post is not None else None,
        PP_H3_COL: post.pp_h3 if post is not None else None,
        PP_H4_COL: post.pp_h4 if post is not None else None,
        PRIOR_P1_COL: method.priors.p1,
        PRIOR_P2_COL: method.priors.p2,
        PRIOR_P12_COL: method.priors.p12,
        METHOD_COL: method.name,
    }
    for reason in DROP_REASONS:
        row[drop_count_col(reason)] = outcome.oriented.drop_counts.get(reason, 0)
    return row


_OUTCOME_SCHEMA = pl.Schema(
    {
        LOCUS_ID_COL: pl.String,
        PROTEIN_ID_COL: pl.String,
        _RESULT_STATUS_COL: pl.String,
        FAILURE_REASON_COL: pl.String,
        N_TRAIT_VARIANTS_COL: pl.Int32,
        N_PROTEIN_VARIANTS_COL: pl.Int32,
        N_ALIGNED_COL: pl.Int32,
        N_DUPLICATE_KEYS_DROPPED_COL: pl.Int32,
        PROTEIN_LEAD_POS_COL: pl.Int32,
        PROTEIN_LEAD_MLOG10P_COL: pl.Float64,
        PP_H0_COL: pl.Float64,
        PP_H1_COL: pl.Float64,
        PP_H2_COL: pl.Float64,
        PP_H3_COL: pl.Float64,
        PP_H4_COL: pl.Float64,
        PRIOR_P1_COL: pl.Float64,
        PRIOR_P2_COL: pl.Float64,
        PRIOR_P12_COL: pl.Float64,
        METHOD_COL: pl.String,
        **{drop_count_col(reason): pl.Int32 for reason in DROP_REASONS},
    }
)


def _assemble_results(
    loci: pl.DataFrame,
    ranked: pl.DataFrame,
    outcomes: Sequence[_PairOutcome],
    method: ColocMethod,
) -> pl.DataFrame:
    """Candidates with their test outcomes and locus fields, plus the excluded loci."""
    outcome_table = pl.DataFrame(
        [_outcome_row(outcome, method) for outcome in outcomes], schema=_OUTCOME_SCHEMA
    )
    locus_fields = loci.drop(STATUS_COL)
    pairs = (
        ranked.join(outcome_table, on=[LOCUS_ID_COL, PROTEIN_ID_COL], how="left")
        .with_columns(pl.coalesce(_RESULT_STATUS_COL, STATUS_COL).alias(STATUS_COL))
        .drop(_RESULT_STATUS_COL)
        .join(locus_fields, on=LOCUS_ID_COL, how="left")
    )
    excluded = loci.filter(pl.col(STATUS_COL) == STATUS_EXCLUDED_MHC)
    results = pl.concat([pairs, excluded], how="diagonal_relaxed")
    return _int_columns_as_int32(results)


def _int_columns_as_int32(frame: pl.DataFrame) -> pl.DataFrame:
    """Tabulator cannot format parquet INT64, so every integer column is stored as Int32."""
    return frame.with_columns(
        [
            pl.col(name).cast(pl.Int32)
            for name, dtype in frame.schema.items()
            if dtype.is_integer()
        ]
    )


@frozen(slots=True)
class TraitProteinColocTask(Task):
    """
    trait_task: the GRCh38 genome-reference harmonized trait GWAS.
    trait_model: quantitative (sdY known or estimated) or case-control.
    hit_table_task: the protein hit-segment table that screens candidate proteins.
    gene_coords_task: the extracted Sun et al. 2023 ST3 sheet (GRCh38 gene coordinates).
    fasta_task, panel_task, panel_af_col: GRCh38 reference genome and allele-frequency panel
        used to orient protein windows.
    fetcher: the dense window source of the proteomics GWAS.
    method: the colocalization method.
    """

    meta: ResultDirectoryMeta
    trait_task: Task
    trait_model: TraitModel
    hit_table_task: PppProteinHitTableTask
    gene_coords_task: Task
    fasta_task: Task
    panel_task: Task
    panel_af_col: str
    fetcher: DenseWindowFetcher
    method: ColocMethod
    config: ColocRunConfig

    @property
    def deps(self) -> list[Task]:
        return [
            self.trait_task,
            self.hit_table_task,
            self.gene_coords_task,
            self.fasta_task,
            self.panel_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        trait = _scan_trait(self.trait_task, fetch)
        _assert_trait_supports_model(trait, self.trait_model)
        loci = find_trait_loci(
            _significant_trait_rows(trait, self.config.lead_p_threshold),
            p_threshold=self.config.lead_p_threshold,
            half_window_bp=self.config.half_window_bp,
            exclude=[extended_mhc_interval("38")] if self.config.exclude_mhc else [],
        )
        candidates = build_candidates(
            loci,
            _read_table(self.hit_table_task, fetch),
            gene_coords=read_gene_coords(_read_table(self.gene_coords_task, fetch)),
            cis_window_bp=self.config.selection.cis_window_bp,
        )
        ranked = select_candidates(candidates, loci, self.config.selection)
        reference = _OrientationReference(
            fasta=load_fasta(fetch, self.fasta_task),
            panel_path=load_panel_path(fetch, self.panel_task),
            panel_af_col=self.panel_af_col,
            options=self.config.harmonization,
        )
        out_dir = scratch_dir / "coloc"
        windows_dir = out_dir / WINDOWS_DIRNAME
        windows_dir.mkdir(parents=True)
        fetch_dir = scratch_dir / "dense"
        fetch_dir.mkdir()
        outcomes = _test_selected_pairs(
            ranked, loci, trait, self, reference, fetch_dir, windows_dir, wf
        )
        results = _assemble_results(loci, ranked, outcomes, self.method)
        _log_summary(loci, results)
        results.write_parquet(out_dir / RAW_RESULTS_FILENAME)
        return DirectoryAsset(out_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        trait_task: Task,
        trait_model: TraitModel,
        hit_table_task: PppProteinHitTableTask,
        gene_coords_task: Task,
        fasta_task: Task,
        panel_task: Task,
        panel_ancestry: PanelAncestry,
        fetcher: DenseWindowFetcher,
        method: ColocMethod = ColocAbf(),
        config: ColocRunConfig = ColocRunConfig(),
    ) -> "TraitProteinColocTask":
        trait_meta = trait_task.meta
        assert isinstance(trait_meta, FilteredGWASDataMeta), (
            f"trait_task must carry FilteredGWASDataMeta, got {type(trait_meta).__name__}"
        )
        assert (
            trait_meta.harmonization_info is not None
            and trait_meta.harmonization_info.build == "38"
        ), "the trait must be genome-reference harmonized to GRCh38"
        assert resolve_harmonized_build(fasta_task, panel_task) == "38", (
            "the FASTA and panel must be GRCh38"
        )
        return cls(
            meta=ResultDirectoryMeta(
                id=asset_id,
                trait=trait_meta.trait,
                project=trait_meta.project,
            ),
            trait_task=trait_task,
            trait_model=trait_model,
            hit_table_task=hit_table_task,
            gene_coords_task=gene_coords_task,
            fasta_task=fasta_task,
            panel_task=panel_task,
            panel_af_col=resolve_panel_af_col(panel_task, panel_ancestry),
            fetcher=fetcher,
            method=method,
            config=config,
        )


def _log_summary(loci: pl.DataFrame, results: pl.DataFrame) -> None:
    counts: Mapping[str, int] = dict(results.group_by(STATUS_COL).len().rows())
    logger.info(
        "trait x protein colocalization summary",
        n_loci=loci.height,
        rows_by_status=counts,
    )
