from pathlib import Path
from typing import Literal

import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_ODDS_RATIO_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)

logger = structlog.get_logger()

PalindromeStrategy = Literal["drop", "resolve"]

ChromPosStrategy = Literal["drop", "raise"]

MismatchedAllelesStrategy = Literal["drop", "raise"]

IS_PALINDROMIC = "__is_palindromic__"

REFERENCE_EFFECT_ALLELE = "__reference_effect_allele__"
REFERENCE_NON_EFFECT_ALLELE = "__reference_non_effect_allele__"
REFERENCE_POS = "__reference_pos__"
REFERENCE_CHROM = "__reference_chrom__"

_MISMATCH_POS_CHROM = "__MISMATCH_POS_CHROM__"
MATCH_REFERENCE = "__MATCH_REFERENCE__"
MATCH_REFERENCE_FLIPPED = "__MATCH_REFERENCE_FLIPPED__"


_REF_COLS = [
    REFERENCE_EFFECT_ALLELE,
    REFERENCE_NON_EFFECT_ALLELE,
    REFERENCE_POS,
    REFERENCE_CHROM,
    MATCH_REFERENCE,
    MATCH_REFERENCE_FLIPPED,
    _MISMATCH_POS_CHROM,
]


@frozen
class HarmonizeGWASWithReferenceViaRSIDTask(Task):
    """
    Given a table of reference genetic variants, harmonize gwas data with that table of reference variants using rsid for matching

    In this context harmonization means:
    For non-palindromic variants:
       - keep only variants where:
          -- effect and non-effect allele match the reference, or
          --effect allele matches reference non-effect allele, and vice versa.  In this case, flip beta and invert odds ratio

    For palindromic variants:
    - drop if we lack frequency information
    - Else, can try to use frequency info to resolve


    """

    _meta: Meta
    gwas_data_task: Task
    reference_task: Task
    palindrome_strategy: PalindromeStrategy
    chrom_pos_strategy: ChromPosStrategy
    mismatched_alleles_strategy: MismatchedAllelesStrategy
    gwas_pipe: DataProcessingPipe = IdentityPipe()
    ref_pipe: DataProcessingPipe = IdentityPipe()

    @property
    def source_gwas_asset_id(self) -> AssetId:
        return self.gwas_data_task.asset_id

    @property
    def source_gwas_meta(self) -> Meta:
        return self.gwas_data_task.meta

    @property
    def reference_asset_id(self) -> AssetId:
        return self.reference_task.asset_id

    @property
    def reference_meta(self) -> Meta:
        return self.reference_task.meta

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.gwas_data_task, self.reference_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        gwas_asset = fetch(self.source_gwas_asset_id)
        gwas_data = (
            self.gwas_pipe.process(
                scan_dataframe_asset(gwas_asset, meta=self.source_gwas_meta)
            )
            .collect()
            .to_polars()
        )
        assert GWASLAB_RSID_COL in gwas_data.columns

        reference_asset = fetch(self.reference_asset_id)
        reference = (
            self.ref_pipe.process(
                scan_dataframe_asset(reference_asset, meta=self.reference_meta)
            )
            .collect()
            .to_polars()
        )
        assert GWASLAB_RSID_COL in reference.columns
        reference = reference.rename(
            {
                GWASLAB_EFFECT_ALLELE_COL: REFERENCE_EFFECT_ALLELE,
                GWASLAB_NON_EFFECT_ALLELE_COL: REFERENCE_NON_EFFECT_ALLELE,
                GWASLAB_POS_COL: REFERENCE_POS,
                GWASLAB_CHROM_COL: REFERENCE_CHROM,
            }
        )

        assert (
            reference[GWASLAB_RSID_COL].n_unique() == reference[GWASLAB_RSID_COL].len()
        )
        assert (
            gwas_data[GWASLAB_RSID_COL].n_unique() == gwas_data[GWASLAB_RSID_COL].len()
        )

        gd = gwas_data.with_columns(
            is_palindromic_expr(
                GWASLAB_EFFECT_ALLELE_COL, nea_col=GWASLAB_NON_EFFECT_ALLELE_COL
            ).alias(IS_PALINDROMIC)
        )

        gd = gd.join(reference, on=GWASLAB_RSID_COL)
        if self.palindrome_strategy == "drop":
            logger.debug(f"Dropping {gd[IS_PALINDROMIC].sum()} palindromic variants\n")
            gd = gd.filter(~pl.col(IS_PALINDROMIC))
        else:
            raise NotImplementedError(
                f"mode {self.palindrome_strategy} not implemented"
            )
        gd = gd.with_columns(
            _mismatch_pos_or_chrom_expr(
                pos_col=GWASLAB_POS_COL,
                chrom_col=GWASLAB_CHROM_COL,
                ref_pos_col=REFERENCE_POS,
                ref_chrom_col=REFERENCE_CHROM,
            ).alias(_MISMATCH_POS_CHROM),
            _match_reference_expr(
                ea_col=GWASLAB_EFFECT_ALLELE_COL,
                nea_col=GWASLAB_NON_EFFECT_ALLELE_COL,
                ref_ea_col=REFERENCE_EFFECT_ALLELE,
                ref_nea_col=REFERENCE_NON_EFFECT_ALLELE,
            ).alias(MATCH_REFERENCE),
            _match_flipped_reference_expr(
                ea_col=GWASLAB_EFFECT_ALLELE_COL,
                nea_col=GWASLAB_NON_EFFECT_ALLELE_COL,
                ref_ea_col=REFERENCE_EFFECT_ALLELE,
                ref_nea_col=REFERENCE_NON_EFFECT_ALLELE,
            ).alias(MATCH_REFERENCE_FLIPPED),
        )

        gd = _handle_chrom_pos(gd, self.chrom_pos_strategy)
        gd = handle_flipped(gd)
        _report_matches(gd)
        gd = handle_mismatched_alleles(gd, strategy=self.mismatched_alleles_strategy)
        gd = gd.drop(_REF_COLS)
        logger.debug(f"harmonized gwas has length: {len(gd)}")
        out_path = scratch_dir / "matched_to_ref.parquet"
        gd.write_parquet(out_path)

        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        gwas_data_task: Task,
        reference_task: Task,
        palindrome_strategy: PalindromeStrategy,
        chrom_pos_strategy: ChromPosStrategy,
        mismatched_alleles_strategy: MismatchedAllelesStrategy,
        gwas_pipe: DataProcessingPipe = IdentityPipe(),
        ref_pipe: DataProcessingPipe = IdentityPipe(),
    ):
        source_meta = gwas_data_task.meta
        read_spec = DataFrameReadSpec(DataFrameParquetFormat())
        meta: Meta
        if isinstance(source_meta, FilteredGWASDataMeta):
            meta = FilteredGWASDataMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                sub_dir=source_meta.sub_dir,
                read_spec=read_spec,
            )
        elif isinstance(source_meta, GWASSummaryDataFileMeta):
            meta = GWASSummaryDataFileMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                sub_dir=source_meta.sub_dir,
                project_path=None,
                read_spec=read_spec,
            )
        return cls(
            meta=meta,
            gwas_data_task=gwas_data_task,
            reference_task=reference_task,
            palindrome_strategy=palindrome_strategy,
            chrom_pos_strategy=chrom_pos_strategy,
            mismatched_alleles_strategy=mismatched_alleles_strategy,
            gwas_pipe=gwas_pipe,
            ref_pipe=ref_pipe,
        )


def complement_reverse_expr(col_name: str) -> pl.Expr:
    return (
        pl.col(col_name)
        .str.replace_all("A", "X")
        .str.replace_all("T", "A")
        .str.replace_all("X", "T")
        .str.replace_all("C", "Y")
        .str.replace_all("G", "C")
        .str.replace_all("Y", "G")
        .str.reverse()
    )


def is_palindromic_expr(ea_col: str, nea_col: str) -> pl.Expr:
    """
    Check if a variant is palindromic.
    """
    return pl.col(ea_col) == complement_reverse_expr(nea_col)


def _mismatch_pos_or_chrom_expr(
    pos_col: str,
    chrom_col: str,
    ref_pos_col: str,
    ref_chrom_col: str,
) -> pl.Expr:
    return (pl.col(pos_col) != pl.col(ref_pos_col)) | (
        pl.col(chrom_col) != pl.col(ref_chrom_col)
    )


def _match_reference_expr(
    ea_col: str,
    nea_col: str,
    ref_ea_col: str,
    ref_nea_col: str,
) -> pl.Expr:
    return (pl.col(ea_col) == pl.col(ref_ea_col)) & (
        pl.col(nea_col) == pl.col(ref_nea_col)
    )


def _match_flipped_reference_expr(
    ea_col: str,
    nea_col: str,
    ref_ea_col: str,
    ref_nea_col: str,
) -> pl.Expr:
    return (pl.col(ea_col) == pl.col(ref_nea_col)) & (
        pl.col(nea_col) == pl.col(ref_ea_col)
    )


def _handle_chrom_pos(df: pl.DataFrame, strategy: ChromPosStrategy) -> pl.DataFrame:
    num_mismatch = df[_MISMATCH_POS_CHROM].sum()
    msg = f"Found {num_mismatch} variants with matching rsid but mismatched chromosome or position"
    if num_mismatch > 0:
        if strategy == "raise":
            raise ValueError(msg)
        logger.debug(msg)
        df = df.filter(~pl.col(_MISMATCH_POS_CHROM))
    return df


def handle_flipped(
    df: pl.DataFrame,
) -> pl.DataFrame:
    df = df.with_columns(
        pl.when(MATCH_REFERENCE_FLIPPED)
        .then(pl.col(GWASLAB_NON_EFFECT_ALLELE_COL))
        .otherwise(pl.col(GWASLAB_EFFECT_ALLELE_COL))
        .alias(GWASLAB_EFFECT_ALLELE_COL),
        pl.when(MATCH_REFERENCE_FLIPPED)
        .then(pl.col(GWASLAB_EFFECT_ALLELE_COL))
        .otherwise(pl.col(GWASLAB_NON_EFFECT_ALLELE_COL))
        .alias(GWASLAB_NON_EFFECT_ALLELE_COL),
        pl.when(MATCH_REFERENCE_FLIPPED)
        .then(-pl.col(GWASLAB_BETA_COL))
        .otherwise(pl.col(GWASLAB_BETA_COL))
        .alias(GWASLAB_BETA_COL),
    )
    if GWASLAB_ODDS_RATIO_COL in df.columns:
        df = df.with_columns(
            pl.when(MATCH_REFERENCE_FLIPPED)
            .then(
                1 / pl.col(GWASLAB_ODDS_RATIO_COL),
            )
            .otherwise(pl.col(GWASLAB_ODDS_RATIO_COL))
            .alias(GWASLAB_ODDS_RATIO_COL),
        )
    if GWASLAB_EFFECT_ALLELE_FREQ_COL in df.columns:
        df = df.with_columns(
            pl.when(MATCH_REFERENCE_FLIPPED)
            .then(1 - pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL))
            .otherwise(pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL))
            .alias(GWASLAB_EFFECT_ALLELE_FREQ_COL),
        )
    return df


def handle_mismatched_alleles(
    df: pl.DataFrame, strategy: MismatchedAllelesStrategy
) -> pl.DataFrame:
    mismatched = (~(df[MATCH_REFERENCE])) & (~(df[MATCH_REFERENCE_FLIPPED]))
    num_mismatch = (mismatched).sum()
    msg = f"Found {num_mismatch} variants with matching rsid but mismatched alleles"
    if num_mismatch > 0:
        if strategy == "raise":
            raise ValueError(msg)
        logger.debug(msg)
        df = df.filter(~mismatched)
    return df


def _report_matches(df: pl.DataFrame) -> None:
    num_matches = df[MATCH_REFERENCE].sum()
    flip_matches = df[MATCH_REFERENCE_FLIPPED].sum()
    logger.debug(
        f"Found {num_matches} variants matching the reference in their base orientation\n\n"
        f"Found {flip_matches} variants matching the reference when effect and non-effect alleles are flipped"
    )
