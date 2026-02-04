from pathlib import Path

import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

from .harmonize_gwas_with_reference_table_via_rsid import (
    IS_PALINDROMIC,
    MATCH_REFERENCE,
    MATCH_REFERENCE_FLIPPED,
    REFERENCE_CHROM,
    REFERENCE_EFFECT_ALLELE,
    REFERENCE_NON_EFFECT_ALLELE,
    REFERENCE_POS,
    PalindromeStrategy,
    handle_flipped,
    is_palindromic_expr,
)

logger = structlog.get_logger()


_REF_COLS_ALLELE_MATCH = [
    IS_PALINDROMIC,
    MATCH_REFERENCE,
    MATCH_REFERENCE_FLIPPED,
]


@frozen
class HarmonizeGWASWithReferenceViaAlleles(Task):
    """
    Given a table of reference genetic variants, harmonize gwas data with that table of reference variants
    using chromsome, position, alleles for matching

    In this context harmonization means:
    For non-palindromic variants:
       - keep only variants where:
          -- effect and non-effect allele match the reference
          --effect allele matches reference non-effect allele, and vice versa.  In this case, flip beta and invert odds ratio

    For palindromic variants:
    - drop if we lack frequency information
    - Else, can try to use frequency info to resolve


    """

    _meta: Meta
    gwas_data_task: Task
    reference_task: Task
    palindrome_strategy: PalindromeStrategy
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

        reference_asset = fetch(self.reference_asset_id)
        reference = (
            self.ref_pipe.process(
                scan_dataframe_asset(reference_asset, meta=self.reference_meta)
            )
            .collect()
            .to_polars()
        )

        assert len(
            gwas_data.unique(
                subset=[
                    GWASLAB_CHROM_COL,
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL,
                ]
            )
        ) == len(gwas_data)

        assert len(
            reference.unique(
                subset=[
                    GWASLAB_CHROM_COL,
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL,
                ]
            )
        ) == len(reference)

        reference = reference.select(
            [
                GWASLAB_CHROM_COL,
                GWASLAB_POS_COL,
                GWASLAB_EFFECT_ALLELE_COL,
                GWASLAB_NON_EFFECT_ALLELE_COL,
            ]
        ).rename(
            {
                GWASLAB_EFFECT_ALLELE_COL: REFERENCE_EFFECT_ALLELE,
                GWASLAB_NON_EFFECT_ALLELE_COL: REFERENCE_NON_EFFECT_ALLELE,
                GWASLAB_POS_COL: REFERENCE_POS,
                GWASLAB_CHROM_COL: REFERENCE_CHROM,
            }
        )

        gd = gwas_data.with_columns(
            is_palindromic_expr(
                GWASLAB_EFFECT_ALLELE_COL, nea_col=GWASLAB_NON_EFFECT_ALLELE_COL
            ).alias(IS_PALINDROMIC)
        )

        if self.palindrome_strategy == "drop":
            logger.debug(f"dropping {gd[IS_PALINDROMIC].sum()} palindromic variants")
            gd = gd.filter(~pl.col(IS_PALINDROMIC))
        else:
            raise NotImplementedError(
                f"mode {self.palindrome_strategy} not implemented"
            )

        gd_match = gd.join(
            reference,
            left_on=[
                GWASLAB_CHROM_COL,
                GWASLAB_POS_COL,
                GWASLAB_EFFECT_ALLELE_COL,
                GWASLAB_NON_EFFECT_ALLELE_COL,
            ],
            right_on=[
                REFERENCE_CHROM,
                REFERENCE_POS,
                REFERENCE_EFFECT_ALLELE,
                REFERENCE_NON_EFFECT_ALLELE,
            ],
        ).with_columns(
            pl.lit(False).alias(MATCH_REFERENCE_FLIPPED),
            pl.lit(True).alias(MATCH_REFERENCE),
        )
        gd_reverse_match = gd.join(
            reference,
            left_on=[
                GWASLAB_CHROM_COL,
                GWASLAB_POS_COL,
                GWASLAB_EFFECT_ALLELE_COL,
                GWASLAB_NON_EFFECT_ALLELE_COL,
            ],
            right_on=[
                REFERENCE_CHROM,
                REFERENCE_POS,
                REFERENCE_NON_EFFECT_ALLELE,
                REFERENCE_EFFECT_ALLELE,
            ],
        ).with_columns(
            pl.lit(True).alias(MATCH_REFERENCE_FLIPPED),
            pl.lit(False).alias(MATCH_REFERENCE),
        )
        gd = pl.concat([gd_match, gd_reverse_match], how="vertical")
        gd = handle_flipped(gd)
        _report_matches(gd)
        gd = gd.drop(_REF_COLS_ALLELE_MATCH)
        out_path = scratch_dir / "matched_to_ref.parquet"
        gd.write_parquet(out_path)

        return FileAsset(out_path)


def _report_matches(df: pl.DataFrame) -> None:
    num_matches = df[MATCH_REFERENCE].sum()
    flip_matches = df[MATCH_REFERENCE_FLIPPED].sum()
    logger.debug(
        f"Found {num_matches} variants matching the reference in their base orientation"
        f"Found {flip_matches} variants matching the reference when effect and non-effect alleles are flipped"
    )
