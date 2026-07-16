"""
Construct the shared variant index used to store UKB-PPP summary statistics
compactly.

The index is a parquet table whose ROW ORDER is the canonical alignment slot for
every per-protein beta/se file: a protein's aligned file stores only beta and se,
one row per index row, in this order. This task builds one index (HapMap3 ~1.20M
rows, or common-1kg ~8.26M rows depending on the membership list supplied).

Construction strategy (settled empirically by the experiments/claude probes): the
index is TEMPLATED off a single PPP protein's stacked regenie sumstats. The probes
showed PPP coverage of both reference lists is protein-invariant (rabgap covers
99.24% of HapMap3 and 93.88% of 1kg-common; a second protein recovers only 121 of
~538k missing common variants, so the gap is panel-wide, not per-protein QC).
Templating loses nothing a union of proteins would recover, and gives us, from one
self-consistent source:
  - hg38 position   (GENPOS)          -> primary POS
  - hg19 position   (parsed from ID)  -> secondary POS_HG19 (complete, no liftover)
  - effect / non-effect alleles       (ALLELE1 / ALLELE0)  -> PPP-native orientation
  - in-sample effect-allele frequency (A1FREQ)             -> EAF

The membership task supplies the variant UNIVERSE and rsID. It must expose a
normalized frame with columns (CHR, POS[hg38], EA, NEA, rsID); matching is
allele-aware on the unordered {EA, NEA} set so a swapped reference orientation
still matches, after which PPP's orientation is adopted. The mode (HapMap3 vs
common-1kg) is not stored on this task: it is a property of the membership task.
"""

from pathlib import Path, PurePath

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_rsid import (
    is_palindromic_expr,
)
from mecfs_bio.build_system.task.ppp_database.allele_key import unordered_allele_key
from mecfs_bio.build_system.task.ppp_database.byte_stream_split_parquet import (
    write_byte_stream_split_parquet,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)
from mecfs_bio.constants.ppp_index_constants import (
    PPP_INDEX_ALLELE_KEY_COL,
    PPP_INDEX_IS_STRAND_AMBIGUOUS_COL,
    PPP_INDEX_POS_HG19_COL,
)
from mecfs_bio.constants.regenie_constants import (
    REGENIE_A1FREQ_COL,
    REGENIE_ALLELE0_COL,
    REGENIE_ALLELE1_COL,
    REGENIE_CHROM_COL,
    REGENIE_GENPOS_COL,
    REGENIE_ID_COL,
)

INDEX_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    PPP_INDEX_POS_HG19_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_RSID_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    PPP_INDEX_IS_STRAND_AMBIGUOUS_COL,
]

# The only float column, and thus the one that benefits from byte-stream-split.
_FLOAT_COLUMNS = [GWASLAB_EFFECT_ALLELE_FREQ_COL]


@frozen
class ConstructPppVariantIndexTask(GeneratingTask):
    """
    Build the shared PPP variant index by intersecting a template protein's
    variants with a reference membership list: the reference list's variant
    universe, with the protein's allele orientation and coordinates.

    template_protein_task: a stacked PPP protein (regenie columns CHROM, GENPOS,
        ID, ALLELE0, ALLELE1, A1FREQ). Supplies coordinates (both builds), alleles
        and in-sample EAF.
    membership_task: a normalized reference list exposing (CHR, POS[hg38], EA, NEA,
        rsID). Supplies the variant universe and rsID, and implicitly the mode.
    """

    meta: Meta
    template_protein_task: Task
    membership_task: Task

    @property
    def template_protein_meta(self) -> Meta:
        return self.template_protein_task.meta

    @property
    def membership_meta(self) -> Meta:
        return self.membership_task.meta

    @property
    def deps(self) -> list[Task]:
        return [self.template_protein_task, self.membership_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        # --- Template protein: coordinates (both builds), alleles, in-sample EAF. ---
        protein_asset = fetch(self.template_protein_task.asset_id)
        protein = (
            scan_dataframe_asset(
                protein_asset,
                meta=self.template_protein_meta,
                parquet_backend="polars",
            )
            .to_native()
            .select(
                pl.col(REGENIE_CHROM_COL).cast(pl.Int32).alias(GWASLAB_CHROM_COL),
                pl.col(REGENIE_GENPOS_COL).cast(pl.Int32).alias(GWASLAB_POS_COL),
                # hg19 position is field 1 of the ID; present for every PPP row.
                pl.col(REGENIE_ID_COL)
                .str.split(":")
                .list.get(1)
                .cast(pl.Int32)
                .alias(PPP_INDEX_POS_HG19_COL),
                pl.col(REGENIE_ALLELE1_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
                pl.col(REGENIE_ALLELE0_COL).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
                pl.col(REGENIE_A1FREQ_COL)
                .cast(pl.Float32)
                .alias(GWASLAB_EFFECT_ALLELE_FREQ_COL),
            )
            .with_columns(
                unordered_allele_key(
                    GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
                ).alias(PPP_INDEX_ALLELE_KEY_COL)
            )
            .unique(
                subset=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, PPP_INDEX_ALLELE_KEY_COL]
            )
        )

        # --- Membership reference list: variant universe + rsID. ---
        membership_asset = fetch(self.membership_task.asset_id)
        membership = (
            scan_dataframe_asset(
                membership_asset,
                meta=self.membership_meta,
                parquet_backend="polars",
            )
            .to_native()
            .select(
                pl.col(GWASLAB_CHROM_COL).cast(pl.Int32),
                pl.col(GWASLAB_POS_COL).cast(pl.Int32),
                pl.col(GWASLAB_RSID_COL),
                unordered_allele_key(
                    GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
                ).alias(PPP_INDEX_ALLELE_KEY_COL),
            )
            .unique(
                subset=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, PPP_INDEX_ALLELE_KEY_COL]
            )
        )

        index = (
            protein.join(
                membership,
                on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, PPP_INDEX_ALLELE_KEY_COL],
                how="inner",
            )
            .with_columns(
                is_palindromic_expr(
                    GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
                ).alias(PPP_INDEX_IS_STRAND_AMBIGUOUS_COL)
            )
            # Row order IS the alignment slot; sort fully deterministically. NEA is
            # included so multiallelic sites (shared position, differing alleles)
            # order stably even without assuming NEA == reference allele.
            .sort(
                [
                    GWASLAB_CHROM_COL,
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL,
                ]
            )
            .select(INDEX_COLUMNS)
            .collect()
        )

        out_path = scratch_dir / "ppp_variant_index.parquet"
        write_byte_stream_split_parquet(index, out_path, float_columns=_FLOAT_COLUMNS)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        template_protein_task: Task,
        membership_task: Task,
        asset_id: str,
    ) -> "ConstructPppVariantIndexTask":
        meta = ReferenceFileMeta(
            group="ukbb_ppp_variant_index",
            # Derived from the membership list so the two modes get distinct paths
            # without a redundant, typo-prone mode argument.
            sub_group=str(membership_task.asset_id),
            sub_folder=PurePath("processed"),
            id=AssetId(asset_id),
            filename="ppp_variant_index",
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
        return cls(
            meta=meta,
            template_protein_task=template_protein_task,
            membership_task=membership_task,
        )
