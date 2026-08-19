"""
Construct the shared variant index used to store Western et al. 2024 CSF pQTL
summary statistics compactly.

The index is a parquet table whose ROW ORDER is the canonical alignment slot for
every per-aptamer beta/se/N file: an aptamer's slim file stores only beta, se and N,
one row per index row, in this order.

Construction mirrors ConstructPppVariantIndexTask: the index is TEMPLATED off a
single aptamer's GWAS-SSF summary statistics, intersected with a reference
membership list (HapMap3). The CSF PLINK2 run is joint across aptamers, so one
template aptamer's variant set is representative; templating gives us, from one
self-consistent source, hg38 coordinates, effect/other alleles (the CSF-native
orientation) and the in-sample effect-allele frequency.

Unlike PPP there is no hg19 position: GWAS-SSF carries only the study build (hg38),
so POS_HG19 is omitted. Add a liftover task later only if something needs it.

The membership task supplies the variant UNIVERSE and rsID, exposing a normalized
frame with columns (CHR, POS[hg38], EA, NEA, rsID); matching is allele-aware on the
unordered {EA, NEA} set so a swapped reference orientation still matches, after which
the CSF orientation is adopted. The mode (currently only HapMap3) is a property of
the membership task, not stored here.
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
from mecfs_bio.constants.csf_database_constants import (
    CSF_INDEX_ALLELE_KEY_COL,
    CSF_INDEX_IS_STRAND_AMBIGUOUS_COL,
)
from mecfs_bio.constants.gwas_ssf_constants import (
    GWAS_SSF_CHROM_COL,
    GWAS_SSF_EFFECT_ALLELE_COL,
    GWAS_SSF_EFFECT_ALLELE_FREQ_COL,
    GWAS_SSF_OTHER_ALLELE_COL,
    GWAS_SSF_POS_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)

INDEX_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_RSID_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    CSF_INDEX_IS_STRAND_AMBIGUOUS_COL,
]

# The only float column, and thus the one that benefits from byte-stream-split.
_FLOAT_COLUMNS = [GWASLAB_EFFECT_ALLELE_FREQ_COL]


@frozen
class ConstructCsfVariantIndexTask(GeneratingTask):
    """
    Build the shared CSF variant index by intersecting a template aptamer's variants
    with a reference membership list: the reference list's variant universe, with the
    aptamer's allele orientation, coordinates and in-sample EAF.

    template_aptamer_task: an aptamer's GWAS-SSF sumstats parquet (columns
        chromosome, base_pair_location, effect_allele, other_allele,
        effect_allele_frequency). Supplies hg38 coordinates, alleles and EAF.
    membership_task: a normalized reference list exposing (CHR, POS[hg38], EA, NEA,
        rsID). Supplies the variant universe and rsID, and implicitly the mode.
    """

    meta: Meta
    template_aptamer_task: Task
    membership_task: Task

    @property
    def template_aptamer_meta(self) -> Meta:
        return self.template_aptamer_task.meta

    @property
    def membership_meta(self) -> Meta:
        return self.membership_task.meta

    @property
    def deps(self) -> list[Task]:
        return [self.template_aptamer_task, self.membership_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        # --- Template aptamer: hg38 coordinates, alleles, in-sample EAF. ---
        template_asset = fetch(self.template_aptamer_task.asset_id)
        template = (
            scan_dataframe_asset(
                template_asset,
                meta=self.template_aptamer_meta,
                parquet_backend="polars",
            )
            .to_native()
            .select(
                pl.col(GWAS_SSF_CHROM_COL).cast(pl.Int32).alias(GWASLAB_CHROM_COL),
                pl.col(GWAS_SSF_POS_COL).cast(pl.Int32).alias(GWASLAB_POS_COL),
                pl.col(GWAS_SSF_EFFECT_ALLELE_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
                pl.col(GWAS_SSF_OTHER_ALLELE_COL).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
                pl.col(GWAS_SSF_EFFECT_ALLELE_FREQ_COL)
                .cast(pl.Float32)
                .alias(GWASLAB_EFFECT_ALLELE_FREQ_COL),
            )
            .with_columns(
                unordered_allele_key(
                    GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
                ).alias(CSF_INDEX_ALLELE_KEY_COL)
            )
            .unique(
                subset=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, CSF_INDEX_ALLELE_KEY_COL]
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
                ).alias(CSF_INDEX_ALLELE_KEY_COL),
            )
            .unique(
                subset=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, CSF_INDEX_ALLELE_KEY_COL]
            )
        )

        index = (
            template.join(
                membership,
                on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, CSF_INDEX_ALLELE_KEY_COL],
                how="inner",
            )
            .with_columns(
                is_palindromic_expr(
                    GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
                ).alias(CSF_INDEX_IS_STRAND_AMBIGUOUS_COL)
            )
            # Row order IS the alignment slot; sort fully deterministically on all four
            # of (CHR, POS, EA, NEA). (CHR, POS) is unique in today's data, but that is
            # a property of one file, not a guarantee, and a silent reordering here
            # corrupts every downstream aptamer file at once.
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

        out_path = scratch_dir / "csf_variant_index.parquet"
        write_byte_stream_split_parquet(index, out_path, float_columns=_FLOAT_COLUMNS)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        template_aptamer_task: Task,
        membership_task: Task,
        asset_id: str,
    ) -> "ConstructCsfVariantIndexTask":
        meta = ReferenceFileMeta(
            group="csf_pqtl_variant_index",
            # Derived from the membership list so distinct modes get distinct paths
            # without a redundant, typo-prone mode argument.
            sub_group=str(membership_task.asset_id),
            sub_folder=PurePath("processed"),
            id=AssetId(asset_id),
            filename="csf_variant_index",
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
        return cls(
            meta=meta,
            template_aptamer_task=template_aptamer_task,
            membership_task=membership_task,
        )
