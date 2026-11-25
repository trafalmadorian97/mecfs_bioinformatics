"""
The shared variant index of a deCODE SomaScan group: the rows of a template protein file whose
variant is in a membership list (for example HapMap3), in file order.

Protein files are matched to the index by NAME (deCODE's chrom:pos:effectAllele:otherAllele),
which is unique in the index. Protein files mostly list the template's rows in the template's
order, but some lack a few rows, so a row's number in the template does not locate it elsewhere.

Matching to the membership list is on (CHR, POS, unordered {EA, NEA}), so a swapped orientation
still matches, and the deCODE orientation is adopted. Rows with the missing-allele marker "!"
are dropped. EAF is the true effect allele frequency from variant_info_iceland (see
effect_allele_frequency); the task fails if any index row's EAF is undetermined.
"""

from pathlib import Path, PurePath

import numpy as np
import polars as pl
import rapidgzip
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
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
from mecfs_bio.build_system.task.dataframe_output import write_parquet_table
from mecfs_bio.build_system.task.decode_database.decode_somascan_protein import (
    DecodeSomascanProtein,
)
from mecfs_bio.build_system.task.decode_database.decode_sumstats_format import (
    chrom_code_expr,
    parse_decode_lines,
)
from mecfs_bio.build_system.task.decode_database.effect_allele_frequency import (
    resolve_effect_allele_frequency,
)
from mecfs_bio.build_system.task.decode_database.line_blocks import line_blocks
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_rsid import (
    is_palindromic_expr,
)
from mecfs_bio.build_system.task.ppp_database.allele_key import (
    assert_all_snv,
    unordered_allele_key,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.decode_somascan_constants import (
    DECODE_CHROM_COL,
    DECODE_EFFECT_ALLELE_COL,
    DECODE_IMP_MAF_COL,
    DECODE_INDEX_NAME_COL,
    DECODE_MISSING_ALLELE,
    DECODE_NAME_COL,
    DECODE_OTHER_ALLELE_COL,
    DECODE_POS_COL,
    VI_CHR_COL,
    VI_POS_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)
from mecfs_bio.constants.ppp_database_constants import (
    PPP_INDEX_ALLELE_KEY_COL,
    PPP_INDEX_IS_STRAND_AMBIGUOUS_COL,
)

INDEX_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_RSID_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    PPP_INDEX_IS_STRAND_AMBIGUOUS_COL,
    DECODE_INDEX_NAME_COL,
]
# Row number in the template, to keep the index in file order.
_ROW_COL = "__template_row__"
_JOIN_KEYS = [GWASLAB_CHROM_COL, GWASLAB_POS_COL, PPP_INDEX_ALLELE_KEY_COL]
_PRESIGN_EXPIRES_S = 6 * 3600
_DECOMPRESS_PARALLELIZATION = 8


@frozen(slots=True)
class ConstructDecodeVariantIndexTask(GeneratingTask):
    """
    template: the protein file whose rows define the index.
    membership_task: a normalized membership list (CHR, POS, EA, NEA, rsID).
    variant_info_task: the deduplicated variant_info parquet (BuildDecodeVariantInfoTask).
    """

    meta: ReferenceFileMeta
    template: DecodeSomascanProtein
    membership_task: Task
    variant_info_task: Task

    @property
    def deps(self) -> list[Task]:
        return [self.membership_task, self.variant_info_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        download = scratch_dir / "template.txt.gz"
        download_protein_file(self.template, download, wf)
        membership = _read_membership(self.membership_task, fetch)
        matched = _match_template_rows(download, membership)
        download.unlink()
        assert matched[DECODE_INDEX_NAME_COL].is_unique().all(), (
            f"{self.template.key}: index NAME values are not unique"
        )
        eaf = _effect_allele_frequency(matched, self.variant_info_task, fetch)
        undetermined = eaf.null_count()
        assert undetermined == 0, (
            f"{undetermined} index rows of {self.template.key} have no determined EAF"
        )
        index = _index_frame(matched, eaf)
        assert_all_snv(index, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL)
        out_path = scratch_dir / "decode_variant_index.parquet"
        _write_index(index, out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        template: DecodeSomascanProtein,
        membership_task: Task,
        variant_info_task: Task,
    ) -> "ConstructDecodeVariantIndexTask":
        return cls(
            meta=ReferenceFileMeta(
                group="decode_somascan_variant_index",
                sub_group=str(membership_task.asset_id),
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                filename=asset_id,
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
            template=template,
            membership_task=membership_task,
            variant_info_task=variant_info_task,
        )


def download_protein_file(
    protein: DecodeSomascanProtein, local_path: Path, wf: WF, n_connections: int = 128
) -> None:
    """Download a protein file after checking that the bucket still has the manifest's bytes."""
    info = wf.s3_access.head(protein.ref)
    assert info.size == protein.size and info.md5 == protein.md5, (
        f"{protein.key} changed since the manifest was made: size {info.size} md5 "
        f"{info.md5}, manifest size {protein.size} md5 {protein.md5}"
    )
    wf.download_url_parallel(
        url=wf.s3_access.presigned_url(protein.ref, _PRESIGN_EXPIRES_S),
        size=protein.size,
        local_path=local_path,
        md5_hash=protein.md5,
        n_connections=n_connections,
    )


def _read_membership(membership_task: Task, fetch: Fetch) -> pl.DataFrame:
    return (
        scan_dataframe_asset(
            fetch(membership_task.asset_id),
            meta=membership_task.meta,
            parquet_backend="polars",
        )
        .to_native()
        .select(
            pl.col(GWASLAB_CHROM_COL).cast(pl.Int8),
            pl.col(GWASLAB_POS_COL).cast(pl.Int64),
            pl.col(GWASLAB_RSID_COL),
            unordered_allele_key(
                GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
            ).alias(PPP_INDEX_ALLELE_KEY_COL),
        )
        .unique(subset=_JOIN_KEYS)
        .collect()
    )


def _match_template_rows(gz_path: Path, membership: pl.DataFrame) -> pl.DataFrame:
    """Parse the whole template block by block, keeping the rows in membership."""
    parts: list[pl.DataFrame] = []
    with rapidgzip.RapidgzipFile(
        str(gz_path), parallelization=_DECOMPRESS_PARALLELIZATION
    ) as f:
        for block in line_blocks(f).blocks:
            parsed = parse_decode_lines(block.text).with_columns(
                pl.Series(
                    _ROW_COL,
                    np.arange(block.first_row, block.first_row + block.n_rows),
                    dtype=pl.UInt64,
                )
            )
            assert parsed.height == block.n_rows
            parts.append(_keep_members(parsed, membership))
    return pl.concat(parts)


def _keep_members(parsed: pl.DataFrame, membership: pl.DataFrame) -> pl.DataFrame:
    return (
        parsed.filter(pl.col(DECODE_OTHER_ALLELE_COL) != DECODE_MISSING_ALLELE)
        .select(
            chrom_code_expr(DECODE_CHROM_COL).alias(GWASLAB_CHROM_COL),
            pl.col(DECODE_POS_COL).alias(GWASLAB_POS_COL),
            pl.col(DECODE_EFFECT_ALLELE_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
            pl.col(DECODE_OTHER_ALLELE_COL).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
            pl.col(DECODE_IMP_MAF_COL),
            pl.col(_ROW_COL),
            pl.col(DECODE_NAME_COL).alias(DECODE_INDEX_NAME_COL),
        )
        .with_columns(
            unordered_allele_key(
                GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
            ).alias(PPP_INDEX_ALLELE_KEY_COL)
        )
        .join(membership, on=_JOIN_KEYS, how="inner")
    )


def _effect_allele_frequency(
    rows: pl.DataFrame, variant_info_task: Task, fetch: Fetch
) -> pl.Series:
    positions = rows.select(
        pl.col(GWASLAB_CHROM_COL).alias(VI_CHR_COL),
        pl.col(GWASLAB_POS_COL).alias(VI_POS_COL),
    ).unique()
    variant_info = (
        scan_dataframe_asset(
            fetch(variant_info_task.asset_id),
            meta=variant_info_task.meta,
            parquet_backend="polars",
        )
        .to_native()
        .join(positions.lazy(), on=[VI_CHR_COL, VI_POS_COL], how="semi")
        .collect()
    )
    return resolve_effect_allele_frequency(rows, variant_info)


def _index_frame(rows: pl.DataFrame, eaf: pl.Series) -> pl.DataFrame:
    return (
        rows.with_columns(
            eaf.cast(pl.Float32).alias(GWASLAB_EFFECT_ALLELE_FREQ_COL),
            is_palindromic_expr(
                GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
            ).alias(PPP_INDEX_IS_STRAND_AMBIGUOUS_COL),
        )
        .sort(_ROW_COL)
        .select(INDEX_COLUMNS)
    )


def _write_index(index: pl.DataFrame, out_path: Path) -> None:
    write_parquet_table(
        table=index.to_arrow(),
        out_path=out_path,
        compression="zstd",
        compression_level=None,
        byte_stream_split_columns=[GWASLAB_EFFECT_ALLELE_FREQ_COL],
    )
