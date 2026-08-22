"""
Build the slim per-aptamer aligned beta/se/N file for one Western et al. 2024 CSF
aptamer.

This is the storage-saving production task, analogous to
BuildSlimProteinParquetTask for UKB-PPP. It downloads the aptamer's GWAS-SSF
summary statistics (~203 MB .tsv.gz) from the GWAS Catalog into a scratch
directory, aligns its variants onto the shared CSF variant index, writes only
beta, se and N (float32) in index row order as a Zstd + byte-stream-split parquet,
and lets the bulky download be discarded with the scratch dir. The full summary
statistics are never materialized as an asset.

Simpler than the PPP analogue in three ways: the source is a single file (no tar,
no per-chromosome loop), it arrives over plain HTTPS (no Synapse), and N varies per
variant (it is stored per row, not recovered separately as a constant).

Unlike PPP, the whole aligned frame is written in one row group. The single-file
read already holds the source in memory, so a per-chromosome streaming write would
not lower peak memory here; it would only add complexity. Peak memory is dominated
by the ~0.6-1 GB single-file read, which caps useful parallelism.
"""

from pathlib import Path, PurePath

import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.task.csf_database.gwas_catalog_url import (
    gwas_catalog_sumstats_url,
)
from mecfs_bio.build_system.task.dataframe_output import write_parquet_table
from mecfs_bio.build_system.task.ppp_database.allele_key import unordered_allele_key
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.csf_database_constants import (
    CSF_INDEX_ALLELE_KEY_COL,
    Analyte,
    GcstAccession,
    SeqId,
    UniProtId,
)
from mecfs_bio.constants.gwas_ssf_constants import (
    GWAS_SSF_BETA_COL,
    GWAS_SSF_CHROM_COL,
    GWAS_SSF_EFFECT_ALLELE_COL,
    GWAS_SSF_N_COL,
    GWAS_SSF_OTHER_ALLELE_COL,
    GWAS_SSF_POS_COL,
    GWAS_SSF_SE_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

logger = structlog.get_logger()

# The GWAS Catalog FTP throttles per connection (~0.5 MB/s), so a single-connection
# pull of one ~193 MB aptamer file takes ~5 min; with 16 connections the same file
# lands in ~12 s. Measured end to end, that is the difference between a ~3-week and a
# ~2-day full build (7,008 files, sequential). See experiments/claude/logs for the
# benchmark.
#
# 16 is also aria2's hard maximum for --max-connection-per-server; a larger value makes
# aria2 error out, not go faster.
_DOWNLOAD_CONNECTIONS = 16

# Alignment output: beta, se, N, all float32. N is conceptually integer but stored as
# float so downstream numeric users need no conversion (mirrors PPP).
_ALIGNED_COLUMNS = [GWASLAB_BETA_COL, GWASLAB_SE_COL, GWASLAB_SAMPLE_SIZE_COLUMN]

# Index columns needed for alignment (rsID, EAF, ... are not read).
_INDEX_ALIGN_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
]

# The GWAS-SSF columns the alignment reads; the rest (EAF, -log10 p, variant_id, rsID)
# are dropped at read time.
_APTAMER_READ_COLUMNS = [
    GWAS_SSF_CHROM_COL,
    GWAS_SSF_POS_COL,
    GWAS_SSF_EFFECT_ALLELE_COL,
    GWAS_SSF_OTHER_ALLELE_COL,
    GWAS_SSF_BETA_COL,
    GWAS_SSF_SE_COL,
    GWAS_SSF_N_COL,
]

_ROW_COL = "__row__"
_APTAMER_EA = "__aptamer_ea__"
_APTAMER_NEA = "__aptamer_nea__"
_APTAMER_BETA = "__aptamer_beta__"
_APTAMER_SE = "__aptamer_se__"
_APTAMER_N = "__aptamer_n__"


def align_aptamer_to_index(
    index_df: pl.DataFrame,
    aptamer_df: pl.DataFrame,
) -> pl.DataFrame:
    """Return a (beta, se, N) float32 frame in index row order.

    index_df needs CHR, POS, EA, NEA. aptamer_df needs the GWAS-SSF columns
    chromosome, base_pair_location, effect_allele, other_allele, beta,
    standard_error, n. Beta is oriented to the index effect allele; se and N are
    orientation-independent. Variants absent from the aptamer get NaN in every
    output column.

    Join misses are represented as float("nan"), not polars null, so the columns
    convert zero-copy to numpy downstream (see BuildSlimProteinParquetTask for the
    same rationale). N is stored as float32 even though it is conceptually integer,
    to save downstream users a conversion.
    """
    index_keyed = index_df.with_row_index(_ROW_COL).with_columns(
        unordered_allele_key(
            GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
        ).alias(CSF_INDEX_ALLELE_KEY_COL)
    )
    aptamer_keyed = (
        aptamer_df.select(
            pl.col(GWAS_SSF_CHROM_COL).cast(pl.Int32).alias(GWASLAB_CHROM_COL),
            pl.col(GWAS_SSF_POS_COL).cast(pl.Int32).alias(GWASLAB_POS_COL),
            pl.col(GWAS_SSF_EFFECT_ALLELE_COL).alias(_APTAMER_EA),
            pl.col(GWAS_SSF_OTHER_ALLELE_COL).alias(_APTAMER_NEA),
            pl.col(GWAS_SSF_BETA_COL).cast(pl.Float64).alias(_APTAMER_BETA),
            pl.col(GWAS_SSF_SE_COL).cast(pl.Float64).alias(_APTAMER_SE),
            pl.col(GWAS_SSF_N_COL).cast(pl.Float64).alias(_APTAMER_N),
        )
        .with_columns(
            unordered_allele_key(_APTAMER_EA, _APTAMER_NEA).alias(
                CSF_INDEX_ALLELE_KEY_COL
            )
        )
        .unique(subset=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, CSF_INDEX_ALLELE_KEY_COL])
    )

    joined = index_keyed.join(
        aptamer_keyed,
        on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, CSF_INDEX_ALLELE_KEY_COL],
        how="left",
    ).sort(_ROW_COL)

    beta = (
        pl.when(pl.col(_APTAMER_EA).is_null())
        .then(pl.lit(float("nan")))
        .when(pl.col(GWASLAB_EFFECT_ALLELE_COL) == pl.col(_APTAMER_EA))
        .then(pl.col(_APTAMER_BETA))
        .otherwise(-pl.col(_APTAMER_BETA))
        .cast(pl.Float32)
        .alias(GWASLAB_BETA_COL)
    )
    se = (
        pl.when(pl.col(_APTAMER_SE).is_null())
        .then(pl.lit(float("nan")))
        .otherwise(pl.col(_APTAMER_SE))
        .cast(pl.Float32)
        .alias(GWASLAB_SE_COL)
    )
    # N keys off the same join-miss indicator as beta: absent variant -> NaN.
    n = (
        pl.when(pl.col(_APTAMER_EA).is_null())
        .then(pl.lit(float("nan")))
        .otherwise(pl.col(_APTAMER_N))
        .cast(pl.Float32)
        .alias(GWASLAB_SAMPLE_SIZE_COLUMN)
    )
    return joined.select([beta, se, n])


def read_aptamer_sumstats(gz_path: Path) -> pl.DataFrame:
    """Read only the alignment columns from an aptamer's GWAS-SSF .tsv.gz.

    Polars infers gzip from the file's magic bytes, so the download need not carry a
    .gz extension.
    """
    return pl.read_csv(gz_path, separator="\t", columns=_APTAMER_READ_COLUMNS)


def write_slim_aptamer_parquet(
    aptamer_df: pl.DataFrame,
    index_df: pl.DataFrame,
    out_path: Path,
) -> None:
    """Align aptamer_df onto index_df and write beta/se/N (float32) in index row order
    as a Zstd level-15 + byte-stream-split parquet.

    Level 15 is the measured sweet spot (see the CSF database plan): it captures
    essentially all the compression gain of the maximum level at a quarter of the
    write cost, and compression level does not affect downstream read speed.
    """
    aligned = align_aptamer_to_index(index_df, aptamer_df).select(_ALIGNED_COLUMNS)
    write_parquet_table(
        table=aligned.to_arrow(),
        out_path=out_path,
        compression="zstd",
        compression_level=15,
        byte_stream_split_columns=_ALIGNED_COLUMNS,
    )


@frozen
class CsfAptamerFile:
    """Structured identity of one CSF aptamer's summary statistics (a manifest row).

    analyte is the SomaScan analyte id (primary key, e.g. X13681.173); seq_id is the
    same aptamer in dash form (13681-173); accession is the GWAS Catalog study whose
    GWAS-SSF file holds the sumstats. uniprot is the stable database key for the protein
    target; entrez_gene_symbol names the target for asset layout, but neither identifies
    the aptamer (both can be shared across aptamers).
    """

    analyte: Analyte
    seq_id: SeqId
    accession: GcstAccession
    uniprot: UniProtId
    entrez_gene_symbol: str

    @property
    def sumstats_url(self) -> str:
        return gwas_catalog_sumstats_url(self.accession)


@frozen
class BuildSlimCsfAptamerParquetTask(GeneratingTask):
    """
    Download one CSF aptamer's GWAS-SSF summary statistics from the GWAS Catalog and
    store only its aligned beta/se/N in the variant index's row order, discarding the
    bulky download.

    index_task: a ConstructCsfVariantIndexTask output (CHR, POS, EA, NEA, ...).
    aptamer: the structured identity of the aptamer (analyte/accession/gene).
    md5_hash: the published md5 of the .tsv.gz, verified on download, or None to skip
        verification.
    """

    meta: Meta
    index_task: Task
    aptamer: CsfAptamerFile
    md5_hash: str | None

    @property
    def index_meta(self) -> Meta:
        return self.index_task.meta

    @property
    def deps(self) -> list[Task]:
        return [self.index_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        # Keep the column projection on the narwhals frame (pushes into scan_parquet),
        # then collect().to_polars() to guarantee a polars frame regardless of the
        # backend a future DataProcessingPipe might introduce (to_native would not).
        index_df = (
            scan_dataframe_asset(
                fetch(self.index_task.asset_id),
                meta=self.index_meta,
                parquet_backend="polars",
            )
            .select(_INDEX_ALIGN_COLUMNS)
            .collect()
            .to_polars()
        )
        download_dir = scratch_dir / "download"
        download_dir.mkdir(parents=True, exist_ok=True)
        gz_path = download_dir / f"{self.aptamer.accession}.tsv.gz"
        wf.download_from_url(
            url=self.aptamer.sumstats_url,
            local_path=gz_path,
            md5_hash=self.md5_hash,
            request_connections=_DOWNLOAD_CONNECTIONS,
        )
        out_path = scratch_dir / f"{self.meta.asset_id}.parquet.zstd"
        write_slim_aptamer_parquet(read_aptamer_sumstats(gz_path), index_df, out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        index_task: Task,
        aptamer: CsfAptamerFile,
        asset_id: str,
        index_name: str,
        md5_hash: str | None = None,
    ) -> "BuildSlimCsfAptamerParquetTask":
        # No aptamer asset dependency (by design, to avoid materializing the full
        # sumstats), so trait/project come from the aptamer's gene, not a dep.
        meta = GWASSummaryDataFileMeta(
            id=AssetId(asset_id),
            trait="western_csf_pqtl",
            project=aptamer.entrez_gene_symbol,
            sub_dir="aligned",
            project_path=PurePath(f"{index_name}_index/{asset_id}.parquet"),
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            extension=".parquet",
        )
        return cls(
            meta=meta,
            index_task=index_task,
            aptamer=aptamer,
            md5_hash=md5_hash,
        )
