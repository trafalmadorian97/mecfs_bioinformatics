"""
Build the slim per-protein aligned beta/se file for one UKB-PPP protein.

This is the storage-saving production task: it downloads the protein's tar from
Synapse to a scratch directory, aligns its variants onto a shared variant index,
writes only beta and se (float32) in index row order as a Zstd + byte-stream-split
parquet, and lets the bulky download be discarded with the scratch dir. The full
summary statistics are never materialized as an asset.

The format choice (parquet + Zstd + byte-stream-split, default level) is baked in:
a benchmark (experiments/claude/ppp_database/aligner_format_bench) showed
byte-stream-split cuts the file ~22% (near-free), higher Zstd levels add almost
nothing, and duckdb/ALP is larger.

Memory is bounded for parallel execution: the per-chromosome regenie files are
processed one chromosome at a time (only the needed columns, downcast to compact
types), and results are streamed to the output a row group at a time, so peak
memory is roughly one chromosome plus the (allele-only) index -- not the whole
15.9M-row protein.
"""

import tarfile
from pathlib import Path, PurePath

import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq
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
from mecfs_bio.build_system.task.ppp_database.allele_key import unordered_allele_key
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_CHROM_NAME_FOR_CODE,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.regenie_constants import (
    REGENIE_ALLELE0_COL,
    REGENIE_ALLELE1_COL,
    REGENIE_BETA_COL,
    REGENIE_CHROM_COL,
    REGENIE_GENPOS_COL,
    REGENIE_SE_COL,
)

logger = structlog.get_logger()

ALIGNED_COLUMNS = [GWASLAB_BETA_COL, GWASLAB_SE_COL]

# Only these regenie columns are needed to align (the ID/hg19 position, A1FREQ,
# INFO, N, ... are dropped at read time to keep memory low).
_PROTEIN_READ_COLUMNS = [
    REGENIE_CHROM_COL,
    REGENIE_GENPOS_COL,
    REGENIE_ALLELE0_COL,
    REGENIE_ALLELE1_COL,
    REGENIE_BETA_COL,
    REGENIE_SE_COL,
]

_ALIGNED_SCHEMA = pa.schema(
    [(GWASLAB_BETA_COL, pa.float32()), (GWASLAB_SE_COL, pa.float32())]
)

_ROW_COL = "__row__"
_ALLELE_KEY = "__allele_key__"
_PROTEIN_EA = "__protein_ea__"
_PROTEIN_NEA = "__protein_nea__"
_PROTEIN_BETA = "__protein_beta__"
_PROTEIN_SE = "__protein_se__"

# Index columns needed for alignment (the rest -- rsID, EAF, ... -- are not read).
_INDEX_ALIGN_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
]


def align_protein_to_index(
    index_df: pl.DataFrame, protein_df: pl.DataFrame
) -> pl.DataFrame:
    """Return a (beta, se) float32 frame in index row order.

    index_df needs CHR, POS, EA, NEA. protein_df needs the regenie columns CHROM,
    GENPOS, ALLELE0, ALLELE1, BETA, SE. Beta is oriented to the index effect allele;
    variants absent from the protein get NaN beta and se.
    """
    index_keyed = index_df.with_row_index(_ROW_COL).with_columns(
        unordered_allele_key(
            GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
        ).alias(_ALLELE_KEY)
    )
    protein_keyed = (
        protein_df.select(
            pl.col(REGENIE_CHROM_COL).cast(pl.Int32).alias(GWASLAB_CHROM_COL),
            pl.col(REGENIE_GENPOS_COL).cast(pl.Int32).alias(GWASLAB_POS_COL),
            pl.col(REGENIE_ALLELE1_COL).alias(_PROTEIN_EA),
            pl.col(REGENIE_ALLELE0_COL).alias(_PROTEIN_NEA),
            pl.col(REGENIE_BETA_COL).cast(pl.Float64).alias(_PROTEIN_BETA),
            pl.col(REGENIE_SE_COL).cast(pl.Float64).alias(_PROTEIN_SE),
        )
        .with_columns(
            unordered_allele_key(_PROTEIN_EA, _PROTEIN_NEA).alias(_ALLELE_KEY)
        )
        .unique(subset=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, _ALLELE_KEY])
    )

    joined = index_keyed.join(
        protein_keyed,
        on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, _ALLELE_KEY],
        how="left",
    ).sort(_ROW_COL)

    beta = (
        pl.when(pl.col(_PROTEIN_EA).is_null())
        .then(pl.lit(float("nan")))
        .when(pl.col(GWASLAB_EFFECT_ALLELE_COL) == pl.col(_PROTEIN_EA))
        .then(pl.col(_PROTEIN_BETA))
        .otherwise(-pl.col(_PROTEIN_BETA))
        .cast(pl.Float32)
        .alias(GWASLAB_BETA_COL)
    )
    se = (
        pl.when(pl.col(_PROTEIN_SE).is_null())
        .then(pl.lit(float("nan")))
        .otherwise(pl.col(_PROTEIN_SE))
        .cast(pl.Float32)
        .alias(GWASLAB_SE_COL)
    )
    return joined.select(beta, se)


def _read_protein_chromosome(gz_path: Path) -> pl.DataFrame:
    """Read only the alignment columns from one per-chromosome regenie file."""
    return pl.read_csv(gz_path, separator=" ", columns=_PROTEIN_READ_COLUMNS)


def write_slim_aligned_parquet(
    extracted_dir: Path, index_df: pl.DataFrame, out_path: Path
) -> None:
    """Align an extracted protein (per-chromosome regenie files under extracted_dir)
    onto index_df, writing beta/se in index row order as Zstd + byte-stream-split
    parquet, one chromosome per row group to bound memory.

    Relies on the index being sorted by chromosome (as ConstructPppVariantIndexTask
    produces it) so that processing chromosomes in first-appearance order reproduces
    the exact index row order.

    Every chromosome in the index must have exactly one regenie file; a missing file
    is an error, not an empty block. Individual index variants the protein does not
    carry are still a normal, expected NaN.
    """
    assert index_df[GWASLAB_CHROM_COL].is_sorted(), "index must be chromosome-sorted"
    chromosomes = index_df[GWASLAB_CHROM_COL].unique(maintain_order=True).to_list()

    with pq.ParquetWriter(
        str(out_path),
        _ALIGNED_SCHEMA,
        compression="zstd",
        use_dictionary=False,
        use_byte_stream_split=True,
    ) as writer:
        for chromosome in chromosomes:
            index_chunk = index_df.filter(
                pl.col(GWASLAB_CHROM_COL) == chromosome
            ).select(_INDEX_ALIGN_COLUMNS)
            # PPP names its files with letters (discovery_chrX_...) while coding the
            # chromosome numerically inside them (CHROM=23), so the filename needs the
            # letter back even though the join below works off the code.
            chromosome_name = GWASLAB_CHROM_NAME_FOR_CODE.get(chromosome, None) or str(
                chromosome
            )
            # Fail fast: a chromosome with no file would NaN out that whole block,
            # wasting a ~545MB download to produce a silently useless file.
            matches = list(extracted_dir.rglob(f"*chr{chromosome_name}_*.gz"))
            assert len(matches) == 1, (
                f"expected exactly one regenie file for chr{chromosome_name} under "
                f"{extracted_dir}, found {len(matches)}: {matches}"
            )
            aligned = align_protein_to_index(
                index_chunk, _read_protein_chromosome(matches[0])
            )
            writer.write_table(
                aligned.select(ALIGNED_COLUMNS).to_arrow().cast(_ALIGNED_SCHEMA)
            )


def extract_protein_tar(tar_path: Path, dest_dir: Path) -> Path:
    """Extract a protein tar (per-chromosome gzipped regenie files) into dest_dir,
    returning the extraction directory."""
    extract_dir = dest_dir / "extracted"
    extract_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r") as tar_object:
        tar_object.extractall(extract_dir)
    return extract_dir


@frozen
class BuildSlimProteinParquetTask(GeneratingTask):
    """
    Download one UKB-PPP protein from Synapse and store only its aligned beta/se in
    the variant index's row order, discarding the bulky download.

    index_task: a ConstructPppVariantIndexTask output (CHR, POS, EA, NEA, ...).
    synid / expected_tar_filename: Synapse entity and tar name for this protein.
    """

    meta: Meta
    index_task: Task
    synid: str
    expected_tar_filename: str

    @property
    def index_meta(self) -> Meta:
        return self.index_task.meta

    @property
    def deps(self) -> list[Task]:
        return [self.index_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        index_df = (
            scan_dataframe_asset(
                fetch(self.index_task.asset_id),
                meta=self.index_meta,
                parquet_backend="polars",
            )
            .to_native()
            .select(_INDEX_ALIGN_COLUMNS)
            .collect()
        )
        download_dir = scratch_dir / "download"
        download_dir.mkdir(parents=True, exist_ok=True)
        tar_path = wf.download_from_synapse(
            self.synid, download_dir, self.expected_tar_filename
        )
        extracted_dir = extract_protein_tar(tar_path, download_dir)
        out_path = scratch_dir / f"{self.meta.asset_id}.parquet.zstd"
        write_slim_aligned_parquet(extracted_dir, index_df, out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        index_task: Task,
        synid: str,
        expected_tar_filename: str,
        gene: str,
        asset_id: str,
        index_name: str,
    ) -> "BuildSlimProteinParquetTask":
        # No protein asset dependency (by design, to avoid materializing the full
        # sumstats), so trait/project come from the manifest-supplied gene, not a dep.
        meta = GWASSummaryDataFileMeta(
            id=AssetId(asset_id),
            trait="ukbb_ppp",
            project=gene,
            sub_dir="aligned",
            project_path=PurePath(f"{index_name}_index/{asset_id}.parquet.zstd"),
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            extension=".parquet",
        )
        return cls(
            meta=meta,
            index_task=index_task,
            synid=synid,
            expected_tar_filename=expected_tar_filename,
        )
