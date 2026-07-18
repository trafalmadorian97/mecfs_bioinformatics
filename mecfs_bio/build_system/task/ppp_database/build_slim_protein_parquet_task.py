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
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.regenie_constants import (
    REGENIE_ALLELE0_COL,
    REGENIE_ALLELE1_COL,
    REGENIE_BETA_COL,
    REGENIE_CHROM_COL,
    REGENIE_GENPOS_COL,
    REGENIE_N_COL,
    REGENIE_SE_COL,
)

logger = structlog.get_logger()

# The columns needed to align (CHROM/GENPOS/alleles), plus beta/se and optionally N.
# The ID/hg19 position, A1FREQ, INFO, ... are dropped at read time to keep memory low; N
# is dropped too unless include_sample_size stores it (it is constant per protein, so it
# run-length-compresses to almost nothing but lets a downstream reader skip the separate
# ranged-read sample-size recovery).
_PROTEIN_READ_COLUMNS = [
    REGENIE_CHROM_COL,
    REGENIE_GENPOS_COL,
    REGENIE_ALLELE0_COL,
    REGENIE_ALLELE1_COL,
    REGENIE_BETA_COL,
    REGENIE_SE_COL,
]

_ROW_COL = "__row__"
_ALLELE_KEY = "__allele_key__"
_PROTEIN_EA = "__protein_ea__"
_PROTEIN_NEA = "__protein_nea__"
_PROTEIN_BETA = "__protein_beta__"
_PROTEIN_SE = "__protein_se__"
_PROTEIN_N = "__protein_n__"


def aligned_columns(include_sample_size: bool) -> list[str]:
    """Output columns of the slim parquet: beta/se, and N when include_sample_size."""
    columns = [GWASLAB_BETA_COL, GWASLAB_SE_COL]
    if include_sample_size:
        columns.append(GWASLAB_SAMPLE_SIZE_COLUMN)
    return columns


def _aligned_schema(include_sample_size: bool) -> pa.Schema:
    fields = [(GWASLAB_BETA_COL, pa.float32()), (GWASLAB_SE_COL, pa.float32())]
    if include_sample_size:
        fields.append((GWASLAB_SAMPLE_SIZE_COLUMN, pa.float32()))
    return pa.schema(fields)


# Index columns needed for alignment (the rest -- rsID, EAF, ... -- are not read).
_INDEX_ALIGN_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
]


def align_protein_to_index(
    index_df: pl.DataFrame,
    protein_df: pl.DataFrame,
    *,
    include_sample_size: bool,
) -> pl.DataFrame:
    """Return a (beta, se[, N]) float32 frame in index row order.

    index_df needs CHR, POS, EA, NEA. protein_df needs the regenie columns CHROM,
    GENPOS, ALLELE0, ALLELE1, BETA, SE (and N when include_sample_size). Beta is
    oriented to the index effect allele; N is allele-independent. Variants absent from
    the protein get NaN in every output column.

    NOTE:
        - The beta and se columns in the output dataframe are stored as float32s.
        - To represent cases in which the protein is missing a variant present in the index, we use float(nan) rather than
          polars' null value.  This deliberate choice allows zero-copy conversion to numpy arrays.  This increases throughput in downstream computations
          (see https://docs.pola.rs/api/python/stable/reference/series/api/polars.Series.to_numpy.html)
        - The sample size column is stored as float, even though it is conceptually an integer. This is because downstream users
          will typically want a float, so by taking this approach we avoid conversions.
    """
    index_keyed = index_df.with_row_index(_ROW_COL).with_columns(
        unordered_allele_key(
            GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
        ).alias(_ALLELE_KEY)
    )
    protein_select = [
        pl.col(REGENIE_CHROM_COL).cast(pl.Int32).alias(GWASLAB_CHROM_COL),
        pl.col(REGENIE_GENPOS_COL).cast(pl.Int32).alias(GWASLAB_POS_COL),
        pl.col(REGENIE_ALLELE1_COL).alias(_PROTEIN_EA),
        pl.col(REGENIE_ALLELE0_COL).alias(_PROTEIN_NEA),
        pl.col(REGENIE_BETA_COL).cast(pl.Float64).alias(_PROTEIN_BETA),
        pl.col(REGENIE_SE_COL).cast(pl.Float64).alias(_PROTEIN_SE),
    ]
    if include_sample_size:
        protein_select.append(pl.col(REGENIE_N_COL).cast(pl.Float64).alias(_PROTEIN_N))
    protein_keyed = (
        protein_df.select(protein_select)
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
    outputs = [beta, se]
    if include_sample_size:
        # N keys off the same join-miss indicator as beta: absent variant -> NaN.
        n = (
            pl.when(pl.col(_PROTEIN_EA).is_null())
            .then(pl.lit(float("nan")))
            .otherwise(pl.col(_PROTEIN_N))
            .cast(pl.Float32)
            .alias(GWASLAB_SAMPLE_SIZE_COLUMN)
        )
        outputs.append(n)
    return joined.select(outputs)


def _read_protein_chromosome(
    gz_path: Path, *, include_sample_size: bool
) -> pl.DataFrame:
    """Read only the alignment columns (plus N when include_sample_size) from one
    per-chromosome regenie file."""
    columns = _PROTEIN_READ_COLUMNS + ([REGENIE_N_COL] if include_sample_size else [])
    return pl.read_csv(gz_path, separator=" ", columns=columns)


def write_slim_aligned_parquet(
    extracted_dir: Path,
    index_df: pl.DataFrame,
    out_path: Path,
    *,
    include_sample_size: bool,
) -> None:
    """Align an extracted protein (per-chromosome regenie files under extracted_dir)
    onto index_df, writing beta/se (and N when include_sample_size) in index row order as
    Zstd + byte-stream-split parquet, one chromosome per row group to bound memory.

    Relies on the index being sorted by chromosome (as ConstructPppVariantIndexTask
    produces it) so that processing chromosomes in first-appearance order reproduces
    the exact index row order.

    Every chromosome in the index must have exactly one regenie file; a missing file
    is an error, not an empty block. Individual index variants the protein does not
    carry are still a normal, expected NaN.
    """
    assert index_df[GWASLAB_CHROM_COL].is_sorted(), "index must be chromosome-sorted"
    chromosomes = index_df[GWASLAB_CHROM_COL].unique(maintain_order=True).to_list()
    schema = _aligned_schema(include_sample_size)
    output_columns = aligned_columns(include_sample_size)

    with pq.ParquetWriter(
        str(out_path),
        schema,
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
                index_chunk,
                _read_protein_chromosome(
                    matches[0], include_sample_size=include_sample_size
                ),
                include_sample_size=include_sample_size,
            )
            writer.write_table(aligned.select(output_columns).to_arrow().cast(schema))


def extract_protein_tar(tar_path: Path, dest_dir: Path) -> Path:
    """Extract a protein tar (per-chromosome gzipped regenie files) into dest_dir,
    returning the extraction directory."""
    extract_dir = dest_dir / "extracted"
    extract_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r") as tar_object:
        tar_object.extractall(extract_dir)
    return extract_dir


@frozen
class PppProteinFile:
    """Structured identity of one UKB-PPP per-protein summary-statistics file (a manifest row).

    The Synapse tar is named gene_uniprot_OID_version_panel.tar, so its filename is derived
    from these fields (tar_filename) rather than stored as, or parsed back out of, a string.
    synid is the Synapse entity id of the tar (cohort-specific).
    """

    gene: str
    uniprot: str
    oid: str
    version: str
    panel: str
    synid: str

    @property
    def tar_filename(self) -> str:
        return f"{self.gene}_{self.uniprot}_{self.oid}_{self.version}_{self.panel}.tar"


@frozen
class BuildSlimProteinParquetTask(GeneratingTask):
    """
    Download one UKB-PPP protein from Synapse and store only its aligned beta/se (and
    optionally N) in the variant index's row order, discarding the bulky download.

    index_task: a ConstructPppVariantIndexTask output (CHR, POS, EA, NEA, ...).
    protein: the structured identity of the protein file (gene/OID/Synapse id/...).
    include_sample_size: also store the per-variant sample size N (constant per protein);
        it run-length-compresses to almost nothing but lets downstream readers skip the
        separate ranged-read N recovery. Kept False for cohorts whose slim files were
        already built without it, so their on-disk assets stay bit-identical.
    """

    meta: Meta
    index_task: Task
    protein: PppProteinFile
    include_sample_size: bool

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
            self.protein.synid, download_dir, self.protein.tar_filename
        )
        extracted_dir = extract_protein_tar(tar_path, download_dir)
        out_path = scratch_dir / f"{self.meta.asset_id}.parquet.zstd"
        write_slim_aligned_parquet(
            extracted_dir,
            index_df,
            out_path,
            include_sample_size=self.include_sample_size,
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        index_task: Task,
        protein: PppProteinFile,
        asset_id: str,
        index_name: str,
        include_sample_size: bool = True,
    ) -> "BuildSlimProteinParquetTask":
        # No protein asset dependency (by design, to avoid materializing the full
        # sumstats), so trait/project come from the protein's gene, not a dep.
        meta = GWASSummaryDataFileMeta(
            id=AssetId(asset_id),
            trait="ukbb_ppp",
            project=protein.gene,
            sub_dir="aligned",
            project_path=PurePath(f"{index_name}_index/{asset_id}.parquet.zstd"),
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            extension=".parquet",
        )
        return cls(
            meta=meta,
            index_task=index_task,
            protein=protein,
            include_sample_size=include_sample_size,
        )
