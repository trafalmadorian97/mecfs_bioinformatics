"""
Synthetic deCODE SomaScan inputs and fakes for the deCODE database tests.

synthetic_rows gives a small protein file in the exact deCODE layout: chromosomes in text order
(chr1, chr10, chr2, chr22, chrX), positions ascending within each. Besides the deCODE columns
it carries fixture-only columns (prefix "_fx_") that the writers drop:
  _fx_true_eaf: the effect allele frequency the variant_info writer encodes; null marks the
      ambiguous site whose variant_info has both f and 1 - f for the effect allele;
  _fx_special: rows that tests rely on staying outside any synthetic index (the "!" row,
      the ambiguous site, and the rows placed around the significant-rows threshold).
"""

import gzip
import hashlib
import shutil
from pathlib import Path

import numpy as np
import polars as pl
from attrs import field, frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.decode_database.build_decode_somascan_protein_task import (
    BuildDecodeSomascanProteinTask,
)
from mecfs_bio.build_system.task.decode_database.build_decode_variant_info_task import (
    write_variant_info_parquet,
)
from mecfs_bio.build_system.task.decode_database.construct_decode_variant_index_task import (
    ConstructDecodeVariantIndexTask,
)
from mecfs_bio.build_system.task.decode_database.decode_somascan_protein import (
    DecodeSomascanProtein,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.protein_coloc.decode_window_reader import (
    RangeSource,
    RangeSourceFactory,
)
from mecfs_bio.build_system.wf.base_wf import WF, make_wf
from mecfs_bio.build_system.wf.parallel_range_downloader import ParallelRangeDownloader
from mecfs_bio.build_system.wf.s3_access import (
    S3Access,
    S3Endpoint,
    S3ObjectInfo,
    S3ObjectRef,
)
from mecfs_bio.constants.csf_database_constants import SeqId, UniProtId
from mecfs_bio.constants.decode_somascan_constants import (
    DECODE_CHROM_COL,
    DECODE_EFFECT_ALLELE_COL,
    DECODE_HEADER,
    DECODE_IMP_MAF_COL,
    DECODE_MISSING_ALLELE,
    DECODE_MLOG10P_COL,
    DECODE_NAME_COL,
    DECODE_OTHER_ALLELE_COL,
    DECODE_POS_COL,
    DECODE_PVAL_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)
from mecfs_bio.util.download.verify import verify_hash

FX_TRUE_EAF_COL = "_fx_true_eaf"
FX_SPECIAL_COL = "_fx_special"
SYNTHETIC_CHROMS = ("chr1", "chr10", "chr2", "chr22", "chrX")
ROWS_PER_CHROM = 4_000
MEAN_GAP_BP = 1_000
# minus_log10_pval values placed around the significant-rows threshold (-log10 5e-8 = 7.30103).
THRESHOLD_PROBE_MLOG10P = (7.29, 7.31, 12.0)
# Chromosomes carrying the threshold probes: chrX and one autosome.
THRESHOLD_PROBE_CHROMS = ("chrX", "chr2")

FAKE_ENDPOINT = S3Endpoint(
    url="https://fake.example", region="us-east-1", credentials_profile="fake"
)
_BASES = np.array(["A", "C", "G", "T"])


def _chrom_rows(rng: np.random.Generator, chrom: str) -> pl.DataFrame:
    n = ROWS_PER_CHROM
    pos = np.cumsum(rng.integers(1, 2 * MEAN_GAP_BP, size=n)) + 10_000
    ea_idx = rng.integers(0, 4, size=n)
    oa_idx = (ea_idx + rng.integers(1, 4, size=n)) % 4
    eaf = np.round(rng.uniform(0.02, 0.98, size=n), 6)
    mlog10p = np.round(rng.exponential(1 / np.log(10), size=n), 4)
    se = np.round(rng.uniform(0.01, 0.05, size=n), 5)
    beta = np.round(rng.normal(0, 1, size=n) * se, 6)
    ea = _BASES[ea_idx]
    oa = _BASES[oa_idx]
    return pl.DataFrame(
        {
            DECODE_CHROM_COL: [chrom] * n,
            DECODE_POS_COL: pos.astype(np.int64),
            DECODE_NAME_COL: [
                f"{chrom}:{p}:{a}:{b}" for p, a, b in zip(pos, ea, oa, strict=True)
            ],
            "rsids": [f"rs{i}" if i % 5 else None for i in rng.integers(1, 10**8, n)],
            DECODE_EFFECT_ALLELE_COL: ea,
            DECODE_OTHER_ALLELE_COL: oa,
            "Beta": beta,
            DECODE_PVAL_COL: 10.0**-mlog10p,
            DECODE_MLOG10P_COL: mlog10p,
            "SE": se,
            "N": np.round(rng.uniform(30_400, 35_700, size=n)),
            DECODE_IMP_MAF_COL: np.round(np.minimum(eaf, 1 - eaf), 6),
            FX_TRUE_EAF_COL: eaf,
            FX_SPECIAL_COL: np.zeros(n, dtype=bool),
        }
    )


def _set_row(frame: pl.DataFrame, i: int, **values: object) -> pl.DataFrame:
    mask = pl.int_range(pl.len()) == i
    return frame.with_columns(
        pl.when(mask).then(pl.lit(v)).otherwise(pl.col(c)).alias(c)
        for c, v in values.items()
    )


def synthetic_rows(seed: int) -> pl.DataFrame:
    """About 20k deCODE rows (1.3 MB of text) with the special rows described above."""
    rng = np.random.default_rng(seed)
    parts = []
    for chrom in SYNTHETIC_CHROMS:
        part = _chrom_rows(rng, chrom)
        if chrom == "chr1":
            part = _set_row(
                part,
                100,
                **{
                    DECODE_OTHER_ALLELE_COL: DECODE_MISSING_ALLELE,
                    FX_SPECIAL_COL: True,
                },
            )
            part = _set_row(
                part,
                200,
                **{
                    DECODE_IMP_MAF_COL: 0.3,
                    FX_TRUE_EAF_COL: None,
                    FX_SPECIAL_COL: True,
                },
            )
        if chrom in THRESHOLD_PROBE_CHROMS:
            for k, mlog10p in enumerate(THRESHOLD_PROBE_MLOG10P):
                part = _set_row(
                    part,
                    1_000 + 500 * k,
                    **{
                        DECODE_MLOG10P_COL: mlog10p,
                        DECODE_PVAL_COL: 10.0**-mlog10p,
                        FX_SPECIAL_COL: True,
                    },
                )
        parts.append(part)
    return pl.concat(parts)


def decode_text(rows: pl.DataFrame) -> bytes:
    """The deCODE file text (header included) of rows."""
    return rows.select(DECODE_HEADER).write_csv(separator="\t", null_value=".").encode()


def write_synthetic_decode_file(path: Path, rows: pl.DataFrame) -> None:
    """Plain single-member gzip at compresslevel 1, like deCODE's files."""
    path.write_bytes(gzip.compress(decode_text(rows), compresslevel=1, mtime=0))


def synthetic_variant_info(rows: pl.DataFrame) -> pl.DataFrame:
    """
    variant_info_iceland rows (Chrom, Pos, A1, ImpFreqA1, ImpInfo) for rows: A1 = EA with the
    true EAF and A1 = other allele with 1 - EAF. The ambiguous site (null true EAF) gets
    A1 = EA twice, with f = ImpMAF and f = 1 - ImpMAF. Every 7th row also gets a disagreeing
    A1 = EA candidate (f = 0.5), as at a multi-allelic site.
    """
    base = rows.with_row_index("_i").filter(
        pl.col(DECODE_OTHER_ALLELE_COL) != DECODE_MISSING_ALLELE
    )
    eaf = pl.col(FX_TRUE_EAF_COL)
    maf = pl.col(DECODE_IMP_MAF_COL)

    def vi(a1: str | pl.Expr, f: pl.Expr, frame: pl.DataFrame) -> pl.DataFrame:
        return frame.select(
            pl.col(DECODE_CHROM_COL).alias("Chrom"),
            pl.col(DECODE_POS_COL).alias("Pos"),
            (pl.col(a1) if isinstance(a1, str) else a1).alias("A1"),
            f.alias("ImpFreqA1"),
            pl.lit(0.95).alias("ImpInfo"),
        )

    known = base.filter(eaf.is_not_null())
    ambiguous = base.filter(eaf.is_null())
    return pl.concat(
        [
            vi(DECODE_EFFECT_ALLELE_COL, eaf, known),
            vi(DECODE_OTHER_ALLELE_COL, (1 - eaf).round(6), known),
            vi(DECODE_EFFECT_ALLELE_COL, maf, ambiguous),
            vi(DECODE_EFFECT_ALLELE_COL, (1 - maf).round(6), ambiguous),
            vi(
                DECODE_EFFECT_ALLELE_COL,
                pl.lit(0.5),
                known.filter(pl.col("_i") % 7 == 0),
            ),
        ]
    )


def write_synthetic_variant_info(path: Path, rows: pl.DataFrame) -> None:
    """Gzipped tab-separated variant_info with every row written twice, as deCODE's has."""
    vi = synthetic_variant_info(rows)
    text = pl.concat([vi, vi]).write_csv(separator="\t")
    path.write_bytes(gzip.compress(text.encode(), compresslevel=1, mtime=0))


def file_md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def fake_ref(key: str) -> S3ObjectRef:
    return S3ObjectRef(endpoint=FAKE_ENDPOINT, bucket="fake-bucket", key=key)


@frozen(slots=True)
class FakeS3Access(S3Access):
    """
    objects: the local file behind each ref.
    etag_only: head reports md5=None and the md5 as the ETag (as for variant_info).
    reported_md5: md5 values head reports instead of the true ones.
    """

    objects: dict[S3ObjectRef, Path]
    etag_only: bool = False
    reported_md5: dict[S3ObjectRef, str] = field(factory=dict)

    def head(self, ref: S3ObjectRef) -> S3ObjectInfo:
        path = self.objects[ref]
        md5 = self.reported_md5.get(ref, file_md5(path))
        return S3ObjectInfo(
            size=path.stat().st_size,
            md5=None if self.etag_only else md5,
            etag=md5,
        )

    def presigned_url(self, ref: S3ObjectRef, expires_s: int) -> str:
        return self.objects[ref].resolve().as_uri()


@frozen(slots=True)
class FakeParallelRangeDownloader(ParallelRangeDownloader):
    """Copies file:// sources, then verifies the md5 like the real downloader."""

    def download(
        self,
        url: str,
        size: int,
        local_path: Path,
        md5_hash: str,
        n_connections: int,
        chunk_bytes: int,
    ) -> None:
        assert url.startswith("file://"), url
        source = Path(url.removeprefix("file://"))
        assert source.stat().st_size == size
        shutil.copyfile(source, local_path)
        verify_hash(local_path, md5_hash)


def membership_rows(rows: pl.DataFrame, every: int = 3) -> pl.DataFrame:
    """
    The rows a synthetic HapMap3-like membership list selects: every every-th autosomal row
    that is not special.
    """
    return (
        rows.with_row_index("_i")
        .filter(
            (pl.col("_i") % every == 0)
            & ~pl.col(FX_SPECIAL_COL)
            & (pl.col(DECODE_CHROM_COL) != "chrX")
        )
        .drop("_i")
    )


def synthetic_membership(selected: pl.DataFrame) -> pl.DataFrame:
    """Membership frame (CHR, POS, EA, NEA, rsID) for selected rows, alternate rows swapped."""
    swap = pl.int_range(pl.len()) % 2 == 1
    ea, oa = pl.col(DECODE_EFFECT_ALLELE_COL), pl.col(DECODE_OTHER_ALLELE_COL)
    return selected.select(
        pl.col(DECODE_CHROM_COL)
        .str.strip_prefix("chr")
        .cast(pl.Int32)
        .alias(GWASLAB_CHROM_COL),
        pl.col(DECODE_POS_COL).cast(pl.Int32).alias(GWASLAB_POS_COL),
        pl.when(swap).then(oa).otherwise(ea).alias(GWASLAB_EFFECT_ALLELE_COL),
        pl.when(swap).then(ea).otherwise(oa).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
        (pl.lit("rs") + pl.col(DECODE_POS_COL).cast(pl.String)).alias(GWASLAB_RSID_COL),
    )


def write_variant_info_asset(path: Path, rows: pl.DataFrame) -> None:
    """The parquet BuildDecodeVariantInfoTask would produce for rows."""
    source = path.with_suffix(".txt.gz")
    write_synthetic_variant_info(source, rows)
    write_variant_info_parquet(source, path)


def synthetic_protein(
    path: Path, key: str, gene: str = "GENE1"
) -> DecodeSomascanProtein:
    """A manifest row for the local file at path, with the file's true size and md5."""
    return DecodeSomascanProtein(
        seq_id=SeqId("10000-1"),
        gene_symbol=gene,
        target_name=gene,
        uniprot=UniProtId("P00001"),
        entrez_gene_id=None,
        ref=fake_ref(key),
        size=path.stat().st_size,
        md5=file_md5(path),
    )


KIB = 1024
# Small rapidgzip chunks and position-map steps give a 2 MB fixture many seek points.
TEST_CHUNK_SIZE = 32 * KIB
TEST_POSMAP_STEP = 8 * KIB
_TEMPLATE_KEY = "final_somascan_smp/template.txt.gz"


def table_meta(asset_id: str) -> ResultTableMeta:
    return ResultTableMeta(
        id=AssetId(asset_id),
        trait="t",
        project="p",
        sub_dir=Path("x"),
        extension=".parquet",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    )


@frozen(slots=True)
class DecodeDatabaseFixture:
    """
    A synthetic group: template rows, the index and variant_info assets, and the fake S3.
    assets maps asset ids to materialized assets, for a fetch.
    """

    root: Path
    rows: pl.DataFrame
    index_task: Task
    variant_info_task: Task
    assets: dict[AssetId, Asset]
    objects: dict[S3ObjectRef, Path]

    def fetch(self, asset_id: AssetId) -> Asset:
        return self.assets[asset_id]

    def wf(self, access: S3Access | None = None) -> WF:
        return make_wf(
            s3_access=FakeS3Access(objects=self.objects) if access is None else access,
            parallel_range_downloader=FakeParallelRangeDownloader(),
        )

    def protein_task(
        self, protein_rows: pl.DataFrame, key: str, md5: str | None = None
    ) -> BuildDecodeSomascanProteinTask:
        """Write protein_rows as the object at key; md5 overrides the manifest's md5."""
        path = self.root / Path(key).name
        write_synthetic_decode_file(path, protein_rows)
        ref = fake_ref(key)
        self.objects[ref] = path
        protein = synthetic_protein(path, key)
        if md5 is not None:
            protein = DecodeSomascanProtein(
                seq_id=protein.seq_id,
                gene_symbol=protein.gene_symbol,
                target_name=protein.target_name,
                uniprot=protein.uniprot,
                entrez_gene_id=protein.entrez_gene_id,
                ref=protein.ref,
                size=protein.size,
                md5=md5,
            )
        return BuildDecodeSomascanProteinTask.create(
            index_task=self.index_task,
            protein=protein,
            asset_id=f"protein_{Path(key).stem.split('.')[0]}",
            index_name="test",
            group="smp",
            parallelization=2,
            chunk_size=TEST_CHUNK_SIZE,
            posmap_step=TEST_POSMAP_STEP,
        )

    def build_protein(self, task: BuildDecodeSomascanProteinTask) -> DirectoryAsset:
        scratch = self.root / f"scratch_{task.asset_id}"
        scratch.mkdir()
        asset = task.execute(scratch, fetch=self.fetch, wf=self.wf())
        assert isinstance(asset, DirectoryAsset)
        self.assets[task.asset_id] = asset
        return asset


def decode_database_fixture(root: Path, seed: int) -> DecodeDatabaseFixture:
    """Build the variant_info parquet and the index of a synthetic template (seed)."""
    rows = synthetic_rows(seed)
    template_path = root / "template.txt.gz"
    write_synthetic_decode_file(template_path, rows)
    membership_path = root / "membership.parquet"
    synthetic_membership(membership_rows(rows)).write_parquet(membership_path)
    vi_path = root / "variant_info.parquet"
    write_variant_info_asset(vi_path, rows)
    membership_task = FakeTask(meta=table_meta("membership"))
    variant_info_task = FakeTask(meta=table_meta("variant_info"))
    assets: dict[AssetId, Asset] = {
        membership_task.asset_id: FileAsset(membership_path),
        variant_info_task.asset_id: FileAsset(vi_path),
    }
    template = synthetic_protein(template_path, _TEMPLATE_KEY)
    objects = {template.ref: template_path}
    index_task = ConstructDecodeVariantIndexTask.create(
        "index",
        template,
        membership_task=membership_task,
        variant_info_task=variant_info_task,
    )
    scratch = root / "scratch_index"
    scratch.mkdir()
    index_asset = index_task.execute(
        scratch,
        fetch=lambda asset_id: assets[asset_id],
        wf=make_wf(
            s3_access=FakeS3Access(objects=objects),
            parallel_range_downloader=FakeParallelRangeDownloader(),
        ),
    )
    assets[index_task.asset_id] = index_asset
    return DecodeDatabaseFixture(
        root=root,
        rows=rows,
        index_task=index_task,
        variant_info_task=variant_info_task,
        assets=assets,
        objects=objects,
    )


@frozen(slots=True)
class CountingRangeSource(RangeSource):
    """Reads ranges of a local file:// object, recording each (start, length) in reads."""

    path: Path
    reads: list[tuple[int, int]]

    def read_range(self, start: int, length: int) -> bytes:
        self.reads.append((start, length))
        with open(self.path, "rb") as f:
            f.seek(start)
            return f.read(length)


@frozen(slots=True)
class CountingRangeSourceFactory(RangeSourceFactory):
    """Opens file:// URLs as CountingRangeSources sharing one reads list."""

    reads: list[tuple[int, int]] = field(factory=list)

    def open(self, url: str, size: int) -> RangeSource:
        assert url.startswith("file://"), url
        path = Path(url.removeprefix("file://"))
        assert path.stat().st_size == size
        return CountingRangeSource(path=path, reads=self.reads)
