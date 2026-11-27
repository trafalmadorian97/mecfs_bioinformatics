from pathlib import Path

import polars as pl
import pytest

from mecfs_bio.build_system.task.decode_database.build_decode_somascan_protein_task import (
    BuildDecodeSomascanProteinTask,
)
from mecfs_bio.build_system.task.decode_database.decode_sumstats_format import (
    chrom_code_expr,
)
from mecfs_bio.build_system.task.protein_coloc.decode_window_reader import (
    DecodeProteinEntry,
    DecodeWindowReader,
    RemoteObjectChanged,
)
from mecfs_bio.build_system.task.protein_coloc.dense_window_fetch import (
    ParallelWindowFetcher,
    ProteinWindowRequest,
)
from mecfs_bio.build_system.wf.s3_access import S3Access, S3ObjectInfo, S3ObjectRef
from mecfs_bio.constants.decode_somascan_constants import (
    DECODE_BETA_COL,
    DECODE_CHROM_COL,
    DECODE_EFFECT_ALLELE_COL,
    DECODE_IMP_MAF_COL,
    DECODE_MISSING_ALLELE,
    DECODE_N_COL,
    DECODE_OTHER_ALLELE_COL,
    DECODE_POS_COL,
    DECODE_SE_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.protein_coloc_constants import WINDOW_SCHEMA
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval
from test_mecfs_bio.unit.build_system.task.decode_database.decode_fixtures import (
    FX_CLEAN_COL,
    FX_TRUE_EAF_COL,
    CountingRangeSourceFactory,
    DecodeDatabaseFixture,
    FakeS3Access,
    decode_database_fixture,
)

_KEY = "final_somascan_smp/protein_w.txt.gz"
_BIG_PIECE = 1 << 30


@pytest.fixture(scope="module")
def database(tmp_path_factory: pytest.TempPathFactory) -> DecodeDatabaseFixture:
    return decode_database_fixture(tmp_path_factory.mktemp("decode_windows"), seed=5)


@pytest.fixture(scope="module")
def built(database: DecodeDatabaseFixture) -> BuildDecodeSomascanProteinTask:
    task = database.protein_task(database.rows, _KEY)
    database.build_protein(task)
    return task


def _reader(
    database: DecodeDatabaseFixture,
    task: BuildDecodeSomascanProteinTask,
    factory: CountingRangeSourceFactory,
) -> DecodeWindowReader:
    return DecodeWindowReader(
        proteins=(DecodeProteinEntry(task),),
        clean_variant_info_task=database.clean_variant_info_task,
        range_source_factory=factory,
        prefetch_connections=2,
        prefetch_piece_bytes=_BIG_PIECE,
    )


def _read(
    database: DecodeDatabaseFixture,
    task: BuildDecodeSomascanProteinTask,
    windows: list[GenomicInterval],
    tmp_path: Path,
    factory: CountingRangeSourceFactory | None = None,
) -> dict[GenomicInterval, pl.DataFrame]:
    reader = _reader(database, task, factory or CountingRangeSourceFactory())
    return reader.fetch_protein(
        str(task.protein.seq_id), windows, tmp_path, database.fetch, database.wf()
    )


def _pos(database: DecodeDatabaseFixture, chrom: str, i: int) -> int:
    return database.rows.filter(pl.col(DECODE_CHROM_COL) == chrom)[DECODE_POS_COL][i]


def _n_clean(database: DecodeDatabaseFixture, chrom: str, lo: int, hi: int) -> int:
    """Clean rows among the chromosome's rows lo..hi-1."""
    rows = database.rows.filter(pl.col(DECODE_CHROM_COL) == chrom)[lo:hi]
    return int(rows[FX_CLEAN_COL].sum())


def _brute_force(rows: pl.DataFrame, window: GenomicInterval) -> pl.DataFrame:
    return (
        rows.with_columns(chrom_code_expr(DECODE_CHROM_COL).alias("_code"))
        .filter(
            (pl.col("_code") == window.chrom)
            & pl.col(DECODE_POS_COL).is_between(window.start, window.end)
            & (pl.col(DECODE_OTHER_ALLELE_COL) != DECODE_MISSING_ALLELE)
            & pl.col(FX_CLEAN_COL)
        )
        .select(
            pl.col("_code").alias(GWASLAB_CHROM_COL),
            pl.col(DECODE_POS_COL).alias(GWASLAB_POS_COL),
            pl.col(DECODE_EFFECT_ALLELE_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
            pl.col(DECODE_OTHER_ALLELE_COL).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
            pl.col(DECODE_BETA_COL).alias(GWASLAB_BETA_COL),
            pl.col(DECODE_SE_COL).alias(GWASLAB_SE_COL),
            pl.col(FX_TRUE_EAF_COL).alias(GWASLAB_EFFECT_ALLELE_FREQ_COL),
            pl.col(DECODE_N_COL).alias(GWASLAB_SAMPLE_SIZE_COLUMN),
        )
        .cast(WINDOW_SCHEMA)
    )


def _assert_matches_brute_force(
    database: DecodeDatabaseFixture,
    frames: dict[GenomicInterval, pl.DataFrame],
    windows: list[GenomicInterval],
) -> None:
    assert set(frames) == set(windows)
    for window in windows:
        expected = _brute_force(database.rows, window)
        assert frames[window].schema == WINDOW_SCHEMA
        assert frames[window].equals(expected), window


def test_windows_match_brute_force(
    database: DecodeDatabaseFixture,
    built: BuildDecodeSomascanProteinTask,
    tmp_path: Path,
) -> None:
    windows = [
        GenomicInterval(1, 1, _pos(database, "chr1", 250)),
        GenomicInterval(23, _pos(database, "chrX", 3_900), 10**9),
        # About 100 KB of text: spans several seek points.
        GenomicInterval(
            10, _pos(database, "chr10", 1_000), _pos(database, "chr10", 2_500)
        ),
        GenomicInterval(22, _pos(database, "chr22", 10), _pos(database, "chr22", 20)),
        GenomicInterval(
            23, _pos(database, "chrX", 1_000), _pos(database, "chrX", 1_100)
        ),
        # No chromosome 3 in the file.
        GenomicInterval(3, 1, 10**9),
    ]
    frames = _read(database, built, windows, tmp_path)
    assert frames[windows[0]].height == _n_clean(database, "chr1", 0, 251) - 1
    assert frames[windows[-1]].height == 0
    _assert_matches_brute_force(database, frames, windows)


def test_window_at_text_order_chromosome_boundary(
    database: DecodeDatabaseFixture,
    built: BuildDecodeSomascanProteinTask,
    tmp_path: Path,
) -> None:
    # chr10 is followed by chr2 in the file, and chr1 by chr10.
    windows = [
        GenomicInterval(10, _pos(database, "chr10", 3_990), 10**9),
        GenomicInterval(2, 1, _pos(database, "chr2", 20)),
        GenomicInterval(1, _pos(database, "chr1", 3_995), 10**9),
    ]
    frames = _read(database, built, windows, tmp_path)
    assert frames[windows[0]].height == _n_clean(database, "chr10", 3_990, 4_000)
    assert frames[windows[1]].height == _n_clean(database, "chr2", 0, 21)
    _assert_matches_brute_force(database, frames, windows)


def test_prefetch_covers_reads(
    database: DecodeDatabaseFixture,
    built: BuildDecodeSomascanProteinTask,
    tmp_path: Path,
) -> None:
    # With prefetch pieces larger than the file, each window's prefetch is one read; any read
    # indexed_gzip makes outside the prefetched span would be another.
    factory = CountingRangeSourceFactory()
    windows = [
        GenomicInterval(2, _pos(database, "chr2", 500), _pos(database, "chr2", 900)),
        GenomicInterval(22, _pos(database, "chr22", 3_000), 10**9),
    ]
    _read(database, built, windows, tmp_path, factory)
    assert len(factory.reads) == len(windows)


def test_eaf_null_where_rule_undetermined(
    database: DecodeDatabaseFixture,
    built: BuildDecodeSomascanProteinTask,
    tmp_path: Path,
) -> None:
    window = GenomicInterval(
        1, _pos(database, "chr1", 150), _pos(database, "chr1", 250)
    )
    frame = _read(database, built, [window], tmp_path)[window]
    ambiguous_pos = _pos(database, "chr1", 200)
    eaf = pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL)
    assert (
        frame.filter(pl.col(GWASLAB_POS_COL) == ambiguous_pos)[
            GWASLAB_EFFECT_ALLELE_FREQ_COL
        ]
        .is_null()
        .all()
    )
    source = database.rows.filter(
        (pl.col(DECODE_CHROM_COL) == "chr1")
        & pl.col(DECODE_POS_COL).is_between(window.start, window.end)
    ).select(
        pl.col(DECODE_POS_COL).alias(GWASLAB_POS_COL),
        pl.col(DECODE_EFFECT_ALLELE_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
        pl.col(DECODE_IMP_MAF_COL),
    )
    major = frame.join(source, on=[GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL]).filter(
        eaf > 0.5
    )
    assert major.height > 10
    assert (major[GWASLAB_EFFECT_ALLELE_FREQ_COL] != major[DECODE_IMP_MAF_COL]).all()


class _CountingHeads(S3Access):
    def __init__(self, inner: S3Access) -> None:
        self.inner = inner
        self.heads = 0

    def head(self, ref: S3ObjectRef) -> S3ObjectInfo:
        self.heads += 1
        return self.inner.head(ref)

    def presigned_url(self, ref: S3ObjectRef, expires_s: int) -> str:
        return self.inner.presigned_url(ref, expires_s)


def test_md5_mismatch_is_not_retried(
    database: DecodeDatabaseFixture,
    built: BuildDecodeSomascanProteinTask,
    tmp_path: Path,
) -> None:
    access = _CountingHeads(
        FakeS3Access(
            objects=database.objects, reported_md5={built.protein.ref: "0" * 32}
        )
    )
    reader = _reader(database, built, CountingRangeSourceFactory())
    request = ProteinWindowRequest(
        protein_id=str(built.protein.seq_id), windows=(GenomicInterval(1, 1, 50_000),)
    )
    results = ParallelWindowFetcher(reader, max_parallel=1).fetch_all(
        [request], tmp_path, database.fetch, database.wf(access)
    )
    with pytest.raises(RemoteObjectChanged):
        next(results)
    assert access.heads == 1
