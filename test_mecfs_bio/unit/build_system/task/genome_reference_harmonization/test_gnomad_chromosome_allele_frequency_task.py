"""Task-level tests of GnomadChromosomeAlleleFrequencyTask on tiny local bgzipped VCFs."""

from pathlib import Path, PurePath
from subprocess import CalledProcessError

import attrs
import polars as pl
import pyarrow.parquet as pq
import pysam
import pytest
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.fasta_meta import FASTAMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_chromosome_allele_frequency_task import (
    GnomadChromosomeAlleleFrequencyTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_release import (
    GnomadRelease,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.panel_batch_checks import (
    DEFAULT_BATCH_ROWS,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.allele_frequency_panel_constants import (
    GnomadGroup,
    panel_af_col,
    panel_an_col,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL
from test_mecfs_bio.unit.build_system.task.genome_reference_harmonization.genome_reference_fixtures import (
    write_fasta,
)

# chr1: ACGT repeated (pos p holds "ACGT"[(p - 1) % 4]), with N at 101-104.
_CHR1 = "ACGT" * 25 + "NNNN" + "ACGT" * 24
_CHRX = "TTGCA" * 20
_MAIN_GROUPS: tuple[GnomadGroup, ...] = ("afr", "nfe")
_EXTRA_GROUPS: tuple[GnomadGroup, ...] = ("nfe_nwe",)
_ALL_GROUPS = _MAIN_GROUPS + _EXTRA_GROUPS
_ASSEMBLY = "gnomAD_GRCh37"
_FASTA_ID = "fasta"
# One row per batch, so every adjacent pair of records straddles a batch boundary.
_ONE_ROW_BATCH = 1
_LONG_DELETION_REF = _CHR1[8:38]  # 30 bases from pos 9


@frozen(slots=True)
class Record:
    pos: int
    ref: str
    alt: str
    filter: str = "PASS"
    afr: str = "0"
    afr_an: int = 100
    nfe: str = "0"
    nfe_an: int = 200
    nfe_nwe: str = "0"
    nfe_nwe_an: int = 100


def _header(groups: tuple[GnomadGroup, ...] = _ALL_GROUPS) -> list[str]:
    lines = [
        "##fileformat=VCFv4.2",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=RF,Description="Failed random forest">',
        f"##contig=<ID=1,length={len(_CHR1)},assembly={_ASSEMBLY}>",
        f"##contig=<ID=X,length={len(_CHRX)},assembly={_ASSEMBLY}>",
    ]
    for group in groups:
        lines.append(
            f'##INFO=<ID={panel_af_col(group)},Number=A,Type=Float,Description="AF {group}">'
        )
        lines.append(
            f'##INFO=<ID={panel_an_col(group)},Number=1,Type=Integer,Description="AN {group}">'
        )
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    return lines


def _line(record: Record) -> str:
    info = ";".join(
        [
            f"AF_afr={record.afr};AN_afr={record.afr_an}",
            f"AF_nfe={record.nfe};AN_nfe={record.nfe_an}",
            f"AF_nfe_nwe={record.nfe_nwe};AN_nfe_nwe={record.nfe_nwe_an}",
        ]
    )
    return f"1\t{record.pos}\t.\t{record.ref}\t{record.alt}\t.\t{record.filter}\t{info}"


def _release(tmp_path: Path) -> GnomadRelease:
    return GnomadRelease(
        name="test_release",
        build="19",
        vcf_url_template=str(tmp_path / "sites.{chrom}.vcf.bgz"),
        contig_prefix="",
        chromosomes=(1, 23),
        main_groups=_MAIN_GROUPS,
        extra_groups=_EXTRA_GROUPS,
        header_assembly=_ASSEMBLY,
    )


def _run(
    tmp_path: Path,
    records: list[Record],
    header: list[str] | None = None,
    batch_rows: int = DEFAULT_BATCH_ROWS,
    truncate: bool = False,
) -> Path:
    plain = tmp_path / "sites.1.vcf"
    plain.write_text(
        "\n".join((header or _header()) + [_line(r) for r in records]) + "\n"
    )
    vcf = tmp_path / "sites.1.vcf.bgz"
    pysam.tabix_compress(str(plain), str(vcf), force=True)
    if truncate:
        vcf.write_bytes(vcf.read_bytes()[:-40])
    fasta_dir = tmp_path / "fasta"
    write_fasta(fasta_dir, {"chr1": _CHR1, "chrX": _CHRX})
    fasta_task = FakeTask(
        FASTAMeta(
            group="genome_sequence",
            sub_group="synthetic",
            sub_folder=PurePath("processed"),
            id=AssetId(_FASTA_ID),
            build="19",
        )
    )
    task = attrs.evolve(
        GnomadChromosomeAlleleFrequencyTask.create(
            release=_release(tmp_path), chrom=1, fasta_task=fasta_task
        ),
        batch_rows=batch_rows,
        max_attempts=1,
    )

    def fetch(asset_id: AssetId) -> Asset:
        assert asset_id == _FASTA_ID
        return DirectoryAsset(fasta_dir)

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    return result.path


def _encodings(path: Path, column: str) -> set[str]:
    row_group = pq.ParquetFile(path).metadata.row_group(0)
    columns = [row_group.column(i) for i in range(row_group.num_columns)]
    return set(next(c for c in columns if c.path_in_schema == column).encodings)


_HAPPY_RECORDS = [
    Record(1, "A", "G", afr="0.1"),  # kept
    Record(1, "A", "T", nfe="0.01"),  # same site, other ALT: kept, not a duplicate
    Record(2, "C", "T", filter="RF", afr="0.1"),  # not PASS: dropped
    Record(3, "G", "A", nfe_nwe="0.01"),  # polymorphic only in an extra group: dropped
    Record(4, "T", "C", afr=".", afr_an=0, nfe="0.2"),  # kept, AF_afr null (AN 0)
    Record(5, "A", "AC", nfe="0.3"),  # insertion: kept
    Record(9, _LONG_DELETION_REF, "A", nfe="0.05"),  # 30-base deletion: kept
    Record(102, "C", "G", afr="0.1"),  # REF over N: dropped and counted, not fatal
]


@pytest.mark.parametrize("batch_rows", [DEFAULT_BATCH_ROWS, _ONE_ROW_BATCH])
def test_polymorphic_pass_records_become_an_allele_frequency_table(
    tmp_path: Path, batch_rows: int
) -> None:
    path = _run(tmp_path, _HAPPY_RECORDS, batch_rows=batch_rows)
    table = pl.read_parquet(path)
    assert table.select(GWASLAB_POS_COL, "ALT").rows() == [
        (1, "G"),
        (1, "T"),
        (4, "C"),
        (5, "AC"),
        (9, "A"),
    ]
    assert table[GWASLAB_CHROM_COL].unique().to_list() == [1]
    assert table.schema[GWASLAB_CHROM_COL] == pl.Int32
    assert table.schema[panel_af_col("afr")] == pl.Float32
    assert table.schema[panel_an_col("afr")] == pl.Int32
    at_4 = table.filter(pl.col(GWASLAB_POS_COL) == 4).row(0, named=True)
    assert at_4[panel_af_col("afr")] is None
    assert at_4[panel_an_col("afr")] == 0
    for group in _ALL_GROUPS:
        assert "BYTE_STREAM_SPLIT" in _encodings(path, panel_an_col(group))
        assert "RLE_DICTIONARY" in _encodings(path, panel_af_col(group))


def test_ref_mismatch_over_acgt_reference_fails(tmp_path: Path) -> None:
    # pos 6 holds C
    with pytest.raises(AssertionError):
        _run(tmp_path, [Record(6, "A", "G", afr="0.1")])


def test_header_missing_a_group_field_fails(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        _run(
            tmp_path,
            [Record(1, "A", "G", afr="0.1")],
            header=_header(groups=("afr", "nfe")),
        )


@pytest.mark.parametrize("batch_rows", [DEFAULT_BATCH_ROWS, _ONE_ROW_BATCH])
def test_duplicate_key_fails(tmp_path: Path, batch_rows: int) -> None:
    records = [Record(1, "A", "G", afr="0.1"), Record(1, "A", "G", afr="0.2")]
    with pytest.raises(AssertionError):
        _run(tmp_path, records, batch_rows=batch_rows)


def test_decreasing_position_fails(tmp_path: Path) -> None:
    records = [Record(5, "A", "G", afr="0.1"), Record(1, "A", "G", afr="0.1")]
    with pytest.raises(AssertionError):
        _run(tmp_path, records)


def test_chromosome_without_polymorphic_pass_records_fails(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        _run(tmp_path, [Record(1, "A", "G", filter="RF", afr="0.1")])


def test_truncated_vcf_fails(tmp_path: Path) -> None:
    with pytest.raises(CalledProcessError):
        _run(tmp_path, _HAPPY_RECORDS, truncate=True)
