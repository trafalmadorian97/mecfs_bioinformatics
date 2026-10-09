"""Task-level tests of PanUkbbAlleleFrequencyPanelTask on a synthetic bgzipped manifest."""

from collections.abc import Sequence
from pathlib import Path, PurePath

import attrs
import polars as pl
import pysam
import pytest
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.fasta_meta import FASTAMeta
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.genome_reference_harmonization.pan_ukbb.pan_ukbb_allele_frequency_panel_task import (
    PAN_UKBB_MANIFEST_READ_SPEC,
    PanUkbbAlleleFrequencyPanelTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.allele_frequency_panel_constants import (
    PanUkbbGroup,
    panel_af_col,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL
from test_mecfs_bio.unit.build_system.task.genome_reference_harmonization.genome_reference_fixtures import (
    write_fasta,
)

# The manifest's real header (full_variant_qc_metrics.txt.bgz, 2020-08-28).
_HEADER = (
    "chrom pos ref alt rsid varid pass_gnomad_genomes n_passing_populations high_quality "
    "nearest_genes info ac_AFR af_AFR an_AFR gnomad_genomes_ac_AFR gnomad_genomes_af_AFR "
    "gnomad_genomes_an_AFR ac_AMR af_AMR an_AMR gnomad_genomes_ac_AMR gnomad_genomes_af_AMR "
    "gnomad_genomes_an_AMR ac_CSA af_CSA an_CSA ac_EAS af_EAS an_EAS gnomad_genomes_ac_EAS "
    "gnomad_genomes_af_EAS gnomad_genomes_an_EAS ac_EUR af_EUR an_EUR gnomad_genomes_ac_EUR "
    "gnomad_genomes_af_EUR gnomad_genomes_an_EUR ac_MID af_MID an_MID"
).split()
_GROUPS: tuple[PanUkbbGroup, ...] = (
    "ukb_afr",
    "ukb_amr",
    "ukb_csa",
    "ukb_eas",
    "ukb_eur",
    "ukb_mid",
)
# chr1 cycles ACGT and ends in ten Ns (1991-2000); chrX is 200 bases.
_CHR1 = ("ACGT" * 500)[:1990] + "N" * 10
_CHRX = "TTGCA" * 40
_NEXT_BASE = {"A": "C", "C": "G", "G": "T", "T": "A"}
_SWAPPED_POSITIONS = (10, 20)
_BGZF_BLOCK_BYTES = 65536
_MANIFEST_ID = "manifest"
_FASTA_ID = "fasta"


@frozen(slots=True)
class ManifestRow:
    chrom: str
    pos: int
    ref: str
    alt: str
    af: float
    info: float = 0.95
    gnomad_af: float | None = None


def _reference_rows(chrom: str, sequence: str, n_positions: int) -> list[ManifestRow]:
    return [
        ManifestRow(
            chrom=chrom,
            pos=pos,
            ref=sequence[pos - 1],
            alt=_NEXT_BASE[sequence[pos - 1]],
            af=pos / 10_000,
        )
        for pos in range(1, n_positions + 1)
    ]


def _swapped(row: ManifestRow) -> ManifestRow:
    return attrs.evolve(row, ref=row.alt, alt=row.ref)


def _chr1_rows() -> list[ManifestRow]:
    return [
        _swapped(row) if row.pos in _SWAPPED_POSITIONS else row
        for row in _reference_rows("1", _CHR1, 1500)
    ]


def _chrx_rows() -> list[ManifestRow]:
    return _reference_rows("X", _CHRX, 100)


def _write_manifest(path: Path, rows: list[ManifestRow], header: Sequence[str]) -> None:
    lines = ["\t".join(header)]
    for row in rows:
        values = {column: "NA" for column in header} | {
            "chrom": row.chrom,
            "pos": str(row.pos),
            "ref": row.ref,
            "alt": row.alt,
            "varid": f"{row.chrom}:{row.pos}_{row.ref}_{row.alt}",
            "info": f"{row.info:.4e}",
            "gnomad_genomes_af_EUR": (
                "NA" if row.gnomad_af is None else f"{row.gnomad_af:.4e}"
            ),
            "nearest_genes": "GENE1,GENE2",
        }
        for column in header:
            if column.startswith("af_"):
                values[column] = f"{row.af:.4e}"
        lines.append("\t".join(values[column] for column in header))
    plain = path.parent / (path.name + ".plain")
    plain.write_text("\n".join(lines) + "\n")
    # Several BGZF blocks, so a reader that stops after the first gzip member is caught.
    assert plain.stat().st_size > 3 * _BGZF_BLOCK_BYTES
    pysam.tabix_compress(str(plain), str(path), force=True)


# Small batches, so the manifest arrives in many batches and chromosome changes and
# swapped rows fall inside and across batch boundaries.
_TEST_BATCH_ROWS = 97


def _run(
    tmp_path: Path,
    rows: list[ManifestRow],
    expected_ref_mismatches: int = len(_SWAPPED_POSITIONS),
    header: Sequence[str] = _HEADER,
) -> pl.DataFrame:
    manifest_path = tmp_path / "manifest.txt.bgz"
    _write_manifest(manifest_path, rows, header)
    fasta_dir = tmp_path / "fasta"
    write_fasta(fasta_dir, {"chr1": _CHR1, "chrX": _CHRX})
    created = PanUkbbAlleleFrequencyPanelTask.create(
        asset_id="pan_ukbb_panel",
        manifest_task=FakeTask(
            ReferenceFileMeta(
                group="pan_ukbb",
                sub_group="variant_manifest",
                sub_folder=PurePath("raw"),
                extension=".txt.bgz",
                id=AssetId(_MANIFEST_ID),
                read_spec=PAN_UKBB_MANIFEST_READ_SPEC,
            )
        ),
        fasta_task=FakeTask(
            FASTAMeta(
                group="genome_sequence",
                sub_group="synthetic",
                sub_folder=PurePath("processed"),
                id=AssetId(_FASTA_ID),
                build="19",
            )
        ),
        groups=_GROUPS,
        chromosomes=(1, 23),
        expected_ref_mismatches=expected_ref_mismatches,
    )
    task = attrs.evolve(created, batch_rows=_TEST_BATCH_ROWS)
    assets: dict[str, Asset] = {
        _MANIFEST_ID: FileAsset(manifest_path),
        _FASTA_ID: DirectoryAsset(fasta_dir),
    }

    def fetch(asset_id: AssetId) -> Asset:
        return assets[asset_id]

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    return pl.read_parquet(result.path)


def test_manifest_becomes_a_panel_without_the_swapped_rows(tmp_path: Path) -> None:
    rows = _chr1_rows() + _chrx_rows()
    panel = _run(tmp_path, rows)
    assert panel.height == len(rows) - len(_SWAPPED_POSITIONS)
    assert panel.columns == [
        GWASLAB_CHROM_COL,
        GWASLAB_POS_COL,
        PANEL_REF_COL,
        PANEL_ALT_COL,
        *[panel_af_col(group) for group in _GROUPS],
    ]
    assert panel.schema[panel_af_col("ukb_eur")] == pl.Float32
    assert panel[GWASLAB_CHROM_COL].unique(maintain_order=True).to_list() == [1, 23]
    chr1 = panel.filter(pl.col(GWASLAB_CHROM_COL) == 1)
    assert not set(_SWAPPED_POSITIONS) & set(chr1[GWASLAB_POS_COL].to_list())
    at_100 = chr1.filter(pl.col(GWASLAB_POS_COL) == 100).row(0, named=True)
    assert at_100[panel_af_col("ukb_eur")] == pytest.approx(0.01)


def test_swap_profile_rows_dropped_only_on_chr21_chr22_and_x(tmp_path: Path) -> None:
    genotyped_no_gnomad = {"info": 1.0, "gnomad_af": None}
    chrx = _chrx_rows()
    # A swapped genotyped chrX row without a gnomAD frequency: the profile drops it before
    # the FASTA check, so it is not a ref mismatch.
    chrx[50] = attrs.evolve(_swapped(chrx[50]), **genotyped_no_gnomad)
    # The same profile on chrX with a gnomAD frequency, and on chr1, is kept.
    chrx[60] = attrs.evolve(chrx[60], info=1.0, gnomad_af=0.2)
    chr1 = _chr1_rows()
    chr1[200] = attrs.evolve(chr1[200], **genotyped_no_gnomad)
    panel = _run(tmp_path, chr1 + chrx)
    kept = set(panel.select(GWASLAB_CHROM_COL, GWASLAB_POS_COL).iter_rows())
    assert (23, chrx[50].pos) not in kept
    assert (23, chrx[60].pos) in kept
    assert (1, chr1[200].pos) in kept
    assert panel.height == len(chr1) + len(chrx) - len(_SWAPPED_POSITIONS) - 1


def test_unexpected_number_of_ref_mismatches_fails(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        _run(
            tmp_path,
            _chr1_rows() + _chrx_rows(),
            expected_ref_mismatches=len(_SWAPPED_POSITIONS) - 1,
        )


def test_row_over_an_n_base_fails(tmp_path: Path) -> None:
    over_n = ManifestRow(chrom="1", pos=1995, ref="A", alt="C", af=0.1)
    with pytest.raises(AssertionError):
        _run(tmp_path, _chr1_rows() + [over_n] + _chrx_rows())


def test_ref_span_beyond_contig_end_fails(tmp_path: Path) -> None:
    beyond = ManifestRow(chrom="X", pos=len(_CHRX) + 1, ref="A", alt="C", af=0.1)
    with pytest.raises(AssertionError):
        _run(tmp_path, _chr1_rows() + _chrx_rows() + [beyond])


def test_missing_frequency_column_fails(tmp_path: Path) -> None:
    header = [column for column in _HEADER if column != "af_MID"]
    with pytest.raises(AssertionError):
        _run(tmp_path, _chr1_rows() + _chrx_rows(), header=header)


def test_chromosomes_out_of_order_fail(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        _run(tmp_path, _chrx_rows() + _chr1_rows())
