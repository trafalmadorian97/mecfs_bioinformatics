"""Task-level tests of GnomadAlleleFrequencyPanelTask on pre-written per-chromosome parts."""

from pathlib import Path, PurePath

import polars as pl
import pyarrow.parquet as pq
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.fasta_meta import FASTAMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_allele_frequency_panel_task import (
    GnomadAlleleFrequencyPanelTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_chromosome_allele_frequency_task import (
    GnomadChromosomeAlleleFrequencyTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_release import (
    GnomadRelease,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.allele_frequency_panel_constants import (
    panel_af_col,
    panel_an_col,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL

_RELEASE = GnomadRelease(
    name="test_release",
    build="19",
    vcf_url_template="unused/sites.{chrom}.vcf.bgz",
    contig_prefix="",
    chromosomes=(1, 23),
    main_groups=("afr", "nfe"),
    extra_groups=(),
    header_assembly="gnomAD_GRCh37",
)
_FASTA = FakeTask(
    FASTAMeta(
        group="genome_sequence",
        sub_group="synthetic",
        sub_folder=PurePath("processed"),
        id=AssetId("fasta"),
        build="19",
    )
)


def _part(chrom: int, positions: list[int], extra_column: bool = False) -> pl.DataFrame:
    n = len(positions)
    frame = pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [chrom] * n,
            GWASLAB_POS_COL: positions,
            "REF": ["A"] * n,
            "ALT": ["G"] * n,
            panel_af_col("afr"): [0.1] * n,
            panel_an_col("afr"): [100] * n,
            panel_af_col("nfe"): [0.2] * n,
            panel_an_col("nfe"): [200] * n,
        },
        schema={
            GWASLAB_CHROM_COL: pl.Int32,
            GWASLAB_POS_COL: pl.Int32,
            "REF": pl.String,
            "ALT": pl.String,
            panel_af_col("afr"): pl.Float32,
            panel_an_col("afr"): pl.Int32,
            panel_af_col("nfe"): pl.Float32,
            panel_an_col("nfe"): pl.Int32,
        },
    )
    if extra_column:
        frame = frame.with_columns(pl.lit(0).alias("UNEXPECTED"))
    return frame


def _part_task(chrom: int) -> GnomadChromosomeAlleleFrequencyTask:
    return GnomadChromosomeAlleleFrequencyTask.create(
        release=_RELEASE, chrom=chrom, fasta_task=_FASTA
    )


def _panel_task() -> GnomadAlleleFrequencyPanelTask:
    return GnomadAlleleFrequencyPanelTask.create(
        asset_id="gnomad_panel",
        release=_RELEASE,
        part_tasks=[_part_task(chrom) for chrom in _RELEASE.chromosomes],
    )


def _run(tmp_path: Path, parts: dict[int, pl.DataFrame]) -> Path:
    task = _panel_task()
    assets: dict[str, Asset] = {}
    for part_task in task.part_tasks:
        path = tmp_path / f"{part_task.asset_id}.parquet"
        parts[part_task.chrom].write_parquet(path)
        assets[part_task.asset_id] = FileAsset(path)

    def fetch(asset_id: AssetId) -> Asset:
        return assets[asset_id]

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    return result.path


def test_parts_are_concatenated_in_chromosome_order(tmp_path: Path) -> None:
    path = _run(tmp_path, {23: _part(23, [5]), 1: _part(1, [3, 7])})
    panel = pl.read_parquet(path)
    assert panel.select(GWASLAB_CHROM_COL, GWASLAB_POS_COL).rows() == [
        (1, 3),
        (1, 7),
        (23, 5),
    ]
    columns = _panel_task().meta.allele_frequency_columns
    assert columns is not None
    assert columns.column_for("nfe") in panel.columns
    row_group = pq.ParquetFile(path).metadata.row_group(0)
    encodings = {
        row_group.column(i).path_in_schema: set(row_group.column(i).encodings)
        for i in range(row_group.num_columns)
    }
    assert "BYTE_STREAM_SPLIT" in encodings[panel_an_col("nfe")]
    assert "RLE_DICTIONARY" in encodings[panel_af_col("nfe")]


def test_parts_with_different_schemas_fail(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        _run(tmp_path, {1: _part(1, [3]), 23: _part(23, [5], extra_column=True)})


@pytest.mark.parametrize("chromosomes", [[1], [23, 1], [1, 23, 23]])
def test_parts_not_matching_the_release_chromosomes_are_rejected(
    chromosomes: list[int],
) -> None:
    # A missing, reordered or repeated part would give an incomplete or unsorted panel.
    with pytest.raises(AssertionError):
        GnomadAlleleFrequencyPanelTask.create(
            asset_id="gnomad_panel",
            release=_RELEASE,
            part_tasks=[_part_task(chrom) for chrom in chromosomes],
        )
