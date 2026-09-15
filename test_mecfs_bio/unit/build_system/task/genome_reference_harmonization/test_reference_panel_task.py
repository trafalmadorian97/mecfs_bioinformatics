"""Task-level test of ReferencePanelAlleleFrequencyTask on a tiny VCF."""

from pathlib import Path, PurePath

import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_COLUMNS,
    PANEL_REF_COL,
    ReferencePanelAlleleFrequencyTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL

_VCF_ID = "panel_vcf"
_VCF = (
    "##fileformat=VCFv4.2\n"
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n'
    "##contig=<ID=1>\n##contig=<ID=X>\n##contig=<ID=GL000191.1>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    "1\t200\t.\tTG\tT\t.\t.\tAF=0.86\tGT\t0|1\n"
    "1\t100\t.\tA\tG\t.\t.\tAF=0.2\tGT\t0|0\n"
    "1\t200\t.\tT\tTG\t.\t.\tAF=0\tGT\t0|0\n"  # mirrored record at the same site: kept
    "1\t300\t.\tC\tT\t.\t.\tAF=0.1\tGT\t0|0\n"
    "1\t300\t.\tC\tT\t.\t.\tAF=0.1\tGT\t0|0\n"  # exact duplicate: collapsed
    "1\t400\t.\tG\tA\t.\t.\tAF=0.3\tGT\t0|0\n"
    "1\t400\t.\tG\tA\t.\t.\tAF=0.4\tGT\t0|0\n"  # conflicting duplicate: removed
    "1\t500\t.\tA\t.\t.\t.\tAF=.\tGT\t0|0\n"  # monomorphic, no ALT: removed
    "X\t50\t.\tC\tG\t.\t.\tAF=0.5\tGT\t0|1\n"
    "GL000191.1\t10\t.\tA\tC\t.\t.\tAF=0.5\tGT\t0|1\n"  # non-main contig: removed
)


def _run(tmp_path: Path) -> pl.DataFrame:
    vcf_path = tmp_path / "panel.vcf"
    vcf_path.write_text(_VCF)
    source = FakeTask(
        ReferenceFileMeta(
            group="thousand_genomes",
            sub_group="synthetic",
            sub_folder=PurePath("raw"),
            extension=".vcf",
            id=AssetId(_VCF_ID),
        )
    )
    task = ReferencePanelAlleleFrequencyTask.create(vcf_task=source, asset_id="panel_af")

    def fetch(asset_id: AssetId) -> Asset:
        assert asset_id == _VCF_ID
        return FileAsset(vcf_path)

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    return pl.read_parquet(result.path)


def test_panel_keeps_main_contig_sites_with_unambiguous_frequencies(tmp_path: Path) -> None:
    panel = _run(tmp_path)
    assert panel.columns == PANEL_COLUMNS
    assert panel.rows() == [
        (1, 100, "A", "G", pytest.approx(0.2)),
        (1, 200, "T", "TG", pytest.approx(0.0)),
        (1, 200, "TG", "T", pytest.approx(0.86)),
        (1, 300, "C", "T", pytest.approx(0.1)),
        (23, 50, "C", "G", pytest.approx(0.5)),
    ]
    assert panel.schema[GWASLAB_CHROM_COL] == pl.Int32
    assert panel.schema[GWASLAB_POS_COL] == pl.Int32
    assert panel.schema[PANEL_AF_COL] == pl.Float32
    assert panel.schema[PANEL_REF_COL] == pl.String
    assert panel.schema[PANEL_ALT_COL] == pl.String
