from pathlib import Path

import polars as pl

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.protein_coloc.coloc_pleiotropy_flag_task import (
    ColocPleiotropyFlagTask,
    near_st36_pleiotropic_smp_col,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.decode_somascan_constants import (
    ST36_CHROM_COL,
    ST36_PLEIOTROPIC_SMP_COL,
    ST36_SMP_POS_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    LOCUS_CHR_COL,
    PROTEIN_ID_COL,
    PROTEIN_LEAD_POS_COL,
    RAW_RESULTS_FILENAME,
)
from test_mecfs_bio.unit.build_system.task.protein_coloc.coloc_result_fixtures import (
    coloc_task,
)

_ST36_POS = 1_000_000
# A pQTL that ST36 lists as pleiotropic for Olink only, not for SMP.
_OLINK_ONLY_POS = 5_000_000

# (protein id, locus chromosome, protein lead position, expected flag)
_ROWS = [
    ("at", 1, _ST36_POS, True),
    ("249kb", 1, _ST36_POS + 249_000, True),
    ("251kb", 1, _ST36_POS - 251_000, False),
    ("other_chrom", 2, _ST36_POS, False),
    ("olink_only", 1, _OLINK_ONLY_POS, False),
    ("no_lead", 1, None, None),
]


def test_flag_marks_leads_near_st36_smp_rows(tmp_path: Path) -> None:
    coloc_dir = tmp_path / "coloc"
    coloc_dir.mkdir()
    pl.DataFrame(
        {
            PROTEIN_ID_COL: [r[0] for r in _ROWS],
            LOCUS_CHR_COL: [r[1] for r in _ROWS],
            PROTEIN_LEAD_POS_COL: [r[2] for r in _ROWS],
        },
        schema={
            PROTEIN_ID_COL: pl.String,
            LOCUS_CHR_COL: pl.Int32,
            PROTEIN_LEAD_POS_COL: pl.Int32,
        },
    ).write_parquet(coloc_dir / RAW_RESULTS_FILENAME)
    st36 = tmp_path / "st36.parquet"
    pl.DataFrame(
        {
            ST36_CHROM_COL: [1, 1],
            ST36_SMP_POS_COL: [_ST36_POS, _OLINK_ONLY_POS],
            ST36_PLEIOTROPIC_SMP_COL: [True, False],
        }
    ).write_parquet(st36)
    upstream = coloc_task()
    st36_task = FakeTask(
        SimpleFileMeta("st36", read_spec=DataFrameReadSpec(DataFrameParquetFormat()))
    )
    task = ColocPleiotropyFlagTask.create(
        asset_id="flagged", coloc_task=upstream, pleiotropic_task=st36_task
    )
    assets = {
        upstream.asset_id: DirectoryAsset(coloc_dir),
        st36_task.asset_id: FileAsset(st36),
    }
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(
        scratch,
        fetch=lambda asset_id: assets[asset_id],
        wf=make_wf(),
    )
    assert isinstance(result, FileAsset)
    flagged = pl.read_parquet(result.path)
    assert flagged[PROTEIN_ID_COL].to_list() == [r[0] for r in _ROWS]
    assert flagged[near_st36_pleiotropic_smp_col(250_000)].to_list() == [
        r[3] for r in _ROWS
    ]
