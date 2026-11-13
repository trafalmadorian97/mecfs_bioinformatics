from pathlib import Path

import polars as pl
import pytest

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.protein_coloc.coloc_display_table_task import (
    ColocDisplayTableTask,
    ranked_display_loci,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.protein_coloc_constants import (
    DISPLAY_FDR_SIGNIFICANT_COL,
    DISPLAY_GENE_COL,
    DISPLAY_LEAD_VARIANT_COL,
    DISPLAY_LOCUS_COL,
    DISPLAY_N_SNP_COL,
    DISPLAY_PP_H4_COL,
    DISPLAY_PROTEIN_P_COL,
    DISPLAY_Q_VALUE_COL,
    DISPLAY_SAME_SIGN_COL,
    DISPLAY_TRAIT_P_COL,
)
from test_mecfs_bio.unit.build_system.task.protein_coloc.coloc_result_fixtures import (
    COLOC_ID,
    FIXTURE_LEAD_MLOG10P,
    FIXTURE_N_ALIGNED,
    FIXTURE_PROTEIN_LEAD_MLOG10P,
    coloc_task,
    write_coloc_dir,
)


def _display(
    tmp_path: Path, h4_threshold: float, fdr_alpha: float = 0.05
) -> pl.DataFrame:
    coloc_dir = write_coloc_dir(tmp_path / "coloc")
    task = ColocDisplayTableTask.create(
        asset_id="display",
        coloc_task=coloc_task(),
        h4_threshold=h4_threshold,
        fdr_alpha=fdr_alpha,
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    assets = {AssetId(COLOC_ID): DirectoryAsset(coloc_dir)}
    result = task.execute(
        scratch, fetch=lambda asset_id: assets[asset_id], wf=make_wf()
    )
    assert isinstance(result, FileAsset)
    return pl.read_parquet(result.path)


def test_display_keeps_only_rows_above_threshold_sorted(tmp_path: Path) -> None:
    display = _display(tmp_path, h4_threshold=0.8)
    assert display[DISPLAY_GENE_COL].to_list() == ["AAA", "BBB"]
    assert display[DISPLAY_PP_H4_COL].to_list() == [0.95, 0.85]
    assert display[DISPLAY_LOCUS_COL].to_list() == ["chr1:1000-3000"] * 2
    assert display[DISPLAY_LEAD_VARIANT_COL].to_list() == ["chr1:2000:G:A"] * 2


def test_display_integer_columns_are_int32(tmp_path: Path) -> None:
    display = _display(tmp_path, h4_threshold=0.8)
    assert pl.Int64 not in display.schema.dtypes()


def test_empty_display_when_nothing_passes(tmp_path: Path) -> None:
    display = _display(tmp_path, h4_threshold=0.99)
    assert display.height == 0
    assert DISPLAY_GENE_COL in display.columns and DISPLAY_PP_H4_COL in display.columns


def test_q_value_is_running_mean_of_one_minus_pp_h4(tmp_path: Path) -> None:
    """Tested PP.H4 0.95, 0.85 give q = 0.05, (0.05 + 0.15) / 2; alpha 0.07 splits them."""
    display = _display(tmp_path, h4_threshold=0.5, fdr_alpha=0.07)
    assert display[DISPLAY_Q_VALUE_COL].to_list() == pytest.approx([0.05, 0.10])
    assert display[DISPLAY_FDR_SIGNIFICANT_COL].to_list() == [True, False]


def test_display_adds_lead_p_values_and_n_snp(tmp_path: Path) -> None:
    display = _display(tmp_path, h4_threshold=0.8)
    assert display[DISPLAY_TRAIT_P_COL].to_list() == pytest.approx(
        [10**-FIXTURE_LEAD_MLOG10P] * 2
    )
    assert display[DISPLAY_PROTEIN_P_COL].to_list() == pytest.approx(
        [10**-FIXTURE_PROTEIN_LEAD_MLOG10P] * 2
    )
    assert display[DISPLAY_N_SNP_COL].to_list() == [FIXTURE_N_ALIGNED] * 2


def test_same_sign_true_and_false(tmp_path: Path) -> None:
    """The fixture flips the trait betas of OID2 (BBB) relative to its protein betas."""
    display = _display(tmp_path, h4_threshold=0.8)
    assert display[DISPLAY_SAME_SIGN_COL].to_list() == [True, False]


def test_ranked_display_loci_orders_by_best_pp_h4_and_truncates() -> None:
    display = pl.DataFrame(
        {
            DISPLAY_LOCUS_COL: ["A", "A", "B", "C"],
            DISPLAY_PP_H4_COL: [0.9, 0.82, 0.95, 0.85],
        }
    )
    assert ranked_display_loci(display, max_loci=2) == ["B", "A"]
    assert ranked_display_loci(display, max_loci=10) == ["B", "A", "C"]
