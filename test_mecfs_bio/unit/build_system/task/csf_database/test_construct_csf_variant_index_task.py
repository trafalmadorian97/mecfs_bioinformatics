from pathlib import Path

import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.csf_database.construct_csf_variant_index_task import (
    INDEX_COLUMNS,
    ConstructCsfVariantIndexTask,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.csf_database_constants import (
    CSF_INDEX_IS_STRAND_AMBIGUOUS_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)


def _parquet_meta(asset_id: str) -> SimpleFileMeta:
    return SimpleFileMeta(
        AssetId(asset_id),
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    )


def test_construct_csf_variant_index_task(tmp_path: Path):
    # Template aptamer (GWAS-SSF columns), rows deliberately unsorted to prove the
    # task sorts. chr1:3000 is absent from membership and must drop.
    template_path = tmp_path / "template.parquet"
    pl.DataFrame(
        {
            "chromosome": [1, 1, 2, 1],
            "base_pair_location": [2000, 1000, 5000, 3000],
            "effect_allele": ["T", "G", "T", "C"],  # effect
            "other_allele": ["C", "A", "A", "G"],  # non-effect
            "effect_allele_frequency": [0.3, 0.2, 0.1, 0.4],
        }
    ).write_parquet(template_path)

    # Membership reference list (normalized gwaslab columns). chr1:2000 is stored in the
    # OPPOSITE orientation to the template to prove the CSF orientation is adopted; rs99
    # is absent from the template and must not appear.
    membership_path = tmp_path / "membership.parquet"
    pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [1, 1, 2, 3],
            GWASLAB_POS_COL: [1000, 2000, 5000, 9999],
            GWASLAB_EFFECT_ALLELE_COL: ["A", "C", "A", "A"],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["G", "T", "T", "G"],
            GWASLAB_RSID_COL: ["rs1", "rs2", "rs3", "rs99"],
        }
    ).write_parquet(membership_path)

    task = ConstructCsfVariantIndexTask.create(
        template_aptamer_task=FakeTask(_parquet_meta("template_task")),
        membership_task=FakeTask(_parquet_meta("membership_task")),
        asset_id="csf_index",
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "template_task":
            return FileAsset(template_path)
        if asset_id == "membership_task":
            return FileAsset(membership_path)
        raise ValueError(f"unknown asset id {asset_id}")

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)

    out = pl.read_parquet(result.path)

    # Exact column set and order; no POS_HG19.
    assert out.columns == INDEX_COLUMNS

    # chr1:3000 (absent from membership) and rs99 (absent from template) both drop; 3
    # variants survive, sorted by (CHR, POS, EA, NEA).
    assert out.height == 3
    assert out[GWASLAB_CHROM_COL].to_list() == [1, 1, 2]
    assert out[GWASLAB_POS_COL].to_list() == [1000, 2000, 5000]
    assert out[GWASLAB_RSID_COL].to_list() == ["rs1", "rs2", "rs3"]

    # Orientation is the template's (SSF effect_allele = EA), even where membership
    # disagreed: chr1:2000 must be EA=T / NEA=C, not membership's EA=C / NEA=T.
    assert out[GWASLAB_EFFECT_ALLELE_COL].to_list() == ["G", "T", "T"]
    assert out[GWASLAB_NON_EFFECT_ALLELE_COL].to_list() == ["A", "C", "A"]

    # EAF is the template's in-sample frequency, carried through (float32).
    assert out[GWASLAB_EFFECT_ALLELE_FREQ_COL].to_list() == pytest.approx(
        [0.2, 0.3, 0.1], abs=1e-6
    )

    # Only the A/T variant (chr2:5000) is strand-ambiguous.
    assert out[CSF_INDEX_IS_STRAND_AMBIGUOUS_COL].to_list() == [False, False, True]
