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
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.ppp_database.construct_ppp_variant_index_task import (
    INDEX_COLUMNS,
    ConstructPppVariantIndexTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)
from mecfs_bio.constants.ppp_index_constants import (
    PPP_INDEX_IS_STRAND_AMBIGUOUS_COL,
    PPP_INDEX_POS_HG19_COL,
)


def _parquet_meta(asset_id: str) -> SimpleFileMeta:
    return SimpleFileMeta(
        AssetId(asset_id),
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    )


def test_construct_ppp_variant_index_task(tmp_path: Path):
    # Template protein (regenie columns). ID = CHROM:POS_hg19:ALLELE0:ALLELE1:imp:v1.
    protein_path = tmp_path / "protein.parquet"
    pl.DataFrame(
        {
            "CHROM": [1, 1, 2, 1],
            "GENPOS": [1000, 2000, 5000, 3000],  # hg38 positions
            "ID": [
                "1:900:A:G:imp:v1",  # matches membership, same orientation
                "1:1900:C:T:imp:v1",  # matches membership but SWAPPED orientation
                "2:4900:A:T:imp:v1",  # matches membership, palindromic (A/T)
                "1:2900:G:C:imp:v1",  # absent from membership -> dropped
            ],
            "ALLELE0": ["A", "C", "A", "G"],  # non-effect
            "ALLELE1": ["G", "T", "T", "C"],  # effect
            "A1FREQ": [0.2, 0.3, 0.1, 0.4],
        }
    ).write_parquet(protein_path)

    # Membership reference list (normalized gwaslab-style columns). Row 2 is stored
    # in the OPPOSITE orientation to the protein to prove PPP orientation is adopted;
    # rs99 is absent from the protein and must not appear.
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

    task = ConstructPppVariantIndexTask.create(
        template_protein_task=FakeTask(_parquet_meta("protein_task")),
        membership_task=FakeTask(_parquet_meta("membership_task")),
        asset_id="ppp_index",
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "protein_task":
            return FileAsset(protein_path)
        if asset_id == "membership_task":
            return FileAsset(membership_path)
        raise ValueError(f"unknown asset id {asset_id}")

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)

    out = pl.read_parquet(result.path)

    # Exact column set and order.
    assert out.columns == INDEX_COLUMNS

    # The unmatched protein variant (chr1:3000) and unmatched membership variant
    # (rs99) are both absent; 3 variants survive, sorted by (CHR, POS, EA, NEA).
    assert out.height == 3
    assert out[GWASLAB_CHROM_COL].to_list() == [1, 1, 2]
    assert out[GWASLAB_POS_COL].to_list() == [1000, 2000, 5000]
    assert out[PPP_INDEX_POS_HG19_COL].to_list() == [900, 1900, 4900]
    assert out[GWASLAB_RSID_COL].to_list() == ["rs1", "rs2", "rs3"]

    # Orientation is the protein's (ALLELE1 = EA), even where membership disagreed:
    # chr1:2000 must be EA=T / NEA=C, not membership's EA=C / NEA=T.
    assert out[GWASLAB_EFFECT_ALLELE_COL].to_list() == ["G", "T", "T"]
    assert out[GWASLAB_NON_EFFECT_ALLELE_COL].to_list() == ["A", "C", "A"]

    # EAF is the in-sample A1FREQ, carried through (float32, so compare approximately).
    assert out[GWASLAB_EFFECT_ALLELE_FREQ_COL].to_list() == pytest.approx(
        [0.2, 0.3, 0.1], abs=1e-6
    )

    # Only the A/T variant is strand-ambiguous.
    assert out[PPP_INDEX_IS_STRAND_AMBIGUOUS_COL].to_list() == [False, False, True]
