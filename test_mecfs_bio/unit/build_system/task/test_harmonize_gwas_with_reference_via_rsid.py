from pathlib import Path

import pandas as pd

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_rsid import (
    HarmonizeGWASWithReferenceViaRSIDTask,
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)


def test_harmonize_gwas_with_reference(tmp_path: Path):
    dummy_gwas = pd.DataFrame(
        [
            {
                GWASLAB_RSID_COL: "rs100",
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 1,
                GWASLAB_EFFECT_ALLELE_COL: "A",
                GWASLAB_NON_EFFECT_ALLELE_COL: "C",
                GWASLAB_BETA_COL: 0.1,  # should be flipped, since representation is reversed in reference
            },
            {
                GWASLAB_RSID_COL: "rs200",
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 2,
                GWASLAB_EFFECT_ALLELE_COL: "A",
                GWASLAB_NON_EFFECT_ALLELE_COL: "C",
                GWASLAB_BETA_COL: 0.1,  # drop, since not present in reference
            },
            {
                GWASLAB_RSID_COL: "rs300",
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 3,
                GWASLAB_EFFECT_ALLELE_COL: "A",
                GWASLAB_NON_EFFECT_ALLELE_COL: "C",
                GWASLAB_BETA_COL: 0.1,  # keep, since present in reference
            },
            {
                GWASLAB_RSID_COL: "rs400",
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 4,
                GWASLAB_EFFECT_ALLELE_COL: "A",
                GWASLAB_NON_EFFECT_ALLELE_COL: "T",
                GWASLAB_BETA_COL: 0.1,  # drop, since palindromic an no frequency info
            },
            {
                GWASLAB_RSID_COL: "rs500",
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 5,
                GWASLAB_EFFECT_ALLELE_COL: "A",
                GWASLAB_NON_EFFECT_ALLELE_COL: "C",
                GWASLAB_BETA_COL: 0.1,
            },  # drop, since neither alleles nor flipped alleles  match reference
            {
                GWASLAB_RSID_COL: "rs600",
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 6,
                GWASLAB_EFFECT_ALLELE_COL: "A",
                GWASLAB_NON_EFFECT_ALLELE_COL: "C",
                GWASLAB_BETA_COL: 0.1,
            },  # drop, position does not match reference
        ]
    )
    reference = pd.DataFrame(
        [
            {
                GWASLAB_RSID_COL: "rs100",
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 1,
                GWASLAB_EFFECT_ALLELE_COL: "C",
                GWASLAB_NON_EFFECT_ALLELE_COL: "A",
            },
            {
                GWASLAB_RSID_COL: "rs300",
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 3,
                GWASLAB_EFFECT_ALLELE_COL: "A",
                GWASLAB_NON_EFFECT_ALLELE_COL: "C",
            },
            {
                GWASLAB_RSID_COL: "rs400",
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 4,
                GWASLAB_EFFECT_ALLELE_COL: "A",
                GWASLAB_NON_EFFECT_ALLELE_COL: "T",
                GWASLAB_BETA_COL: 0.1,
            },
            {
                GWASLAB_RSID_COL: "rs500",
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 5,
                GWASLAB_EFFECT_ALLELE_COL: "A",
                GWASLAB_NON_EFFECT_ALLELE_COL: "G",
                GWASLAB_BETA_COL: 0.1,
            },
            {
                GWASLAB_RSID_COL: "rs600",
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 60,
                GWASLAB_EFFECT_ALLELE_COL: "A",
                GWASLAB_NON_EFFECT_ALLELE_COL: "C",
                GWASLAB_BETA_COL: 0.1,
            },  # drop, position does not match reference
        ]
    )
    gwas_path = tmp_path / "gwas.parquet"
    reference_path = tmp_path / "reference.parquet"

    dummy_gwas.to_parquet(gwas_path)
    reference.to_parquet(reference_path)

    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir()
    tsk = HarmonizeGWASWithReferenceViaRSIDTask(
        meta=SimpleFileMeta(AssetId("my_file")),
        gwas_data_task=FakeTask(
            SimpleFileMeta(
                AssetId("gwas"), read_spec=DataFrameReadSpec(DataFrameParquetFormat())
            ),
        ),
        reference_task=FakeTask(
            SimpleFileMeta(
                AssetId("reference"),
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
        ),
        palindrome_strategy="drop",
        chrom_pos_strategy="drop",
        mismatched_alleles_strategy="drop",
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "gwas":
            return FileAsset(gwas_path)
        if asset_id == "reference":
            return FileAsset(reference_path)
        raise ValueError(f"Unknown asset {asset_id}")

    result = tsk.execute(scratch_dir=scratch_dir, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
    result_df = pd.read_parquet(result.path)
    assert set(result_df[GWASLAB_RSID_COL].tolist()) == {"rs100", "rs300"}
    assert result_df.loc[result_df[GWASLAB_RSID_COL] == "rs100"][
        GWASLAB_BETA_COL
    ].tolist() == [-0.1]
    assert result_df.loc[result_df[GWASLAB_RSID_COL] == "rs300"][
        GWASLAB_BETA_COL
    ].tolist() == [0.1]
