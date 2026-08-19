from pathlib import Path

import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.sbayesrc.gctb_gwfm_constants import (
    COJO_MA_COLUMNS,
)
from mecfs_bio.build_system.task.sbayesrc.sumstats_to_cojo_ma_task import (
    SumstatsToCojoMaTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_P_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)


def test_sumstats_to_cojo_ma_task_maps_columns_and_drops_extras(tmp_path: Path):
    src_df = pl.DataFrame(
        {
            GWASLAB_RSID_COL: ["rs1", "rs2", "rs3"],
            GWASLAB_EFFECT_ALLELE_COL: ["A", "C", "G"],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["T", "G", "C"],
            GWASLAB_EFFECT_ALLELE_FREQ_COL: [0.1, 0.2, 0.3],
            GWASLAB_BETA_COL: [0.01, -0.02, 0.03],
            GWASLAB_SE_COL: [0.001, 0.002, 0.003],
            GWASLAB_P_COL: [0.5, 0.4, 0.3],
            GWASLAB_SAMPLE_SIZE_COLUMN: [1000, 1001, 1002],
            # Extra column not part of the COJO .ma spec; selection must drop it.
            "CHR": [1, 2, 3],
        }
    )
    src_parquet_path = tmp_path / "src.parquet"
    src_df.write_parquet(src_parquet_path)

    source_meta = FilteredGWASDataMeta(
        id=AssetId("src"),
        trait="t",
        project="p",
        sub_dir="x",
        read_spec=DataFrameReadSpec(format=DataFrameParquetFormat()),
        extension=".parquet",
    )
    stub_dep = FakeTask(meta=source_meta)

    task = SumstatsToCojoMaTask.create(id="out_ma", sumstats_task=stub_dep)

    def fetch(asset_id: AssetId) -> Asset:
        assert asset_id == stub_dep.asset_id
        return FileAsset(src_parquet_path)

    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir()
    result = task.execute(scratch_dir=scratch_dir, fetch=fetch, wf=make_wf())

    assert isinstance(result, FileAsset)
    result_df = pl.read_csv(result.path, separator="\t")

    assert result_df.columns == list(COJO_MA_COLUMNS)
    assert result_df["SNP"].to_list() == ["rs1", "rs2", "rs3"]
    assert result_df["A1"].to_list() == ["A", "C", "G"]
    assert result_df["A2"].to_list() == ["T", "G", "C"]
    assert result_df["freq"].to_list() == [0.1, 0.2, 0.3]
    assert result_df["b"].to_list() == [0.01, -0.02, 0.03]
    assert result_df["se"].to_list() == [0.001, 0.002, 0.003]
    assert result_df["p"].to_list() == [0.5, 0.4, 0.3]
    assert result_df["N"].to_list() == [1000, 1001, 1002]
