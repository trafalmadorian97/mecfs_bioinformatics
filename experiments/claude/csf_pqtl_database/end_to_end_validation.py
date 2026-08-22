"""
End-to-end validation of the CSF pQTL database production code on one real aptamer.

Drives the actual Task code paths (ConstructCsfVariantIndexTask.execute,
Hapmap3MembershipTask.execute, write_slim_aptamer_parquet) against the real
GCST90421540 summary statistics cached by the storage probe, with an injected fetch
that returns local file assets. Checks the plan's expectations: the HapMap3 index has
~1,025,155 rows and the slim file is ~6.8 MB with byte-stream-split on all three
float columns.

Run: pixi r python experiments/claude/csf_pqtl_database/end_to_end_validation.py
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pyarrow.parquet as pq

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.csf_database.build_slim_aptamer_parquet_task import (
    read_aptamer_sumstats,
    write_slim_aptamer_parquet,
)
from mecfs_bio.build_system.task.csf_database.construct_csf_variant_index_task import (
    INDEX_COLUMNS,
    ConstructCsfVariantIndexTask,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.ppp_database.hapmap3_membership_task import (
    Hapmap3MembershipTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.csf_database_constants import CSF_INDEX_IS_STRAND_AMBIGUOUS_COL
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

CACHED_SUMSTATS = Path.home() / ".cache" / "csf_pqtl_probe" / "GCST90421540.tsv.gz"
EXPECTED_INDEX_ROWS = 1_025_155  # from the plan, measured on this exact file

_PARQUET_META = SimpleFileMeta(
    AssetId("x"), read_spec=DataFrameReadSpec(DataFrameParquetFormat())
)


def build_membership(scratch: Path) -> Path:
    task = Hapmap3MembershipTask.create("hapmap_3_membership_list")
    asset = task.execute(scratch_dir=scratch, fetch=lambda _: None, wf=make_wf())
    assert isinstance(asset, FileAsset)
    return asset.path


def build_index(scratch: Path, template_parquet: Path, membership_parquet: Path) -> Path:
    task = ConstructCsfVariantIndexTask.create(
        template_aptamer_task=FakeTask(SimpleFileMeta(AssetId("template"),
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()))),
        membership_task=FakeTask(SimpleFileMeta(AssetId("membership"),
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()))),
        asset_id="csf_index",
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "template":
            return FileAsset(template_parquet)
        if asset_id == "membership":
            return FileAsset(membership_parquet)
        raise ValueError(asset_id)

    asset = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(asset, FileAsset)
    return asset.path


def main() -> None:
    assert CACHED_SUMSTATS.exists(), f"missing cached sumstats: {CACHED_SUMSTATS}"
    scratch = Path("experiments/claude/csf_pqtl_database/_e2e_scratch")
    scratch.mkdir(parents=True, exist_ok=True)

    print("reading real sumstats (7.3M rows) ...")
    aptamer = read_aptamer_sumstats(CACHED_SUMSTATS)
    print(f"  {aptamer.height:,} variants (alignment columns)")

    # The real template asset (CSF_TEMPLATE_APTAMER via PipeDataFrameTask) keeps ALL
    # GWAS-SSF columns -- including effect_allele_frequency, which the index needs for
    # EAF -- so build the template parquet from a full-column read, not the thinned one.
    print("writing full-column template parquet (mirrors CSF_TEMPLATE_APTAMER) ...")
    template_parquet = scratch / "template.parquet"
    pl.read_csv(CACHED_SUMSTATS, separator="\t").write_parquet(template_parquet)

    print("building HapMap3 membership list ...")
    membership_parquet = build_membership(scratch)

    print("building CSF variant index (real template x HapMap3) ...")
    index_parquet = build_index(scratch, template_parquet, membership_parquet)
    index = pl.read_parquet(index_parquet)
    print(f"  index rows: {index.height:,} (expected {EXPECTED_INDEX_ROWS:,})")
    assert index.columns == INDEX_COLUMNS, index.columns
    assert index.height == EXPECTED_INDEX_ROWS, index.height
    # Sort key is fully deterministic on (CHR, POS, EA, NEA): re-sorting is a no-op.
    key = ["CHR", "POS", "EA", "NEA"]
    assert index.select(key).equals(index.select(key).sort(key)), "index not sorted"
    strand_ambiguous = index[CSF_INDEX_IS_STRAND_AMBIGUOUS_COL].sum()
    print(f"  strand-ambiguous variants: {strand_ambiguous:,}")

    print("aligning aptamer to index and writing slim file ...")
    slim_path = scratch / "slim.parquet.zstd"
    write_slim_aptamer_parquet(aptamer, index, slim_path)
    slim = pl.read_parquet(slim_path)
    size_mb = slim_path.stat().st_size / 1e6
    print(f"  slim rows: {slim.height:,}")
    print(f"  slim file size: {size_mb:.2f} MB (plan: ~6.8 MB)")
    assert slim.height == index.height
    assert slim.columns == [
        GWASLAB_BETA_COL,
        GWASLAB_SE_COL,
        GWASLAB_SAMPLE_SIZE_COLUMN,
    ]

    # Non-null coverage: the index is templated off THIS aptamer, so every index
    # variant must be present (no NaN) in its own slim file.
    non_null_beta = slim[GWASLAB_BETA_COL].is_not_nan().sum()
    print(f"  non-NaN beta: {non_null_beta:,} / {slim.height:,}")
    assert non_null_beta == slim.height, "self-alignment should have no missing variants"

    # N carried per variant, and varies (the whole point vs PPP's constant N).
    n_distinct = slim[GWASLAB_SAMPLE_SIZE_COLUMN].n_unique()
    print(f"  distinct N values: {n_distinct}")
    assert n_distinct > 1

    # Byte-stream-split on all three columns.
    for col in range(3):
        enc = pq.ParquetFile(slim_path).metadata.row_group(0).column(col).encodings
        assert "BYTE_STREAM_SPLIT" in enc, (col, enc)
    print("all checks passed")


if __name__ == "__main__":
    main()
