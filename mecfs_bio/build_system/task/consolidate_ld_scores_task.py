"""
Task to read LD scores in the standard format defined by the authors of LD score regression,
and write them out as a parquet file

Alongside the per-variant LD scores, each chromosome's .l2.M_5_50 file holds a single count: the
number of reference variants with minor allele frequency above 5% on that chromosome. LD score
regression needs the genome-wide total of those counts as the denominator that turns a regression
slope into a heritability. It is a per-chromosome scalar rather than a per-variant value, so we
carry it denormalized: an M_5_50 column that is constant within each chromosome. Consumers recover
the total with total_m_5_50. Keeping it in the same asset as the LD scores means the count and the
variants it describes can never come from different sets of chromosomes.
"""

from pathlib import Path

import narwhals
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

LD_SCORE_CHROM_COL = "CHR"
LD_SCORE_POS_COL = "BP"
LD_SCORE_RSID_COL = "SNP"
LD_SCORE_LD_SCORE_COL = "L2"
LD_SCORE_M_5_50_COL = "M_5_50"

# Standard LDSC filenames: <prefix>.l2.ldscore.gz holds the per-variant scores, and the sibling
# <prefix>.l2.M_5_50 the chromosome's common-variant count.
_LD_SCORE_SUFFIX = ".l2.ldscore.gz"
_M_5_50_SUFFIX = ".l2.M_5_50"


def read_m_5_50(m_path: Path) -> float:
    """The common-variant count from one .l2.M_5_50 file. Partitioned LD scores store one count
    per annotation on a single whitespace-separated line; their sum is the count for the whole
    chromosome, and for unpartitioned scores the file holds that single number already."""
    values = m_path.read_text().split()
    assert values, f"{m_path} is empty; expected at least one variant count"
    return float(sum(float(value) for value in values))


def total_m_5_50(ld_scores: narwhals.DataFrame) -> float:
    """The genome-wide common-variant count from a consolidated LD-score table.

    M_5_50 is constant within a chromosome, so the total is the sum over chromosomes of that
    constant -- taken over exactly the chromosomes present in the table, which is what keeps the
    count in step with the variants it normalizes."""
    per_chromosome = ld_scores.group_by(LD_SCORE_CHROM_COL).agg(
        narwhals.col(LD_SCORE_M_5_50_COL).n_unique().alias("__n_unique__"),
        narwhals.col(LD_SCORE_M_5_50_COL).max().alias(LD_SCORE_M_5_50_COL),
    )
    assert (per_chromosome["__n_unique__"] == 1).all(), (
        f"{LD_SCORE_M_5_50_COL} must be constant within each chromosome, but at least one "
        "chromosome carries several distinct counts"
    )
    return float(per_chromosome[LD_SCORE_M_5_50_COL].sum())


@frozen
class ConsolidateLDScoresTask(Task):
    """
    Task to read LD scores in the standard format defined by the authors of LD score regression,
    and write them out as a parquet file
    """

    meta: Meta
    extracted_ld_scores_task: Task

    @property
    def deps(self) -> list["Task"]:
        return [self.extracted_ld_scores_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        asset = fetch(self.extracted_ld_scores_task.asset_id)
        assert isinstance(asset, DirectoryAsset)
        frames = []
        for ld_file in sorted(asset.path.glob(f"*{_LD_SCORE_SUFFIX}")):
            # Deriving the M_5_50 path from the score file we are reading is what ties the two
            # together: a chromosome can never contribute variants without its count, or vice versa.
            m_path = ld_file.parent / (
                ld_file.name[: -len(_LD_SCORE_SUFFIX)] + _M_5_50_SUFFIX
            )
            assert m_path.exists(), (
                f"{ld_file} has no matching {m_path.name}; LD scores must come with the "
                "common-variant count for the same chromosome"
            )
            frames.append(
                scan_dataframe(
                    ld_file,
                    DataFrameReadSpec(DataFrameTextFormat(separator="\t")),
                ).with_columns(
                    narwhals.lit(read_m_5_50(m_path)).alias(LD_SCORE_M_5_50_COL)
                )
            )
        assert frames, (
            f"no *{_LD_SCORE_SUFFIX} files under {asset.path}; nothing to consolidate"
        )
        result = narwhals.concat(frames, how="vertical").sort(
            by=[LD_SCORE_CHROM_COL, LD_SCORE_POS_COL]
        )
        out_path = scratch_dir / "out.parquet"
        result.sink_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(cls, asset_id: str, extracted_ld_score_task: Task):
        source_meta = extracted_ld_score_task.meta
        assert isinstance(source_meta, ReferenceDataDirectoryMeta)
        meta = ReferenceFileMeta(
            group=source_meta.group,
            sub_group=source_meta.sub_group,
            sub_folder=source_meta.sub_folder,
            id=AssetId(asset_id),
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            extension=".parquet",
        )
        return cls(meta=meta, extracted_ld_scores_task=extracted_ld_score_task)
