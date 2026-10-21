"""
Adapter Task: convert ldsc baseline LD scores (chromosome-split, @ pattern) and regression-weight
LD scores into a prepared S-LDSC baseline (see prepared_sldsc_reference). Zero-variance baseline
columns are dropped here, as ldsc's _check_variance does on the full reference, and the baseline
is restricted to SNPs that have a weight LD score (ldsc inner-joins the weights; the HapMap3
weights omit the MHC region).
"""

from pathlib import Path

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.sldsc.ld_score_files import (
    read_ld_score_prefix,
    read_ld_score_sources,
)
from mecfs_bio.build_system.task.sldsc.prepare_ldcts_annotation_set_task import (
    PREPARED_SLDSC_SUB_FOLDER,
)
from mecfs_bio.build_system.task.sldsc.prepared_sldsc_reference import (
    BP_COL,
    CHR_COL,
    SNP_COL,
    W_LD_COL,
    LdScoreTable,
    PreparedBaseline,
    write_prepared_baseline,
)
from mecfs_bio.build_system.wf.base_wf import WF

_ROW_COL = "_row"


@frozen(slots=True)
class PrepareSLDSCBaselineTask(Task):
    """
    baseline_pattern: path inside the baseline directory with @ standing for the chromosome,
        e.g. baseline_v1.2/baseline.@
    weights_pattern: same, for the regression-weight LD scores (exactly one LD column).
    """

    meta: Meta
    baseline_directory_task: Task
    baseline_pattern: str
    weights_directory_task: Task
    weights_pattern: str

    @property
    def deps(self) -> list[Task]:
        return [self.baseline_directory_task, self.weights_directory_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        baseline_dir = _fetch_directory(fetch, self.baseline_directory_task.asset_id)
        weights_dir = _fetch_directory(fetch, self.weights_directory_task.asset_id)
        baseline = _restrict_to_weighted_snps(
            ld=read_ld_score_sources(
                root=baseline_dir,
                prefixes=[self.baseline_pattern],
                drop_zero_variance=True,
            ),
            w_ld=_read_weight_ld(root=weights_dir, pattern=self.weights_pattern),
        )
        out_dir = scratch_dir / "prepared_baseline"
        write_prepared_baseline(out_dir, baseline)
        return DirectoryAsset(out_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        baseline_directory_task: Task,
        baseline_pattern: str,
        weights_directory_task: Task,
        weights_pattern: str,
    ) -> "PrepareSLDSCBaselineTask":
        source_meta = baseline_directory_task.meta
        assert isinstance(source_meta, ReferenceDataDirectoryMeta)
        return cls(
            meta=ReferenceDataDirectoryMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=PREPARED_SLDSC_SUB_FOLDER,
                id=AssetId(asset_id),
            ),
            baseline_directory_task=baseline_directory_task,
            baseline_pattern=baseline_pattern,
            weights_directory_task=weights_directory_task,
            weights_pattern=weights_pattern,
        )


def _fetch_directory(fetch: Fetch, asset_id: AssetId) -> Path:
    asset = fetch(asset_id)
    assert isinstance(asset, DirectoryAsset)
    return asset.path


def _read_weight_ld(root: Path, pattern: str) -> pl.DataFrame:
    frame = read_ld_score_prefix(root, pattern)
    value_cols = [c for c in frame.columns if c not in (SNP_COL, CHR_COL, BP_COL)]
    assert len(value_cols) == 1, (
        f"weight LD scores must have one column, got {value_cols}"
    )
    return frame.select(pl.col(SNP_COL), pl.col(value_cols[0]).alias(W_LD_COL))


def _restrict_to_weighted_snps(
    ld: LdScoreTable, w_ld: pl.DataFrame
) -> PreparedBaseline:
    """Baseline rows that have a weight LD score, in baseline order, with w_ld aligned to them."""
    aligned = (
        ld.snps.select(SNP_COL)
        .with_row_index(_ROW_COL)
        .join(w_ld, on=SNP_COL, how="inner")
        .sort(_ROW_COL)
    )
    return PreparedBaseline(
        ld=ld.take_rows(aligned[_ROW_COL].to_numpy()),
        w_ld=aligned.select(SNP_COL, W_LD_COL),
    )
