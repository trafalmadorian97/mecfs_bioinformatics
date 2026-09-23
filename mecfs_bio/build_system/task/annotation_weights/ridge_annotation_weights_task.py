"""Fit a ridge surrogate of the polyfun prior on the baseline-LF annotations.

Regresses snpvar_bin (the polyfun binned per-SNP heritability prior actually used
in fine mapping, not the original S-LDSC tau_c) on the 187 annotations, genome
wide. The fit is done from per-chromosome cross-product sufficient statistics, so
the full design matrix is never held in memory (see chromosome_blocked_ridge);
alpha is chosen by leave-one-chromosome-out. Outputs raw-scale coefficients gamma_raw (used by the
explainability contrast) and standardized coefficients gamma_standardized (for
global importance ranking), plus each annotation's family.


NOTE: the annotations are standardized before fitting the ridge regression model.  The resulting annotation ridge regression
coefficients are then un-standardized.  This approach was chosen so that the coefficients of all annotations are shrunk
by the ridge penalty on the same standardized scale.
"""

import json
from pathlib import Path

import narwhals
import numpy as np
import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.annotation_weights.build_baseline_lf_annotation_parquet_task import (
    ANNOT_KEY_COLUMNS,
)
from mecfs_bio.build_system.task.annotation_weights.chromosome_blocked_ridge import (
    ChromRidgeBlock,
    accumulate_block,
    combine,
    fit,
    select_alpha_loco,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import (
    ParquetOutFormat,
    write_df_according_to_format,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.polyfun_annotation_families import family_for_annotation
from mecfs_bio.constants.polyfun_constants import (
    POLYFUN_A1_COL,
    POLYFUN_A2_COL,
    POLYFUN_BP_COL,
    POLYFUN_CHR_COL,
)

WEIGHTS_PARQUET_FILENAME = "weights.parquet"
DIAGNOSTICS_JSON_FILENAME = "diagnostics.json"
ANNOTATION_COL = "annotation"
GAMMA_RAW_COL = "gamma_raw"
GAMMA_STANDARDIZED_COL = "gamma_standardized"
FAMILY_COL = "family"
SNPVAR_COL = "snpvar_bin"
# Both the annotation matrix and snpvar_meta carry alleles (A1/A2), so the
# annotation<->snpvar join is exact on (CHR, BP, A1, A2), both sides
# reference-oriented (A1 == REF). This pairs each allele of a multiallelic site --
# and each orientation of a mirrored indel -- with its own snpvar_bin.
_JOIN_KEYS = [POLYFUN_CHR_COL, POLYFUN_BP_COL, POLYFUN_A1_COL, POLYFUN_A2_COL]

_DEFAULT_ALPHAS: tuple[float, ...] = (0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0)


@frozen
class RidgeAnnotationWeightsTask(Task):
    meta: Meta
    annotation_parquet_task: Task
    snpvar_meta_task: Task
    alphas: tuple[float, ...] = _DEFAULT_ALPHAS

    @property
    def deps(self) -> list["Task"]:
        return [self.annotation_parquet_task, self.snpvar_meta_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        annot_asset = fetch(self.annotation_parquet_task.asset_id)
        assert isinstance(annot_asset, FileAsset)
        meta_asset = fetch(self.snpvar_meta_task.asset_id)
        meta = (
            scan_dataframe_asset(meta_asset, self.snpvar_meta_task.meta)
            .select(
                POLYFUN_CHR_COL,
                POLYFUN_BP_COL,
                POLYFUN_A1_COL,
                POLYFUN_A2_COL,
                SNPVAR_COL,
            )
            .collect()
            .to_polars()
            .unique(subset=_JOIN_KEYS)
        )

        annot_columns = _annotation_columns(annot_asset.path)
        per_chrom = _accumulate_per_chromosome(annot_asset.path, annot_columns, meta)
        selection = select_alpha_loco(per_chrom, self.alphas)
        ridge = fit(combine(list(per_chrom.values())), alpha=selection.alpha)

        weights = pl.DataFrame(
            {
                ANNOTATION_COL: annot_columns,
                GAMMA_RAW_COL: ridge.beta_raw.tolist(),
                GAMMA_STANDARDIZED_COL: ridge.beta_std.tolist(),
                # Any parquet column not classifiable by family_for_annotation
                # hard-fails here by design: this is the de-facto guard that the
                # built parquet's columns match the known baseline-LF annotation set.
                FAMILY_COL: [family_for_annotation(c) for c in annot_columns],
            }
        )
        write_df_according_to_format(
            df=narwhals.from_native(weights).lazy(),
            out_path=scratch_dir / WEIGHTS_PARQUET_FILENAME,
            out_format=ParquetOutFormat(),
        )
        diagnostics = {
            "alpha": selection.alpha,
            "intercept": ridge.intercept,
            "heldout_r2_per_chrom": {
                str(c): r2 for c, r2 in selection.r2_per_chrom.items()
            },
            "mean_heldout_r2": selection.mean_r2,
            # Unweighted blocks: each block's total weight is its variant count.
            "n_variants": round(sum(b.sw for b in per_chrom.values())),
        }
        (scratch_dir / DIAGNOSTICS_JSON_FILENAME).write_text(
            json.dumps(diagnostics, indent=2, sort_keys=True)
        )
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        annotation_parquet_task: Task,
        snpvar_meta_task: Task,
        alphas: tuple[float, ...] = _DEFAULT_ALPHAS,
    ) -> "RidgeAnnotationWeightsTask":
        # Derive the output directory meta from the primary dependency (the
        # annotation parquet), reusing group/sub_group/sub_folder - the
        # CompressedCSVToParquetTask.create pattern applied to a DirMeta.
        source_meta = annotation_parquet_task.meta
        if not isinstance(source_meta, ReferenceFileMeta):
            raise ValueError(f"Unknown meta for annotation parquet task: {source_meta}")
        meta = ReferenceDataDirectoryMeta(
            group=source_meta.group,
            sub_group=source_meta.sub_group,
            sub_folder=source_meta.sub_folder,
            id=AssetId(asset_id),
        )
        return cls(
            meta=meta,
            annotation_parquet_task=annotation_parquet_task,
            snpvar_meta_task=snpvar_meta_task,
            alphas=alphas,
        )


def _annotation_columns(annot_path: Path) -> list[str]:
    schema = pl.scan_parquet(annot_path).collect_schema()
    return [c for c in schema.names() if c not in ANNOT_KEY_COLUMNS]


def _accumulate_per_chromosome(
    annot_path: Path, annot_columns: list[str], meta: pl.DataFrame
) -> dict[int, ChromRidgeBlock]:
    chroms = (
        pl.scan_parquet(annot_path)
        .select(POLYFUN_CHR_COL)
        .unique()
        .collect()
        .to_series()
        .to_list()
    )
    per_chrom: dict[int, ChromRidgeBlock] = {}
    for chrom in sorted(chroms):
        annot_chrom = (
            pl.scan_parquet(annot_path)
            .filter(pl.col(POLYFUN_CHR_COL) == chrom)
            .collect()
        )
        frame = annot_chrom.join(meta, on=_JOIN_KEYS, how="inner")
        # meta is unique on the join key and the annotation matrix is unique on
        # (CHR, BP, A1, A2), so the inner join must not multiply rows.
        assert frame.height <= annot_chrom.height, (
            f"annotation<->snpvar join multiplied rows on chr{chrom}: "
            f"{annot_chrom.height} -> {frame.height}"
        )
        x = frame.select(annot_columns).to_numpy().astype(np.float64)
        y = frame.select(SNPVAR_COL).to_numpy().ravel().astype(np.float64)
        per_chrom[chrom] = accumulate_block(x=x, y=y)
    return per_chrom
