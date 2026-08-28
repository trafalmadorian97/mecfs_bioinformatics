"""Explain a polyfun-vs-uniform SUSIE result in annotation terms.

Computes the prior lift m*pi_i and the local annotation contrast
C_c(i) = gamma_raw_c * (a_ic - abar_c), where abar_c is the uniform-run
PIP-weighted mean of annotation c over all locus variants. Aggregates the
contrast to families, selects the top families at the focal (max-PIP-polyfun)
variant, and writes the docs-facing display table plus detail tables.
"""

import json
from pathlib import Path, PurePath

import narwhals as nw
import numpy as np
import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    ANNOTATION_COL,
    FAMILY_COL,
    GAMMA_RAW_COL,
    WEIGHTS_PARQUET_FILENAME,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    COMBINED_CS_FILENAME,
    CS_COLUMN,
    FILTERED_GWAS_FILENAME,
    PIP_COLUMN,
    PIP_FILENAME,
    PRIOR_FILENAME,
    PRIOR_WEIGHT_COLUMN,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)
from mecfs_bio.constants.polyfun_annotation_families import (
    FAMILY_SHORT_LABELS,
    AnnotationFamily,
)

DISPLAY_TABLE_FILENAME = "display_table.parquet"
PER_ANNOTATION_CONTRAST_FILENAME = "per_annotation_contrast.parquet"
PER_FAMILY_CONTRAST_FILENAME = "per_family_contrast.parquet"
PRIOR_LIFT_FILENAME = "prior_lift.parquet"
SELECTION_JSON_FILENAME = "selection.json"

DISP_CHR = "chr"
DISP_POS = "pos"
DISP_EA = "ea"
DISP_NEA = "nea"
DISP_CS_PF = "cs_pf"
DISP_CS_U = "cs_u"
DISP_PIP_PF = "pip_pf"
DISP_PIP_U = "pip_u"
DISP_LIFT = "lift"

FAMILY_CONTRAST_COL = "family_contrast"
FAMILY_SCALED_COL = "family_scaled"  # sum_c gamma_raw_c * a_ic (NOT the contrast)
CONTRAST_COL = "contrast"
_ANNOT_BP_COL = "BP"

_KEY = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
]
# The annotation source carries no alleles (see _load_annotations), so it can
# only be joined to a run's variants on (CHR, POS); the run side of the join
# is what supplies EA/NEA to the result.
_ANNOT_KEY = [GWASLAB_CHROM_COL, GWASLAB_POS_COL]


@frozen
class PolyfunExplainContrastTask(Task):
    meta: Meta
    susie_uniform_task: Task
    susie_polyfun_task: Task
    ridge_weights_task: Task
    annotation_parquet_task: Task
    n_important_families: int = 3

    @property
    def deps(self) -> list["Task"]:
        return [
            self.susie_uniform_task,
            self.susie_polyfun_task,
            self.ridge_weights_task,
            self.annotation_parquet_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        uni_dir = _dir(fetch, self.susie_uniform_task)
        pf_dir = _dir(fetch, self.susie_polyfun_task)

        uni_variants = _load_run_variants(uni_dir)
        pf_variants = _load_run_variants(pf_dir)
        pf_prior = pl.read_parquet(pf_dir / PRIOR_FILENAME)

        weights = _load_weights(fetch, self.ridge_weights_task)
        annot_cols = weights[ANNOTATION_COL].to_list()
        gamma = dict(zip(weights[ANNOTATION_COL], weights[GAMMA_RAW_COL]))
        family = dict(zip(weights[ANNOTATION_COL], weights[FAMILY_COL]))

        chrom = int(pf_variants[GWASLAB_CHROM_COL][0])
        bp_min = int(pf_variants[GWASLAB_POS_COL].to_numpy().min())
        bp_max = int(pf_variants[GWASLAB_POS_COL].to_numpy().max())
        annot = _load_annotations(
            fetch, self.annotation_parquet_task, chrom, bp_min, bp_max, annot_cols
        )

        # abar_c: uniform PIP-weighted mean of each annotation over all uniform vars.
        # If the uniform run found no signal (all PIPs ~0), fall back to an
        # unweighted mean so every locus variant contributes equally.
        uni_annot = uni_variants.join(annot, on=_ANNOT_KEY, how="inner")
        w = uni_annot[PIP_COLUMN].to_numpy()
        if w.sum() <= 0.0:
            w = None
        abar = {
            c: float(np.average(uni_annot[c].to_numpy(), weights=w)) for c in annot_cols
        }

        # prior lift on the polyfun run.
        prior_w = pf_prior[PRIOR_WEIGHT_COLUMN].to_numpy()
        m = pf_variants.height
        lift = m * prior_w / prior_w.sum()
        pf_variants = pf_variants.with_columns(pl.Series(name=DISP_LIFT, values=lift))

        # attribution row set: union of the two runs' credible-set variants.
        cs_pf = _load_cs_numbers(pf_dir)
        cs_u = _load_cs_numbers(uni_dir)
        union_keys = pl.concat(
            [cs_pf.select(_KEY), cs_u.select(_KEY)], how="vertical"
        ).unique()

        pf_annot = pf_variants.join(annot, on=_ANNOT_KEY, how="inner")

        per_annot, per_family = _contrasts(
            pf_annot, union_keys, annot_cols, gamma, family, abar
        )
        family_scaled = _family_scaled(pf_annot, annot_cols, gamma, family)

        focal = pf_variants.sort(PIP_COLUMN, descending=True).head(1)
        focal_key = {k: focal[k][0] for k in _KEY}
        focal_families = _select_families(
            per_family, focal_key, self.n_important_families
        )

        display = _display_table(
            union_keys=union_keys,
            pf_variants=pf_variants,
            uni_variants=uni_variants,
            cs_pf=cs_pf,
            cs_u=cs_u,
            family_scaled=family_scaled,
            focal_families=focal_families,
        )

        per_annot.write_parquet(scratch_dir / PER_ANNOTATION_CONTRAST_FILENAME)
        per_family.write_parquet(scratch_dir / PER_FAMILY_CONTRAST_FILENAME)
        pf_variants.select(*_KEY, DISP_LIFT).write_parquet(
            scratch_dir / PRIOR_LIFT_FILENAME
        )
        display.write_parquet(scratch_dir / DISPLAY_TABLE_FILENAME)
        (scratch_dir / SELECTION_JSON_FILENAME).write_text(
            json.dumps(
                {
                    "focal_variant": {
                        "chr": int(focal_key[GWASLAB_CHROM_COL]),
                        "pos": int(focal_key[GWASLAB_POS_COL]),
                        "ea": focal_key[GWASLAB_EFFECT_ALLELE_COL],
                        "nea": focal_key[GWASLAB_NON_EFFECT_ALLELE_COL],
                    },
                    "important_families": focal_families,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        susie_uniform_task: Task,
        susie_polyfun_task: Task,
        ridge_weights_task: Task,
        annotation_parquet_task: Task,
        n_important_families: int = 3,
    ) -> "PolyfunExplainContrastTask":
        source_meta = susie_polyfun_task.meta
        if not isinstance(source_meta, ResultDirectoryMeta):
            raise ValueError(f"Unknown meta for polyfun susie task: {source_meta}")
        meta = ResultDirectoryMeta(
            id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=PurePath("analysis"),
        )
        return cls(
            meta=meta,
            susie_uniform_task=susie_uniform_task,
            susie_polyfun_task=susie_polyfun_task,
            ridge_weights_task=ridge_weights_task,
            annotation_parquet_task=annotation_parquet_task,
            n_important_families=n_important_families,
        )


def _dir(fetch: Fetch, task: Task) -> Path:
    asset = fetch(task.asset_id)
    assert isinstance(asset, DirectoryAsset)
    return asset.path


def _load_run_variants(run_dir: Path) -> pl.DataFrame:
    """filtered_gwas keyed rows + the run's PIP, in the same order."""
    gwas = pl.read_parquet(run_dir / FILTERED_GWAS_FILENAME).select(_KEY)
    pip = pl.read_parquet(run_dir / PIP_FILENAME).select(PIP_COLUMN)
    return gwas.hstack(pip)


def _load_cs_numbers(run_dir: Path) -> pl.DataFrame:
    """One row per credible-set variant with its 1-based L-index (lowest if many)."""
    cs = pl.read_parquet(run_dir / COMBINED_CS_FILENAME)
    if cs.height == 0:
        return pl.DataFrame(schema={**{k: cs.schema.get(k, pl.Int64) for k in _KEY}})
    return (
        cs.with_columns(
            pl.col(CS_COLUMN).str.replace("L", "").cast(pl.Int32).alias("cs_number")
        )
        .group_by(_KEY)
        .agg(pl.col("cs_number").min())
    )


def _load_weights(fetch: Fetch, task: Task) -> pl.DataFrame:
    asset = fetch(task.asset_id)
    assert isinstance(asset, (FileAsset, DirectoryAsset))
    path = (
        asset.path
        if isinstance(asset, FileAsset)
        else asset.path / WEIGHTS_PARQUET_FILENAME
    )
    return pl.read_parquet(path)


def _load_annotations(
    fetch: Fetch,
    task: Task,
    chrom: int,
    bp_min: int,
    bp_max: int,
    annot_cols: list[str],
) -> pl.DataFrame:
    """Locus-windowed annotation slice, keyed like the run variants (BP -> POS)."""
    asset = fetch(task.asset_id)
    frame = (
        scan_dataframe_asset(asset, task.meta)
        .filter(
            (nw.col("CHR") == chrom)
            & (nw.col(_ANNOT_BP_COL) >= bp_min)
            & (nw.col(_ANNOT_BP_COL) <= bp_max)
        )
        .select("CHR", _ANNOT_BP_COL, *annot_cols)
        .collect()
        .to_polars()
    )
    # annotation has no alleles; broadcast to the run's alleles via the join key by
    # renaming BP->POS and letting the (CHR,POS) match carry EA/NEA from the run.
    return frame.rename({_ANNOT_BP_COL: GWASLAB_POS_COL})


def _contrasts(
    pf_annot: pl.DataFrame,
    union_keys: pl.DataFrame,
    annot_cols: list[str],
    gamma: dict[str, float],
    family: dict[str, str],
    abar: dict[str, float],
) -> tuple[pl.DataFrame, pl.DataFrame]:
    rows = pf_annot.join(union_keys, on=_KEY, how="inner")
    long = rows.unpivot(
        on=annot_cols, index=_KEY, variable_name=ANNOTATION_COL, value_name="a_ic"
    ).with_columns(
        (
            pl.col(ANNOTATION_COL).replace_strict(gamma)
            * (pl.col("a_ic") - pl.col(ANNOTATION_COL).replace_strict(abar))
        ).alias(CONTRAST_COL),
        pl.col(ANNOTATION_COL).replace_strict(family).alias(FAMILY_COL),
    )
    per_annotation = long.select(*_KEY, ANNOTATION_COL, FAMILY_COL, CONTRAST_COL)
    per_family = (
        long.group_by([*_KEY, FAMILY_COL])
        .agg(pl.col(CONTRAST_COL).sum().alias(FAMILY_CONTRAST_COL))
        .sort([*_KEY, FAMILY_COL])
    )
    return per_annotation, per_family


def _family_scaled(
    pf_annot: pl.DataFrame,
    annot_cols: list[str],
    gamma: dict[str, float],
    family: dict[str, str],
) -> pl.DataFrame:
    """Per variant per family: sum_c gamma_raw_c * a_ic (raw scaled value)."""
    long = pf_annot.unpivot(
        on=annot_cols, index=_KEY, variable_name=ANNOTATION_COL, value_name="a_ic"
    ).with_columns(
        (pl.col(ANNOTATION_COL).replace_strict(gamma) * pl.col("a_ic")).alias("scaled"),
        pl.col(ANNOTATION_COL).replace_strict(family).alias(FAMILY_COL),
    )
    return (
        long.group_by([*_KEY, FAMILY_COL])
        .agg(pl.col("scaled").sum().alias(FAMILY_SCALED_COL))
        .sort([*_KEY, FAMILY_COL])
    )


def _select_families(
    per_family: pl.DataFrame, focal_key: dict, n: int
) -> list[AnnotationFamily]:
    focal = per_family
    for k, v in focal_key.items():
        focal = focal.filter(pl.col(k) == v)
    return (
        focal.sort(FAMILY_CONTRAST_COL, descending=True).head(n)[FAMILY_COL].to_list()
    )


def _display_table(
    union_keys: pl.DataFrame,
    pf_variants: pl.DataFrame,
    uni_variants: pl.DataFrame,
    cs_pf: pl.DataFrame,
    cs_u: pl.DataFrame,
    family_scaled: pl.DataFrame,
    focal_families: list[AnnotationFamily],
) -> pl.DataFrame:
    out = (
        union_keys.join(
            pf_variants.select(*_KEY, pl.col(PIP_COLUMN).alias(DISP_PIP_PF), DISP_LIFT),
            on=_KEY,
            how="left",
        )
        .join(
            uni_variants.select(*_KEY, pl.col(PIP_COLUMN).alias(DISP_PIP_U)),
            on=_KEY,
            how="left",
        )
        .join(
            cs_pf.select(*_KEY, pl.col("cs_number").alias(DISP_CS_PF)),
            on=_KEY,
            how="left",
        )
        .join(
            cs_u.select(*_KEY, pl.col("cs_number").alias(DISP_CS_U)),
            on=_KEY,
            how="left",
        )
    )
    for fam in focal_families:
        col = FAMILY_SHORT_LABELS[fam]
        fam_col = family_scaled.filter(pl.col(FAMILY_COL) == fam).select(
            *_KEY, pl.col(FAMILY_SCALED_COL).alias(col)
        )
        out = out.join(fam_col, on=_KEY, how="left")
    out = out.rename(
        {
            GWASLAB_CHROM_COL: DISP_CHR,
            GWASLAB_POS_COL: DISP_POS,
            GWASLAB_EFFECT_ALLELE_COL: DISP_EA,
            GWASLAB_NON_EFFECT_ALLELE_COL: DISP_NEA,
        }
    ).with_columns(pl.col(DISP_CHR).cast(pl.Int32), pl.col(DISP_POS).cast(pl.Int32))
    ordered = [
        DISP_CHR,
        DISP_POS,
        DISP_EA,
        DISP_NEA,
        DISP_CS_PF,
        DISP_CS_U,
        DISP_PIP_PF,
        DISP_PIP_U,
        DISP_LIFT,
    ] + [FAMILY_SHORT_LABELS[f] for f in focal_families]
    return out.select(ordered).sort(DISP_PIP_PF, descending=True, nulls_last=True)
