"""Explain a polyfun-vs-uniform SUSIE result in annotation terms.

Computes the prior lift m*pi_i and the local annotation contrast
C_c(i) = gamma_raw_c * (a_ic - abar_c), where abar_c is the uniform-run
PIP-weighted mean of annotation c over all locus variants. Aggregates the
contrast to families, selects the top families at the focal (max-PIP-polyfun)
variant, and writes two docs-facing display tables (a top-line table with no
annotation columns and a wide detailed table carrying a per-family contrast
column) plus detail tables.
"""

import json
from pathlib import Path, PurePath
from typing import cast

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
from mecfs_bio.build_system.task.ppp_database.allele_key import unordered_allele_key
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

# Two docs-facing display tables. The top-line table is the headline result
# (identifier + per-run credible set / PIP / prior lift, no annotation columns);
# the detailed table adds one column per annotation family carrying the local
# contrast gamma_raw_c * (a_ic - abar_c) so a reader can see which families drove
# each PIP change. The detailed table is wide and is meant to be scrolled.
TOP_LINE_DISPLAY_TABLE_FILENAME = "top_line_display_table.parquet"
DETAILED_DISPLAY_TABLE_FILENAME = "detailed_display_table.parquet"
PER_ANNOTATION_CONTRAST_FILENAME = "per_annotation_contrast.parquet"
PER_FAMILY_CONTRAST_FILENAME = "per_family_contrast.parquet"
PRIOR_LIFT_FILENAME = "prior_lift.parquet"
SELECTION_JSON_FILENAME = "selection.json"
# Keys inside selection.json. The plot task reads SELECTION_IMPORTANT_FAMILIES_KEY
# back to decide which family panels to draw, so both sides share these constants.
SELECTION_FOCAL_VARIANT_KEY = "focal_variant"
SELECTION_IMPORTANT_FAMILIES_KEY = "important_families"

CALLOUTS_FILENAME = "callouts.parquet"
CALLOUT_CS_COL = "cs"
CALLOUT_PIP_PF_COL = "pip_pf"
CALLOUT_PIP_U_COL = "pip_u"
CALLOUT_LABEL_COL = "label"
# Fixed schema so an empty callout set still round-trips through parquet and the
# plot task can read a well-typed (possibly zero-row) frame.
_CALLOUT_SCHEMA: dict[str, pl.DataType] = {
    GWASLAB_CHROM_COL: pl.Int64(),
    GWASLAB_POS_COL: pl.Int64(),
    GWASLAB_EFFECT_ALLELE_COL: pl.String(),
    GWASLAB_NON_EFFECT_ALLELE_COL: pl.String(),
    CALLOUT_CS_COL: pl.Int32(),
    CALLOUT_PIP_PF_COL: pl.Float64(),
    CALLOUT_PIP_U_COL: pl.Float64(),
    CALLOUT_LABEL_COL: pl.String(),
}
# Selection thresholds (see design doc). Change-based only; no absolute PIP floor.
_DOMINANCE_MARGIN = 0.05
_PRIOR_EFFECT_MARGIN = 0.10
_MAX_CALLOUT_FAMILIES = 3

DISP_CHR = "chr"
DISP_POS = "pos"
DISP_EA = "ea"
DISP_NEA = "nea"
DISP_CS_PF = "cs_pf"
DISP_CS_U = "cs_u"
DISP_PIP_PF = "pip_pf"
DISP_PIP_U = "pip_u"
DISP_LIFT = "lift"
# Prefix on the detailed table's per-family contrast columns, e.g. annot_coding.
DISP_ANNOT_PREFIX = "annot_"

FAMILY_CONTRAST_COL = "family_contrast"
FAMILY_SCALED_COL = "family_scaled"  # sum_c gamma_raw_c * a_ic (NOT the contrast)
CONTRAST_COL = "contrast"
_ANNOT_BP_COL = "BP"
_ANNOT_A1_COL = "A1"
_ANNOT_A2_COL = "A2"
_ALLELE_KEY_COL = "allele_key"

_KEY = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
]
_CS_NUMBER_COL = "cs_number"
# Canonical key dtypes. A SUSIE run that finds no credible sets writes an empty
# combined_cs whose parquet columns default to Float64; without this the empty
# frame's CHR (f64) would poison the union_keys concat and break the i64 join in
# _contrasts. Casting every keyed frame to these types keeps all _KEY joins
# consistent regardless of whether a run found a credible set.
_KEY_SCHEMA: dict[str, pl.DataType] = {
    GWASLAB_CHROM_COL: pl.Int64(),
    GWASLAB_POS_COL: pl.Int64(),
    GWASLAB_EFFECT_ALLELE_COL: pl.String(),
    GWASLAB_NON_EFFECT_ALLELE_COL: pl.String(),
}
# The annotation source carries alleles (A1/A2), so it is joined to a run's
# variants allele-aware on (CHR, POS, unordered-allele-key): each allele of a
# multiallelic site matches its own annotation row. The run side supplies EA/NEA
# (and hence the allele key) to the result.
_ANNOT_KEY = [GWASLAB_CHROM_COL, GWASLAB_POS_COL, _ALLELE_KEY_COL]


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

        focal = pf_variants.sort(PIP_COLUMN, descending=True).head(1)
        focal_key = {k: focal[k][0] for k in _KEY}
        focal_families = _select_families(
            per_family, focal_key, self.n_important_families
        )

        top_line = _top_line_display_table(
            union_keys=union_keys,
            pf_variants=pf_variants,
            uni_variants=uni_variants,
            cs_pf=cs_pf,
            cs_u=cs_u,
        )
        detailed = _detailed_display_table(
            union_keys=union_keys,
            pf_variants=pf_variants,
            uni_variants=uni_variants,
            cs_pf=cs_pf,
            cs_u=cs_u,
            per_family=per_family,
        )

        per_annot.write_parquet(scratch_dir / PER_ANNOTATION_CONTRAST_FILENAME)
        per_family.write_parquet(scratch_dir / PER_FAMILY_CONTRAST_FILENAME)
        pf_variants.select(*_KEY, DISP_LIFT).write_parquet(
            scratch_dir / PRIOR_LIFT_FILENAME
        )
        top_line.write_parquet(scratch_dir / TOP_LINE_DISPLAY_TABLE_FILENAME)
        detailed.write_parquet(scratch_dir / DETAILED_DISPLAY_TABLE_FILENAME)
        (scratch_dir / SELECTION_JSON_FILENAME).write_text(
            json.dumps(
                {
                    SELECTION_FOCAL_VARIANT_KEY: {
                        "chr": int(focal_key[GWASLAB_CHROM_COL]),
                        "pos": int(focal_key[GWASLAB_POS_COL]),
                        "ea": focal_key[GWASLAB_EFFECT_ALLELE_COL],
                        "nea": focal_key[GWASLAB_NON_EFFECT_ALLELE_COL],
                    },
                    SELECTION_IMPORTANT_FAMILIES_KEY: focal_families,
                },
                indent=2,
                sort_keys=True,
            )
        )

        family_sd = _family_background_sd(uni_annot, annot_cols, gamma, family)
        callouts = _build_callouts(
            pf_variants=pf_variants,
            uni_variants=uni_variants,
            cs_pf=cs_pf,
            per_family=per_family,
            family_sd=family_sd,
        )
        callouts.write_parquet(scratch_dir / CALLOUTS_FILENAME)
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
    """filtered_gwas keyed rows + the run's PIP, in the same order, with the
    unordered allele key used to join the annotation matrix allele-aware."""
    gwas = (
        pl.read_parquet(run_dir / FILTERED_GWAS_FILENAME)
        .select(_KEY)
        .with_columns(pl.col(k).cast(dt) for k, dt in _KEY_SCHEMA.items())
    )
    pip = pl.read_parquet(run_dir / PIP_FILENAME).select(PIP_COLUMN)
    return gwas.hstack(pip).with_columns(
        unordered_allele_key(
            GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
        ).alias(_ALLELE_KEY_COL)
    )


def _load_cs_numbers(run_dir: Path) -> pl.DataFrame:
    """One row per credible-set variant with its 1-based L-index (lowest if many)."""
    cs = pl.read_parquet(run_dir / COMBINED_CS_FILENAME)
    if cs.height == 0:
        return pl.DataFrame(schema={**_KEY_SCHEMA, _CS_NUMBER_COL: pl.Int32()})
    return (
        cs.with_columns(
            *(pl.col(k).cast(dt) for k, dt in _KEY_SCHEMA.items()),
            pl.col(CS_COLUMN).str.replace("L", "").cast(pl.Int32).alias(_CS_NUMBER_COL),
        )
        .group_by(_KEY)
        .agg(pl.col(_CS_NUMBER_COL).min())
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
        .select("CHR", _ANNOT_BP_COL, _ANNOT_A1_COL, _ANNOT_A2_COL, *annot_cols)
        .collect()
        .to_polars()
    )
    # Rename BP->POS and derive the unordered allele key so each annotation row
    # is matched to the run variant with the same alleles (regardless of which
    # allele each side labels "effect"). A1/A2 are dropped once the key is built.
    result = (
        frame.rename({_ANNOT_BP_COL: GWASLAB_POS_COL})
        .with_columns(
            unordered_allele_key(_ANNOT_A1_COL, _ANNOT_A2_COL).alias(_ALLELE_KEY_COL)
        )
        .drop(_ANNOT_A1_COL, _ANNOT_A2_COL)
    )
    _assert_annotation_keys_unique(result, chrom, bp_min, bp_max)
    return result


def _assert_annotation_keys_unique(
    annot: pl.DataFrame, chrom: int, bp_min: int, bp_max: int
) -> None:
    """Fail fast if the annotation slice has more than one row per
    (CHR, POS, allele-key) within this locus window.

    The annotation matrix is built unique on (CHR, BP, unordered-allele-key) (see
    BuildBaselineLFAnnotationParquetTask), so this holds by construction; an
    undetected duplicate would silently cross-multiply a run's variant rows into
    doubled or misattributed contrast/family_scaled values. Asserting on the
    annotation slice itself, rather than on a join result, localizes the cause.
    """
    n_rows = annot.height
    n_unique = annot.select(_ANNOT_KEY).n_unique()
    if n_unique != n_rows:
        raise ValueError(
            f"Annotation source has {n_rows - n_unique} duplicate "
            f"(CHR, POS, allele-key) row(s) within locus "
            f"chr{chrom}:{bp_min}-{bp_max}; the annotation join keys on that "
            "tuple, so duplicates would silently cross-multiply variant rows; "
            "refusing to proceed."
        )


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


def _format_callout_label(
    focal_key: dict, families: list[tuple[AnnotationFamily, str]]
) -> str:
    """Render one callout string: pos:nea:ea, then the key families with their
    strength markers. No families -> just the variant id (no parentheses)."""
    head = (
        f"{focal_key[GWASLAB_POS_COL]}:"
        f"{focal_key[GWASLAB_NON_EFFECT_ALLELE_COL]}:"
        f"{focal_key[GWASLAB_EFFECT_ALLELE_COL]}"
    )
    if not families:
        return head
    # Full family names (underscores -> spaces) so the callout reads at a glance,
    # e.g. "conserved ++, coding +" rather than the compact display-table forms.
    inner = ", ".join(f"{fam.replace('_', ' ')} {marker}" for fam, marker in families)
    return f"{head} ({inner})"


def _callout_families(
    per_family: pl.DataFrame,
    focal_key: dict,
    family_sd: dict[str, float],
    max_families: int,
) -> list[tuple[AnnotationFamily, str]]:
    """The key families for one flagged variant: those whose per-family contrast
    is positive (elevated at this variant) and exceeds one background SD, bucketed
    1-2 SD -> '+', >2 SD -> '++'. Top max_families by z, z descending. Families
    with a degenerate (<= 0) background SD are skipped."""
    focal = per_family
    for k, v in focal_key.items():
        focal = focal.filter(pl.col(k) == v)
    scored: list[tuple[float, AnnotationFamily, str]] = []
    for row in focal.iter_rows(named=True):
        fam = cast(AnnotationFamily, row[FAMILY_COL])
        diff = row[FAMILY_CONTRAST_COL]
        sd = family_sd.get(fam, 0.0)
        if diff <= 0.0 or sd <= 0.0:
            continue
        z = diff / sd
        if z <= 1.0:
            continue
        scored.append((z, fam, "++" if z > 2.0 else "+"))
    scored.sort(key=lambda t: t[0], reverse=True)
    return [(fam, marker) for _, fam, marker in scored[:max_families]]


def _family_background_sd(
    uni_annot: pl.DataFrame,
    annot_cols: list[str],
    gamma: dict[str, float],
    family: dict[str, str],
) -> dict[str, float]:
    """Per family, the uniform-PIP-weighted standard deviation of family_scaled
    over the uniform-run variants. Falls back to equal weights when the uniform
    run carried no signal (total PIP <= 0), matching the abar fallback."""
    fs = _family_scaled(uni_annot, annot_cols, gamma, family)
    weight = uni_annot.select(*_KEY, pl.col(PIP_COLUMN).alias("w"))
    if uni_annot[PIP_COLUMN].sum() <= 0.0:
        weight = weight.with_columns(pl.lit(1.0).alias("w"))
    stats = (
        fs.join(weight, on=_KEY, how="inner")
        .group_by(FAMILY_COL)
        .agg(
            (pl.col("w") * pl.col(FAMILY_SCALED_COL)).sum().alias("wx"),
            (pl.col("w") * pl.col(FAMILY_SCALED_COL) ** 2).sum().alias("wx2"),
            pl.col("w").sum().alias("wsum"),
        )
    )
    out: dict[str, float] = {}
    for row in stats.iter_rows(named=True):
        wsum = row["wsum"]
        if wsum <= 0.0:
            continue
        mean = row["wx"] / wsum
        var = max(row["wx2"] / wsum - mean * mean, 0.0)
        out[row[FAMILY_COL]] = float(np.sqrt(var))
    return out


def _build_callouts(
    pf_variants: pl.DataFrame,
    uni_variants: pl.DataFrame,
    cs_pf: pl.DataFrame,
    per_family: pl.DataFrame,
    family_sd: dict[str, float],
) -> pl.DataFrame:
    """One callout row per polyfun credible set whose top-PIP variant clears both
    gates: PIP >= _DOMINANCE_MARGIN above the next-highest PIP in the same CS, and
    PIP >= _PRIOR_EFFECT_MARGIN above the same variant's uniform PIP (0 if absent
    from the uniform run)."""
    cs = cs_pf.join(pf_variants.select(*_KEY, PIP_COLUMN), on=_KEY, how="inner")
    uni_pip = {
        tuple(row[k] for k in _KEY): row[PIP_COLUMN]
        for row in uni_variants.select(*_KEY, PIP_COLUMN).iter_rows(named=True)
    }
    rows: list[dict] = []
    for (cs_number,), grp in cs.group_by(_CS_NUMBER_COL, maintain_order=True):
        grp = grp.sort(PIP_COLUMN, descending=True)
        top = grp.row(0, named=True)
        top_pip = top[PIP_COLUMN]
        next_pip = grp[PIP_COLUMN][1] if grp.height > 1 else 0.0
        if top_pip - next_pip < _DOMINANCE_MARGIN:
            continue
        focal_key = {k: top[k] for k in _KEY}
        u = uni_pip.get(tuple(top[k] for k in _KEY), 0.0)
        if top_pip - u < _PRIOR_EFFECT_MARGIN:
            continue
        families = _callout_families(
            per_family, focal_key, family_sd, _MAX_CALLOUT_FAMILIES
        )
        rows.append(
            {
                **{k: focal_key[k] for k in _KEY},
                CALLOUT_CS_COL: int(cs_number),
                CALLOUT_PIP_PF_COL: float(top_pip),
                CALLOUT_PIP_U_COL: float(u),
                CALLOUT_LABEL_COL: _format_callout_label(focal_key, families),
            }
        )
    return pl.DataFrame(rows, schema=_CALLOUT_SCHEMA)


_DISPLAY_BASE_COLS = [
    DISP_CHR,
    DISP_POS,
    DISP_EA,
    DISP_NEA,
    DISP_CS_PF,
    DISP_CS_U,
    DISP_PIP_PF,
    DISP_PIP_U,
    DISP_LIFT,
]


def _display_base(
    union_keys: pl.DataFrame,
    pf_variants: pl.DataFrame,
    uni_variants: pl.DataFrame,
    cs_pf: pl.DataFrame,
    cs_u: pl.DataFrame,
) -> pl.DataFrame:
    """The identifier + per-run credible-set / PIP / prior-lift columns shared by
    both display tables, keyed (still under _KEY names) on the union of the two
    runs' credible-set variants so family columns can be joined on before the
    final rename."""
    return (
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
            cs_pf.select(*_KEY, pl.col(_CS_NUMBER_COL).alias(DISP_CS_PF)),
            on=_KEY,
            how="left",
        )
        .join(
            cs_u.select(*_KEY, pl.col(_CS_NUMBER_COL).alias(DISP_CS_U)),
            on=_KEY,
            how="left",
        )
    )


def _finalize_display(out: pl.DataFrame, extra_cols: list[str]) -> pl.DataFrame:
    """Rename the key columns to their display names, order base-then-extra, and
    sort by descending polyfun PIP."""
    out = out.rename(
        {
            GWASLAB_CHROM_COL: DISP_CHR,
            GWASLAB_POS_COL: DISP_POS,
            GWASLAB_EFFECT_ALLELE_COL: DISP_EA,
            GWASLAB_NON_EFFECT_ALLELE_COL: DISP_NEA,
        }
    ).with_columns(pl.col(DISP_CHR).cast(pl.Int32), pl.col(DISP_POS).cast(pl.Int32))
    return out.select(_DISPLAY_BASE_COLS + extra_cols).sort(
        DISP_PIP_PF, descending=True, nulls_last=True
    )


def _top_line_display_table(
    union_keys: pl.DataFrame,
    pf_variants: pl.DataFrame,
    uni_variants: pl.DataFrame,
    cs_pf: pl.DataFrame,
    cs_u: pl.DataFrame,
) -> pl.DataFrame:
    """The headline result table: identifiers, per-run credible-set numbers, both
    PIPs, and the prior lift. No annotation-family columns."""
    base = _display_base(union_keys, pf_variants, uni_variants, cs_pf, cs_u)
    return _finalize_display(base, [])


def _families_in_canonical_order() -> list[AnnotationFamily]:
    """All eleven families in the fixed taxonomy order (the key order of
    FAMILY_SHORT_LABELS), so the detailed table always has the same column set. A
    family with no annotations at this locus is emitted as an all-null column."""
    return list(FAMILY_SHORT_LABELS)


def _detailed_display_table(
    union_keys: pl.DataFrame,
    pf_variants: pl.DataFrame,
    uni_variants: pl.DataFrame,
    cs_pf: pl.DataFrame,
    cs_u: pl.DataFrame,
    per_family: pl.DataFrame,
) -> pl.DataFrame:
    """The wide table: the top-line columns plus one column per annotation family
    carrying that family's local contrast gamma_raw_c * (a_ic - abar_c) (summed
    over the family's annotations), so a reader can see which families pushed each
    variant's PIP up or down. Family columns use the full family name with an
    annot_ prefix (e.g. annot_coding)."""
    out = _display_base(union_keys, pf_variants, uni_variants, cs_pf, cs_u)
    families = _families_in_canonical_order()
    for fam in families:
        col = f"{DISP_ANNOT_PREFIX}{fam}"
        fam_col = per_family.filter(pl.col(FAMILY_COL) == fam).select(
            *_KEY, pl.col(FAMILY_CONTRAST_COL).alias(col)
        )
        out = out.join(fam_col, on=_KEY, how="left")
    return _finalize_display(out, [f"{DISP_ANNOT_PREFIX}{fam}" for fam in families])
