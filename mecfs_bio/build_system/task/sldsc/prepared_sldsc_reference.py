"""
Format-neutral on-disk contract for S-LDSC reference data.

A prepared annotation set is a directory holding:
    ld_scores.parquet: SNP, CHR, BP, then one float32 column per unique LD-score source,
        rows sorted by CHR then BP with unique SNPs.
    m_5_50.parquet: one row per LD column: column id, genome-wide M_5_50, source.
    models.parquet: one row per analysis label: label, target column id, control column ids.

A prepared baseline holds ld_scores.parquet and m_5_50.parquet for the baseline annotations
plus w_ld.parquet (SNP and the regression-weight LD score). Its rows are restricted to SNPs that
have a weight LD score (the HapMap3 weights exclude the MHC region), so w_ld and the baseline LD
scores share one SNP column.

M_5_50: for each annotation, the sum of its values over reference-panel SNPs with minor allele
frequency above 5% (for a binary annotation, the number of common SNPs it contains). In the LDSC
model E[chi2_j] = N_j * sum_c tau_c * l(j, c) + N_j * a + 1, with l(j, c) the LD score of SNP j
with annotation c and tau_c the per-SNP heritability contributed by annotation c, M_5_50
converts per-SNP effects into the heritability attributed to an annotation
(h2_c = tau_c * M_5_50_c).

Adapter Tasks (one per source format) write this contract; the vectorized S-LDSC Task reads it
and knows nothing about the source format.
"""

from pathlib import Path

import numpy as np
import polars as pl
from attrs import frozen

SNP_COL = "SNP"
CHR_COL = "CHR"
BP_COL = "BP"
W_LD_COL = "w_ld"

LD_SCORES_FILENAME = "ld_scores.parquet"
M_5_50_FILENAME = "m_5_50.parquet"
MODELS_FILENAME = "models.parquet"
W_LD_FILENAME = "w_ld.parquet"

_COLUMN_ID_COL = "column_id"
_M_5_50_COL = "m_5_50"
_SOURCE_COL = "source"
_LABEL_COL = "label"
_TARGET_COL = "target_column"
_CONTROLS_COL = "control_columns"
_KEY_COLS = (SNP_COL, CHR_COL, BP_COL)


@frozen(slots=True)
class LdScoreTable:
    """
    LD scores for a fixed SNP universe.

    snps: SNP, CHR, BP in row order (sorted by CHR then BP, SNP unique).
    column_ids: internal id of each LD column (matrix column order).
    sources: human-readable provenance of each column (e.g. file prefix and column name).
    matrix: (n_snps, n_columns) float32 LD scores.
    m_5_50: (n_columns,) genome-wide M_5_50 per column.
    """

    snps: pl.DataFrame
    column_ids: tuple[str, ...]
    sources: tuple[str, ...]
    matrix: np.ndarray
    m_5_50: np.ndarray

    def __attrs_post_init__(self) -> None:
        assert self.snps.columns == list(_KEY_COLS), (
            f"snps must have columns {_KEY_COLS}, got {self.snps.columns}"
        )
        assert self.snps[SNP_COL].is_unique().all(), "SNP ids must be unique"
        n_cols = len(self.column_ids)
        assert len(set(self.column_ids)) == n_cols, "column ids must be unique"
        assert len(self.sources) == n_cols, "one source per column"
        assert self.matrix.ndim == 2, f"matrix must be 2-D, got {self.matrix.shape}"
        assert self.matrix.shape == (self.snps.height, n_cols), (
            f"matrix shape {self.matrix.shape} != ({self.snps.height}, {n_cols})"
        )
        assert self.matrix.dtype == np.float32, (
            f"matrix must be float32, got {self.matrix.dtype}"
        )
        assert self.m_5_50.ndim == 1 and self.m_5_50.shape[0] == n_cols, (
            f"m_5_50 must have shape ({n_cols},), got {self.m_5_50.shape}"
        )
        assert self.m_5_50.dtype.kind == "f", (
            f"m_5_50 must be float, got {self.m_5_50.dtype}"
        )

    def column_index(self, column_id: str) -> int:
        return self.column_ids.index(column_id)

    def take_rows(self, rows: np.ndarray) -> "LdScoreTable":
        """The table restricted to the given row indices, in that order."""
        return LdScoreTable(
            snps=self.snps.select(pl.all().gather(rows)),
            column_ids=self.column_ids,
            sources=self.sources,
            matrix=self.matrix[rows],
            m_5_50=self.m_5_50,
        )


@frozen(slots=True)
class SLDSCAnnotationModel:
    """One S-LDSC analysis: the reported target column plus the control columns it is fit with."""

    label: str
    target_column: str
    control_columns: tuple[str, ...]

    def __attrs_post_init__(self) -> None:
        assert self.target_column not in self.control_columns, (
            f"model {self.label}: target {self.target_column} is among its own controls"
        )
        assert len(set(self.control_columns)) == len(self.control_columns), (
            f"model {self.label}: duplicate control columns"
        )


@frozen(slots=True)
class PreparedAnnotationSet:
    ld: LdScoreTable
    models: tuple[SLDSCAnnotationModel, ...]

    def __attrs_post_init__(self) -> None:
        labels = [m.label for m in self.models]
        assert len(set(labels)) == len(labels), "model labels must be unique"
        known = set(self.ld.column_ids)
        for model in self.models:
            missing = {model.target_column, *model.control_columns} - known
            assert not missing, (
                f"model {model.label} references unknown columns {missing}"
            )


@frozen(slots=True)
class PreparedBaseline:
    """
    Baseline LD scores plus the regression-weight LD score.

    w_ld: columns SNP and w_ld, row-aligned with ld.snps (same SNPs in the same order).
    """

    ld: LdScoreTable
    w_ld: pl.DataFrame

    def __attrs_post_init__(self) -> None:
        assert self.w_ld.columns == [SNP_COL, W_LD_COL], (
            f"w_ld must have columns {[SNP_COL, W_LD_COL]}, got {self.w_ld.columns}"
        )
        assert self.w_ld[SNP_COL].equals(self.ld.snps[SNP_COL]), (
            "w_ld SNPs must match the baseline LD-score SNPs row for row"
        )


def write_prepared_annotation_set(
    directory: Path, annotation_set: PreparedAnnotationSet
) -> None:
    _write_ld_score_table(directory, annotation_set.ld)
    pl.DataFrame(
        {
            _LABEL_COL: [m.label for m in annotation_set.models],
            _TARGET_COL: [m.target_column for m in annotation_set.models],
            _CONTROLS_COL: [list(m.control_columns) for m in annotation_set.models],
        },
        schema={
            _LABEL_COL: pl.String,
            _TARGET_COL: pl.String,
            _CONTROLS_COL: pl.List(pl.String),
        },
    ).write_parquet(directory / MODELS_FILENAME)


def read_prepared_annotation_set(directory: Path) -> PreparedAnnotationSet:
    models = pl.read_parquet(directory / MODELS_FILENAME)
    return PreparedAnnotationSet(
        ld=_read_ld_score_table(directory),
        models=tuple(
            SLDSCAnnotationModel(
                label=row[_LABEL_COL],
                target_column=row[_TARGET_COL],
                control_columns=tuple(row[_CONTROLS_COL]),
            )
            for row in models.iter_rows(named=True)
        ),
    )


def write_prepared_baseline(directory: Path, baseline: PreparedBaseline) -> None:
    _write_ld_score_table(directory, baseline.ld)
    baseline.w_ld.write_parquet(directory / W_LD_FILENAME)


def read_prepared_baseline(directory: Path) -> PreparedBaseline:
    return PreparedBaseline(
        ld=_read_ld_score_table(directory),
        w_ld=pl.read_parquet(directory / W_LD_FILENAME),
    )


def _write_ld_score_table(directory: Path, table: LdScoreTable) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    values = pl.from_numpy(table.matrix, schema=list(table.column_ids), orient="row")
    pl.concat([table.snps, values], how="horizontal").write_parquet(
        directory / LD_SCORES_FILENAME, compression="zstd"
    )
    pl.DataFrame(
        {
            _COLUMN_ID_COL: list(table.column_ids),
            _M_5_50_COL: table.m_5_50,
            _SOURCE_COL: list(table.sources),
        }
    ).write_parquet(directory / M_5_50_FILENAME)


def _read_ld_score_table(directory: Path) -> LdScoreTable:
    m_frame = pl.read_parquet(directory / M_5_50_FILENAME)
    column_ids = tuple(m_frame[_COLUMN_ID_COL].to_list())
    frame = pl.read_parquet(directory / LD_SCORES_FILENAME)
    return LdScoreTable(
        snps=frame.select(_KEY_COLS),
        column_ids=column_ids,
        sources=tuple(m_frame[_SOURCE_COL].to_list()),
        matrix=frame.select(list(column_ids)).to_numpy().astype(np.float32, copy=False),
        m_5_50=m_frame[_M_5_50_COL].to_numpy().astype(np.float64),
    )
