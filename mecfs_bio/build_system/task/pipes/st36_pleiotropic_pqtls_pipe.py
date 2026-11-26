"""
Normalize sheet ST36_pleiotropic_pQTLs of Eldjarn et al. 2023 (highly pleiotropic pQTLs, one
row per pQTL found by Olink, SomaScan SMP or SomaScan PC0) as read from the workbook.

Keeps the chromosome as a code (chrX = DECODE_CHR_X_CODE), the Olink, SMP and PC0 lead
positions (hg38, null when that platform has no lead), the three pleiotropy flags as booleans
(Y/N in the sheet, empty read as false), and the SMP lead's name, rsName and number of
proteins. The sheet's headers wrap onto several lines; they are matched with newlines replaced
by spaces.
"""

import narwhals
import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.decode_somascan_constants import (
    DECODE_CHR_X_CODE,
    ST36_CHROM_COL,
    ST36_OLINK_POS_COL,
    ST36_PC0_POS_COL,
    ST36_PLEIOTROPIC_OLINK_COL,
    ST36_PLEIOTROPIC_PC0_COL,
    ST36_PLEIOTROPIC_SMP_COL,
    ST36_SMP_NAME_COL,
    ST36_SMP_NPROT_COL,
    ST36_SMP_POS_COL,
    ST36_SMP_RSNAME_COL,
)

_SOURCE_COLUMNS = {
    "Chrom": ST36_CHROM_COL,
    "Olink Pos (hg38)": ST36_OLINK_POS_COL,
    "Highly Pleiotropic Olink": ST36_PLEIOTROPIC_OLINK_COL,
    "Highly Pleiotropic Soma SMP": ST36_PLEIOTROPIC_SMP_COL,
    "Highly Pleiotropic Soma PC0": ST36_PLEIOTROPIC_PC0_COL,
    "Soma SMP Pos": ST36_SMP_POS_COL,
    "Soma SMP Name": ST36_SMP_NAME_COL,
    "Soma SMP rsName": ST36_SMP_RSNAME_COL,
    "Soma SMP nprot": ST36_SMP_NPROT_COL,
    "Soma PC0 Pos": ST36_PC0_POS_COL,
}
_FLAGS = (
    ST36_PLEIOTROPIC_OLINK_COL,
    ST36_PLEIOTROPIC_SMP_COL,
    ST36_PLEIOTROPIC_PC0_COL,
)
_POSITIONS = (ST36_OLINK_POS_COL, ST36_SMP_POS_COL, ST36_PC0_POS_COL)
_FLAG_VALUES = {"Y", "N"}


def _as_int(col: str, dtype: pl.DataType) -> pl.Expr:
    """A numeric sheet column as integers; a float column's NaN (an empty cell) becomes null."""
    return pl.col(col).cast(pl.Float64).fill_nan(None).cast(dtype)


@frozen(slots=True)
class St36PleiotropicPqtlsPipe(DataProcessingPipe):
    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        sheet = x.collect().to_polars()
        renamed = (
            sheet.rename(lambda c: c.replace("\n", " "))
            .select(list(_SOURCE_COLUMNS))
            .rename(_SOURCE_COLUMNS)
            .with_columns(pl.col(c).cast(pl.String) for c in _FLAGS)
        )
        flag_values = set(
            renamed.select(pl.concat_list(_FLAGS).explode()).to_series().drop_nulls()
        )
        assert flag_values <= _FLAG_VALUES, f"unexpected ST36 flag values {flag_values}"
        table = renamed.with_columns(
            pl.col(ST36_CHROM_COL)
            .str.strip_prefix("chr")
            .replace("X", str(DECODE_CHR_X_CODE))
            .cast(pl.Int32),
            *(_as_int(c, pl.Int64()) for c in _POSITIONS),
            _as_int(ST36_SMP_NPROT_COL, pl.Int32()),
            *((pl.col(c) == "Y").fill_null(False) for c in _FLAGS),
        )
        smp = table.filter(pl.col(ST36_PLEIOTROPIC_SMP_COL))
        assert smp[ST36_SMP_POS_COL].null_count() == 0, (
            "an SMP-pleiotropic row has no SMP Pos"
        )
        return narwhals.from_native(table).lazy()
