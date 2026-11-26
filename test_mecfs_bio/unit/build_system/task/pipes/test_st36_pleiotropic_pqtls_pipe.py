import math

import pandas as pd
import polars as pl
import pytest

from mecfs_bio.build_system.task.pipes.st36_pleiotropic_pqtls_pipe import (
    St36PleiotropicPqtlsPipe,
)
from mecfs_bio.constants.decode_somascan_constants import (
    DECODE_CHR_X_CODE,
    ST36_CHROM_COL,
    ST36_PC0_POS_COL,
    ST36_PLEIOTROPIC_OLINK_COL,
    ST36_PLEIOTROPIC_SMP_COL,
)


def _sheet(smp_flags: list[str | None]) -> pl.DataFrame:
    """
    Two ST36 rows as the extracted sheet holds them: wrapped headers, and empty numeric cells as
    NaN (pandas reads the workbook), not null.
    """
    frame = pd.DataFrame(
        {
            "Chrom": ["chr1", "chrX"],
            "Olink\nPos (hg38)": [100.0, math.nan],
            "Highly\nPleiotropic\nOlink": ["Y", None],
            "Highly\nPleiotropic\nSoma\nSMP": smp_flags,
            "Highly\nPleiotropic\nSoma\nPC0": ["N", "Y"],
            "Soma\nSMP\nPos": [110, 220],
            "Soma\nSMP\nName": ["chr1:110:SG", "chrX:220:SG"],
            "Soma\nSMP\nrsName": ["rs1", "."],
            "Soma\nSMP\nnprot": [50.0, 60.0],
            "Soma\nPC0\nPos": [math.nan, 220.0],
            "Unused\ncolumn": [0.5, 0.25],
        }
    )
    return pl.from_pandas(frame, nan_to_null=False)


def test_types_columns_and_maps_empty_cells() -> None:
    out = St36PleiotropicPqtlsPipe().process_eager_polars(_sheet(["Y", "N"]))
    assert out[ST36_CHROM_COL].to_list() == [1, DECODE_CHR_X_CODE]
    assert out[ST36_PC0_POS_COL].to_list() == [None, 220]
    assert out[ST36_PLEIOTROPIC_OLINK_COL].to_list() == [True, False]
    assert out[ST36_PLEIOTROPIC_SMP_COL].to_list() == [True, False]


def test_rejects_unknown_flag_value() -> None:
    with pytest.raises(AssertionError):
        St36PleiotropicPqtlsPipe().process_eager_polars(_sheet(["Y", "maybe"]))
