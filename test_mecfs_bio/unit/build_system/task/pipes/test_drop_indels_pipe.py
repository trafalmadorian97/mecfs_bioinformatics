import narwhals
import numpy as np
import pandas as pd
import polars as pl

from mecfs_bio.build_system.task.pipes.drop_indels_pipe import DropIndelsPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)


def test_drop_indels_pipe_keeps_single_base_variants_only() -> None:
    frame = pl.DataFrame(
        {
            GWASLAB_POS_COL: [1, 2, 3, 4],
            GWASLAB_EFFECT_ALLELE_COL: ["A", "AT", "G", "GC"],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["C", "A", "GT", "TA"],
        }
    ).with_columns(pl.col(GWASLAB_EFFECT_ALLELE_COL).cast(pl.Categorical))
    result = DropIndelsPipe().process_eager_polars(frame)
    assert result[GWASLAB_POS_COL].to_list() == [1]


def test_drop_indels_pipe_on_pandas_categorical_with_very_long_allele() -> None:
    """
    gwaslab stores alleles as pandas categoricals. Casting one to string allocates rows x
    (longest category) characters, which for this frame would be about 4 TB.
    """
    n_rows = 1_000_000
    codes = np.zeros(n_rows, dtype=np.int8)
    codes[1] = 2
    alleles = pd.Categorical.from_codes(
        codes, categories=pd.Index(["A", "C", "G" * 1_000_000])
    )
    frame = pd.DataFrame(
        {
            GWASLAB_POS_COL: np.arange(n_rows),
            GWASLAB_EFFECT_ALLELE_COL: alleles,
            GWASLAB_NON_EFFECT_ALLELE_COL: pd.Categorical(["C"] * n_rows),
        }
    )
    result = (
        DropIndelsPipe()
        .process(narwhals.from_native(frame).lazy())
        .collect()
        .to_polars()
    )
    assert result.height == n_rows - 1
    assert 1 not in result[GWASLAB_POS_COL]
