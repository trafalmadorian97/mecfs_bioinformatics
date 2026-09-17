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
