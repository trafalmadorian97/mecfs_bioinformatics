import narwhals
import polars as pl
import polars.testing

from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL
from mecfs_bio.util.genomic_data_processing.genomic_interval_ops import exclude_mhc


def test_exclude_mhc():
    df = pl.DataFrame(
        {GWASLAB_POS_COL: [1, 10, 30726291, 30726291], GWASLAB_CHROM_COL: [6, 6, 6, 1]}
    )
    expected = pl.DataFrame(
        {GWASLAB_POS_COL: [1, 10, 30726291], GWASLAB_CHROM_COL: [6, 6, 1]}
    )
    df_nw = narwhals.from_native(df).lazy()
    result = exclude_mhc(df_nw).collect().to_polars()
    polars.testing.assert_frame_equal(result, expected)
