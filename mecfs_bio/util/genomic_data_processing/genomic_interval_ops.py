import narwhals
import structlog

from mecfs_bio.constants.genomic_coordinate_constants import (
    GenomeBuild,
    MHCRegion,
    extended_mhc_interval,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

logger = structlog.get_logger()


def exclude_genomic_interval(
    df: narwhals.LazyFrame, interval: GenomicInterval
) -> narwhals.LazyFrame:
    return df.filter(
        ~(
            (narwhals.col(GWASLAB_CHROM_COL) == interval.chrom)
            & (narwhals.col(GWASLAB_POS_COL) >= interval.start)
            & (narwhals.col(GWASLAB_POS_COL) < interval.end)
        )
    )


def exclude_mhc(
    df: narwhals.LazyFrame,
    build: GenomeBuild = "19",
    region: MHCRegion | None = "extended",
) -> narwhals.LazyFrame:
    if region is None:
        return df
    if region == "extended":
        logger.debug("Excluding the extended MHC region", build=build)
        interval = extended_mhc_interval(build)
    else:
        raise ValueError("Region/build combination not implemented")

    return exclude_genomic_interval(
        df=df,
        interval=interval,
    )
