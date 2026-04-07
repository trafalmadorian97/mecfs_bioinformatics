import narwhals
import structlog

from mecfs_bio.constants.genomic_coordinate_constants import (
    EXTENDED_MHC_BUILD_37,
    GenomeBuild,
    MHCRegion,
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
    if region == "extended" and build == "19":
        logger.debug("Excluding the build-19 extended MHC region")
        interval = EXTENDED_MHC_BUILD_37
    else:
        raise ValueError("Region/build combination not implemented")

    return exclude_genomic_interval(
        df=df,
        interval=interval,
    )
