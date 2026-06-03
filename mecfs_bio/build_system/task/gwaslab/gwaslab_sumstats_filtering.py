"""
Shared SNP-filtering helpers for GWASLab Sumstats objects.

These were originally written for CT-LDSC (cross-trait genetic correlation), where filtering
roughly consistent with the methods of Bulik-Sullivan et al. is applied before the regression.
They are shared here so that single-trait S-LDSC tasks can apply the same filtering, which is
useful both for methodological consistency and for experiments that probe how SNP-set choices
(e.g. retaining vs. removing strand-ambiguous palindromic variants) affect results.
"""

import gwaslab
import structlog
from attrs import frozen

from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL

logger = structlog.get_logger()


@frozen
class FilterSettings:
    """
    Options for SNP filtering.
    """

    remove_indels: bool = True
    remove_palindromic: bool = True
    remove_hla: bool = True
    keep_only_hapmap: bool = True


def filter_sumstats(sumstats: gwaslab.Sumstats, settings: FilterSettings, build: str):
    """
    Performing filtering roughly consistent with the procedure described in the methods section of Bulik-Sullivan et al.
    """
    if settings.remove_indels:
        logger.debug("filtering indels")
        sumstats.filter_indel(inplace=True, mode="out")
    if settings.remove_palindromic:
        logger.debug("filtering palindromes")
        sumstats.filter_palindromic(inplace=True, mode="out")
    if settings.remove_hla:
        logger.debug("filtering hla region")
        sumstats.exclude_hla(inplace=True)
    if settings.keep_only_hapmap:
        sumstats.filter_hapmap3(inplace=True, build=build)
    logger.debug("dropping duplicate rsids")
    len_before = len(sumstats.data)
    sumstats.data = sumstats.data.drop_duplicates(subset=[GWASLAB_RSID_COL], keep=False)
    len_after = len(sumstats.data)
    logger.debug(f"dropped {len_before - len_after} variants with identical rsids")
