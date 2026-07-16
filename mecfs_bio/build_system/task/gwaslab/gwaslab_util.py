import json
from pathlib import Path

import gwaslab as gl
import pandas as pd
import structlog
from attrs import frozen

from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)
from mecfs_bio.util.download.verify import calc_md5_checksum
from mecfs_bio.util.plotting.save_fig import normalize_filename

logger = structlog.get_logger()

_REFERENCE_MD5_KEY = "md5sum"


@frozen
class Variant:
    """
    Represents a genetic variant.
    """

    chromosome: int
    position: int
    effect_allele: str
    non_effect_allele: str

    @property
    def id(self) -> str:
        return f"{self.chromosome}:{self.position}:{self.non_effect_allele}:{self.effect_allele}"

    @property
    def id_normalized(self) -> str:
        return normalize_filename(self.id).replace(":", "_")


def df_to_variants(df: pd.DataFrame) -> list[Variant]:
    return [
        Variant(
            chromosome=df[GWASLAB_CHROM_COL].iloc[i],
            position=df[GWASLAB_POS_COL].iloc[i],
            effect_allele=df[GWASLAB_EFFECT_ALLELE_COL].iloc[i],
            non_effect_allele=df[GWASLAB_NON_EFFECT_ALLELE_COL].iloc[i],
        )
        for i in range(len(df))
    ]


def _reference_catalogue() -> dict:
    """The reference metadata bundled with the installed gwaslab."""
    catalogue_path = Path(gl.__file__).parent / "data" / "reference.json"
    with open(catalogue_path) as handle:
        return json.load(handle)


def expected_reference_md5(ref: str) -> str | None:
    """The checksum gwaslab records for a reference, or None if it records none.

    gwaslab leaves this blank for many entries (the genome fastas, the GTFs), so a
    None here means unverifiable, not invalid.
    """
    catalogue = _reference_catalogue()
    assert ref in catalogue, f"unknown gwaslab reference {ref}"
    return catalogue[ref].get(_REFERENCE_MD5_KEY) or None


def _local_reference_path(ref: str) -> Path | None:
    """Where gwaslab has this reference on disk, or None if it does not have it."""
    result = gl.get_path(ref)
    if not result:
        return None
    assert isinstance(result, str), f"unexpected gwaslab path for {ref}: {result!r}"
    return Path(result)


def gwaslab_download_ref_if_missing(ref: str) -> Path:
    """Return the local path to a gwaslab reference, downloading it when absent and
    re-downloading it when its checksum does not match the one gwaslab records.
    """
    expected_md5 = expected_reference_md5(ref)

    local_path = _local_reference_path(ref)
    if local_path is None:
        # NOTE: Gwaslab store paths of downloaded files in a "config" object stored within the gwaslab package
        # This object can grow out of date, meaning that the file exists but is not referenced in the config
        # Besides actually downloading missing files, calling "download ref" can also sometimes update the config to point to
        # existing local file that are not in the config
        gl.download_ref(ref)
        local_path = _local_reference_path(ref)
        assert local_path is not None, f"gwaslab failed to download reference {ref}"

    if expected_md5 is None:
        logger.debug(f"ref {ref} has no expected checksum")
        return local_path
    if calc_md5_checksum(local_path) == expected_md5:
        logger.debug(f"ref {ref} matches the expected checksum")
        return local_path

    logger.warning(
        "cached gwaslab reference does not match its recorded checksum, re-downloading",
        ref=ref,
        path=str(local_path),
    )
    gl.download_ref(ref, overwrite=True)

    refreshed_path = _local_reference_path(ref)
    assert refreshed_path is not None, f"gwaslab failed to download reference {ref}"
    actual_md5 = calc_md5_checksum(refreshed_path)
    assert actual_md5 == expected_md5, (
        f"gwaslab reference {ref} at {refreshed_path} has md5 {actual_md5}, but "
        f"gwaslab records {expected_md5}, even after a fresh download"
    )
    return refreshed_path
