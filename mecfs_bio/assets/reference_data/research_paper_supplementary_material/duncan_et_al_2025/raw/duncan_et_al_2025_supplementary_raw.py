"""
Task to download the supplementary material from the paper Duncan, Laramie E., et al. "Mapping the cellular etiology of schizophrenia and complex brain phenotypes."Nature Neuroscience" 28.2 (2025): 248-258.

May be useful for interpreting the results of MAGMA runs that use the human brain atlas reference data.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

DUNCAN_ET_AL_2025_SUPPLEMENTARY_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="research_paper_supplementary_material",
        sub_group="duncan_et_al_2025",
        sub_folder=PurePath("raw"),
        asset_id="duncan_el_al_2025_raw",
        filename="41593_2024_1834_MOESM2_ESM",
        extension=".xlsx",
    ),
    url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-024-01834-w/MediaObjects/41593_2024_1834_MOESM2_ESM.xlsx",
    md5_hash="5774e0a767514b086a49f1b0fe59893b",
)
