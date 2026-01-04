"""
Download a MAGMA gene specific matrix derived from the human brain cell atlas

This matrix was used in the paper:
Duncan, Laramie E., et al. "Mapping the Cellular Etiology of Schizophrenia and Diverse Brain Phenotypes." medRxiv (2024): 2024-10.

source github repo: https://github.com/Integrative-Mental-Health-Lab/linking_cell_types_to_brain_phenotypes?tab=readme-ov-file
"""
from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import ReferenceFileMeta
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MAGMA_ENTREZ_SPECIFICITY_MATRIX_HBCA_RNA_DUNCAN = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="magma_reference_data",
        sub_group="gene_specificity_matrix",
        sub_folder=PurePath("raw"),
        asset_id="magma_specificity_matrix_entrez_hbca_rna_duncan",
        extension=".txt",
        filename="conti_specificity_matrix",
    ),
    url="https://www.dropbox.com/scl/fi/3p5qoyfw5c3q8yf38s0di/conti_specificity_matrix.txt?rlkey=1aoza02mqxn9il5aj3bjn6shq&e=1&dl=1",
    md5_hash="6a422064ef5ac68d3d4beccec8ab230a",
)
