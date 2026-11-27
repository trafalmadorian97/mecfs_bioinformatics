"""
Load the committed deCODE SomaScan manifests (see regenerate_decode_manifest.py) into
DecodeSomascanProtein rows.
"""

import polars as pl

from mecfs_bio.assets.gwas.decode_somascan.raw.decode_s3 import (
    DECODE_BUCKET,
    DECODE_S3_ENDPOINT,
)
from mecfs_bio.assets.reference_data.decode_somascan_sumstats.regenerate_decode_manifest import (
    MANIFEST_ENTREZ_COL,
    MANIFEST_GENE_COL,
    MANIFEST_KEY_COL,
    MANIFEST_MD5_COL,
    MANIFEST_SEQ_ID_COL,
    MANIFEST_SIZE_COL,
    MANIFEST_TARGET_COL,
    MANIFEST_UNIPROT_COL,
    manifest_path,
)
from mecfs_bio.build_system.task.decode_database.decode_somascan_protein import (
    DecodeSomascanProtein,
)
from mecfs_bio.build_system.wf.s3_access import S3ObjectRef
from mecfs_bio.constants.csf_database_constants import EntrezGeneId, SeqId, UniProtId
from mecfs_bio.constants.decode_somascan_constants import DecodeSomascanGroup


def _protein(row: dict) -> DecodeSomascanProtein:
    entrez = row[MANIFEST_ENTREZ_COL]
    return DecodeSomascanProtein(
        seq_id=SeqId(row[MANIFEST_SEQ_ID_COL]),
        gene_symbol=row[MANIFEST_GENE_COL],
        target_name=row[MANIFEST_TARGET_COL],
        uniprot=UniProtId(row[MANIFEST_UNIPROT_COL]),
        entrez_gene_id=None if entrez is None else EntrezGeneId(entrez),
        ref=S3ObjectRef(
            endpoint=DECODE_S3_ENDPOINT, bucket=DECODE_BUCKET, key=row[MANIFEST_KEY_COL]
        ),
        size=row[MANIFEST_SIZE_COL],
        md5=row[MANIFEST_MD5_COL],
    )


def load_decode_somascan_proteins(
    group: DecodeSomascanGroup,
) -> tuple[DecodeSomascanProtein, ...]:
    """Every assay of the group's manifest, in manifest order (sorted by key)."""
    manifest = pl.read_csv(
        manifest_path(group),
        schema_overrides={MANIFEST_ENTREZ_COL: pl.String, MANIFEST_SIZE_COL: pl.Int64},
    )
    return tuple(_protein(row) for row in manifest.iter_rows(named=True))
