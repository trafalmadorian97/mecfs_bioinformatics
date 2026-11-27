"""
Regenerate the deCODE SomaScan summary-statistics manifests, one per group (SMP-normalized and
non-normalized).

Each manifest row is one SomaScan assay of ST29 (see st29_soma_assays.py): its seq id, the gene
symbol and target name from the file name, ST29's UniProt accession, the NCBI Entrez gene id(s)
of that accession, and the S3 key, size and md5 of its summary statistics file.

Steps:
  1. list both group prefixes of the deCODE bucket (AWS profile decode-download);
  2. download the Eldjarn et al. 2023 supplementary workbook into a temporary directory
     (checking its md5) and read ST29 from it; keep the files whose seq id is in ST29 (the
     bucket also holds assays the paper did not report), asserting one file per seq id and
     group;
  3. HEAD each kept object for its size and md5 (x-amz-meta-md5chksum);
  4. parse file names, Proteomics_[SMP_]PC0_<seq>_<seq>_<gene>_<target>_<date>.txt.gz;
  5. resolve Entrez gene ids from the UniProt accession through the HGNC complete set (never
     by gene symbol); several ids are joined with "|";
  6. write the CSVs, asserting 4,907 rows per group;
  7. print where ST29's gene label differs from the file name's, and where the CSF aptamer
     manifest gives the same seq id another UniProt accession or gene symbol. These are
     printed for review, not resolved.

The HGNC complete set TSV is an asset; materialize it first and pass its path:
  pixi r python -c "from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER; \
from mecfs_bio.assets.reference_data.hgnc.hgnc_complete_set import HGNC_COMPLETE_SET as T; \
print(DEFAULT_RUNNER.run([T])[T.asset_id].path)"
  pixi r python mecfs_bio/assets/reference_data/decode_somascan_sumstats/regenerate_decode_manifest.py \
      <hgnc_complete_set.tsv>
"""

import base64
import re
import sys
import tempfile
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import boto3
import polars as pl
from botocore.config import Config

from mecfs_bio.assets.gwas.decode_somascan.raw.decode_s3 import (
    DECODE_BUCKET,
    DECODE_PREFIX,
    DECODE_S3_ENDPOINT,
)
from mecfs_bio.assets.reference_data.csf_pqtl_sumstats.regenerate_csf_manifest import (
    MANIFEST_PATH as CSF_MANIFEST_PATH,
)
from mecfs_bio.assets.reference_data.decode_somascan_sumstats.st29_soma_assays import (
    EXPECTED_SOMA_ASSAYS,
    ST29_GENE_COL,
    ST29_SEQ_ID_COL,
    ST29_UNIPROT_COL,
    st29_soma_assays,
)
from mecfs_bio.assets.reference_data.research_paper_supplementary_material.eldjarn_et_al_2023.raw.eldjarn_et_al_2023_supplementary_raw import (
    ELDJARN_ET_AL_2023_SUPPLEMENTARY_MD5,
    ELDJARN_ET_AL_2023_SUPPLEMENTARY_URL,
)
from mecfs_bio.constants.decode_somascan_constants import DecodeSomascanGroup
from mecfs_bio.constants.hgnc_constants import (
    HGNC_LIST_SEPARATOR,
    HGNC_LOCUS_TYPE_COL,
    HGNC_READTHROUGH_LOCUS_TYPE,
    HGNC_UNIPROT_IDS_COL,
)
from mecfs_bio.util.download.verify import verify_hash

_DIR = Path(__file__).parent

# Manifest columns.
MANIFEST_SEQ_ID_COL = "seq_id"
MANIFEST_GENE_COL = "gene_symbol"
MANIFEST_TARGET_COL = "target_name"
MANIFEST_UNIPROT_COL = "uniprot"
# NCBI Entrez gene id(s) of the UniProt accession, "|"-joined when several; empty when HGNC
# gives none.
MANIFEST_ENTREZ_COL = "entrez_gene_id"
MANIFEST_KEY_COL = "key"
MANIFEST_SIZE_COL = "size"
MANIFEST_MD5_COL = "md5"
MANIFEST_COLUMNS = [
    MANIFEST_SEQ_ID_COL,
    MANIFEST_GENE_COL,
    MANIFEST_TARGET_COL,
    MANIFEST_UNIPROT_COL,
    MANIFEST_ENTREZ_COL,
    MANIFEST_KEY_COL,
    MANIFEST_SIZE_COL,
    MANIFEST_MD5_COL,
]


def manifest_path(group: DecodeSomascanGroup) -> Path:
    return _DIR / f"decode_somascan_{group}_manifest.csv"


_FILE_NAME = re.compile(
    r"Proteomics_(?:SMP_)?PC0_(\d+)_(\d+)_([^_]+)_(.+)_\d{8}\.txt\.gz"
)
_HGNC_ENTREZ_COL = "entrez_id"
_HEAD_WORKERS = 32


def _client():  # noqa: ANN202
    session = boto3.Session(profile_name=DECODE_S3_ENDPOINT.credentials_profile)
    return session.client(
        "s3",
        endpoint_url=DECODE_S3_ENDPOINT.url,
        region_name=DECODE_S3_ENDPOINT.region,
        config=Config(
            signature_version="s3v4",
            s3={"addressing_style": "path"},
            max_pool_connections=_HEAD_WORKERS,
        ),
    )


def list_keys(client, prefix: str) -> list[str]:  # noqa: ANN001
    keys: list[str] = []
    for page in client.get_paginator("list_objects_v2").paginate(
        Bucket=DECODE_BUCKET, Prefix=prefix
    ):
        keys.extend(o["Key"] for o in page.get("Contents", []))
    return keys


def parse_keys(keys: list[str]) -> pl.DataFrame:
    rows = []
    for key in keys:
        match = _FILE_NAME.fullmatch(Path(key).name)
        assert match is not None, f"unexpected file name {key}"
        rows.append(
            {
                MANIFEST_SEQ_ID_COL: f"{match.group(1)}-{match.group(2)}",
                MANIFEST_GENE_COL: match.group(3),
                MANIFEST_TARGET_COL: match.group(4),
                MANIFEST_KEY_COL: key,
            }
        )
    return pl.DataFrame(rows)


def head_size_md5(client, keys: list[str]) -> pl.DataFrame:  # noqa: ANN001
    def one(key: str) -> dict[str, object]:
        head = client.head_object(Bucket=DECODE_BUCKET, Key=key)
        md5 = base64.b64decode(head["Metadata"]["md5chksum"]).hex()
        return {
            MANIFEST_KEY_COL: key,
            MANIFEST_SIZE_COL: int(head["ContentLength"]),
            MANIFEST_MD5_COL: md5,
        }

    with ThreadPoolExecutor(_HEAD_WORKERS) as pool:
        return pl.DataFrame(list(pool.map(one, keys)))


def entrez_by_uniprot(hgnc_complete_set: Path) -> pl.DataFrame:
    """One row per UniProt accession: its genes' Entrez ids, sorted and joined with "|"."""
    return (
        pl.read_csv(hgnc_complete_set, separator="\t", infer_schema=False)
        .filter(pl.col(HGNC_LOCUS_TYPE_COL) != HGNC_READTHROUGH_LOCUS_TYPE)
        .select(
            pl.col(HGNC_UNIPROT_IDS_COL)
            .str.split(HGNC_LIST_SEPARATOR)
            .alias(MANIFEST_UNIPROT_COL),
            pl.col(_HGNC_ENTREZ_COL),
        )
        .explode(MANIFEST_UNIPROT_COL)
        .drop_nulls()
        .unique()
        .group_by(MANIFEST_UNIPROT_COL)
        .agg(pl.col(_HGNC_ENTREZ_COL).sort().str.join("|").alias(MANIFEST_ENTREZ_COL))
    )


def group_manifest(
    client,  # noqa: ANN001
    group: DecodeSomascanGroup,
    st29: pl.DataFrame,
    entrez: pl.DataFrame,
) -> pl.DataFrame:
    files = parse_keys(list_keys(client, DECODE_PREFIX[group]))
    kept = files.join(
        st29.select(ST29_SEQ_ID_COL, ST29_UNIPROT_COL), on=MANIFEST_SEQ_ID_COL
    )
    assert kept[MANIFEST_SEQ_ID_COL].is_unique().all(), "several files for one seq id"
    print(f"{group}: {files.height} files listed, {kept.height} in ST29")
    sizes = head_size_md5(client, kept[MANIFEST_KEY_COL].to_list())
    manifest = (
        kept.join(sizes, on=MANIFEST_KEY_COL)
        .join(entrez, on=MANIFEST_UNIPROT_COL, how="left")
        .select(MANIFEST_COLUMNS)
        .sort(MANIFEST_KEY_COL)
    )
    assert manifest.height == EXPECTED_SOMA_ASSAYS, manifest.height
    print(
        f"{group}: {manifest[MANIFEST_ENTREZ_COL].null_count()} rows without an Entrez id"
    )
    return manifest


def report_disagreements(manifest: pl.DataFrame, st29: pl.DataFrame) -> None:
    gene = manifest.join(
        st29.select(MANIFEST_SEQ_ID_COL, pl.col(ST29_GENE_COL).alias("st29_gene")),
        on=MANIFEST_SEQ_ID_COL,
    ).filter(pl.col(MANIFEST_GENE_COL) != pl.col("st29_gene"))
    print(f"ST29 gene label differs from the file name's gene: {gene.height}")
    print(gene.select(MANIFEST_SEQ_ID_COL, MANIFEST_GENE_COL, "st29_gene").head(30))
    csf = pl.read_csv(CSF_MANIFEST_PATH, infer_schema=False).select(
        pl.col("seq_id").alias(MANIFEST_SEQ_ID_COL),
        pl.col("uniprot").alias("csf_uniprot"),
        pl.col("entrez_gene_symbol").alias("csf_gene"),
    )
    both = manifest.join(csf, on=MANIFEST_SEQ_ID_COL)
    uniprot = both.filter(pl.col(MANIFEST_UNIPROT_COL) != pl.col("csf_uniprot"))
    csf_gene = both.filter(pl.col(MANIFEST_GENE_COL) != pl.col("csf_gene"))
    print(f"seq ids also in the CSF manifest: {both.height}")
    print(f"  UniProt differs from CSF: {uniprot.height}")
    print(
        uniprot.select(MANIFEST_SEQ_ID_COL, MANIFEST_UNIPROT_COL, "csf_uniprot").head(
            30
        )
    )
    print(f"  gene symbol differs from CSF: {csf_gene.height}")
    print(csf_gene.select(MANIFEST_SEQ_ID_COL, MANIFEST_GENE_COL, "csf_gene").head(30))


def download_st29() -> pl.DataFrame:
    """ST29's SomaScan assays, from the supplementary workbook downloaded and then discarded."""
    with tempfile.TemporaryDirectory() as tmp:
        workbook = Path(tmp) / "eldjarn_et_al_2023_supplementary.xlsx"
        urllib.request.urlretrieve(ELDJARN_ET_AL_2023_SUPPLEMENTARY_URL, workbook)
        verify_hash(workbook, ELDJARN_ET_AL_2023_SUPPLEMENTARY_MD5)
        return st29_soma_assays(workbook)


def main() -> None:
    pl.Config.set_tbl_rows(40)
    entrez = entrez_by_uniprot(Path(sys.argv[1]))
    st29 = download_st29()
    client = _client()
    manifests: dict[DecodeSomascanGroup, pl.DataFrame] = {
        group: group_manifest(client, group, st29, entrez) for group in ("smp", "raw")
    }
    assert set(manifests["smp"][MANIFEST_SEQ_ID_COL]) == set(
        manifests["raw"][MANIFEST_SEQ_ID_COL]
    )
    for group, manifest in manifests.items():
        manifest.write_csv(manifest_path(group))
        print(f"wrote {manifest.height} rows to {manifest_path(group)}")
    report_disagreements(manifests["smp"], st29)


if __name__ == "__main__":
    main()
