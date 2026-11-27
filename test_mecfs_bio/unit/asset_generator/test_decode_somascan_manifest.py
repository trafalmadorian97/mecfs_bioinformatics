import re

from mecfs_bio.asset_generator.decode_somascan_manifest import (
    load_decode_somascan_proteins,
)
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_s3 import DECODE_PREFIX

_MD5 = re.compile(r"[0-9a-f]{32}")


def test_load_manifest_rows_and_refs() -> None:
    seq_ids = {}
    for group in ("smp", "raw"):
        proteins = load_decode_somascan_proteins(group)
        assert len(proteins) == 4_907
        ids = [p.seq_id for p in proteins]
        assert len(set(ids)) == len(ids)
        assert all(p.ref.key.startswith(DECODE_PREFIX[group]) for p in proteins)
        assert all(_MD5.fullmatch(p.md5) for p in proteins)
        seq_ids[group] = set(ids)
    assert seq_ids["smp"] == seq_ids["raw"]
