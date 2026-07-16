"""
Regenerate ppp_manifest.csv by walking the UKB-PPP European (discovery) folder on
Synapse.

The manifest is a static list of every per-protein summary-statistics tar file in the
UK Biobank Pharma Proteomics Project (Sun et al. 2023). Each file is a flat FileEntity
directly under the parent folder, named gene_uniprot_OID_version_panel.tar, so the
filename supplies every column except the Synapse entity id, which comes from the
folder listing.

OID is the primary key: a handful of control proteins are assayed on more than one
Olink panel and share a gene symbol, distinguished only by their OID.

Run via: pixi r python mecfs_bio/assets/reference_data/ukbb_ppp_sumstats/regenerate_ppp_manifest.py

Requires Synapse credentials (see get_file_from_synapse_task). Downloads no data; it
only lists folder metadata. Re-run to pick up any upstream additions, then diff the
result against the committed manifest.
"""

import csv
import re
from pathlib import Path

import synapseclient

# European (discovery) folder holding the 2,940 per-protein tar files.
PPP_DISCOVERY_FOLDER_SYNID = "syn51365303"

MANIFEST_PATH = Path(__file__).parent / "ppp_manifest.csv"

# gene_uniprot_OID_version_panel.tar; uniprot may carry an isoform suffix (e.g. O75882-2)
# and panel names may contain an underscore (e.g. Inflammation_II).
FILENAME_PATTERN = re.compile(r"^(.+)_([A-Za-z0-9-]+)_(OID\d+)_(v\d+)_(.+)\.tar$")

COLUMNS = ["gene", "uniprot", "oid", "version", "panel", "synapse_id", "filename"]


def main() -> None:
    syn = synapseclient.login()
    children = list(syn.getChildren(PPP_DISCOVERY_FOLDER_SYNID, includeTypes=["file"]))

    rows = []
    unparsed = []
    for child in children:
        match = FILENAME_PATTERN.match(child["name"])
        if match is None:
            unparsed.append(child["name"])
            continue
        gene, uniprot, oid, version, panel = match.groups()
        rows.append((gene, uniprot, oid, version, panel, child["id"], child["name"]))

    assert not unparsed, f"filenames did not match expected pattern: {unparsed}"
    oids = {row[2] for row in rows}
    assert len(oids) == len(rows), "OID is not unique across files; it must be the key"

    with MANIFEST_PATH.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(COLUMNS)
        writer.writerows(sorted(rows))

    print(f"wrote {len(rows)} rows to {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
