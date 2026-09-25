"""THROWAWAY spike: inspect the schema of the polyfun bundle's
.l2.ldscore.parquet members.

Gates decision 7 of the Approach-2 spec: do the LD-score members carry
(chrom, pos, nea, ea) so we can join on the exact four-part key, or only rsid?

gzip is not seekable, so we stream the ~30GB bundle sequentially (tarfile r|gz),
log member names/sizes in archive order, extract ONLY the first
.l2.ldscore.parquet member to a scratch file, read its parquet schema + a couple
of rows, then STOP (break) so we never download the full 30GB.

Run:
  pixi r python experiments/claude/polyfun_approach2_ldscore_schema_spike/probe_ldscore_schema.py \
    2>&1 | tee experiments/claude/polyfun_approach2_ldscore_schema_spike/probe.log
"""

import re
import shutil
import tarfile
import urllib.request
from pathlib import Path

import pyarrow.parquet as pq

URL = "https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.polyfun.tar.gz"
LDSCORE_RE = re.compile(r"baselineLF2\.2\.UKB\.(\d+)\.l2\.ldscore\.parquet$")

# Scratch file for the single extracted member (session scratchpad, not the repo).
SCRATCH = Path(
    "/tmp/claude-1000/-home-paiforsyth-src-traf-biostatistics/"
    "6cf5e872-6825-4859-919f-3c44bcae1278/scratchpad"
)
SCRATCH.mkdir(parents=True, exist_ok=True)


def main() -> None:
    print(f"streaming {URL}", flush=True)
    extracted: Path | None = None
    with urllib.request.urlopen(URL) as raw:
        with tarfile.open(fileobj=raw, mode="r|gz") as tar:
            for member in tar:
                if not member.isfile():
                    continue
                print(f"  member: {member.name}  size={member.size}", flush=True)
                if LDSCORE_RE.search(member.name) is not None:
                    dest = SCRATCH / Path(member.name).name
                    src = tar.extractfile(member)
                    assert src is not None
                    with open(dest, "wb") as out:
                        shutil.copyfileobj(src, out)
                    extracted = dest
                    print(f"  -> extracted {dest} ({dest.stat().st_size} bytes)", flush=True)
                    break  # stop: we only need one member's footer

    assert extracted is not None, "no .l2.ldscore.parquet member found in bundle"

    pf = pq.ParquetFile(extracted)
    print("\n=== SCHEMA ===", flush=True)
    print(pf.schema_arrow, flush=True)
    names = pf.schema_arrow.names
    print(f"\nn_columns = {len(names)}", flush=True)
    print(f"first 12 column names: {names[:12]}", flush=True)
    # Key-column presence check (case-insensitive contains).
    lowered = [n.lower() for n in names]
    for want in ("chr", "bp", "pos", "snp", "a1", "a2", "allele"):
        hits = [n for n, lo in zip(names, lowered) if want in lo]
        print(f"  columns containing '{want}': {hits}", flush=True)

    print("\n=== HEAD (first 3 rows, first 8 cols) ===", flush=True)
    head = pf.read_row_group(0, columns=names[: min(8, len(names))]).slice(0, 3)
    print(head.to_pydict(), flush=True)


if __name__ == "__main__":
    main()
