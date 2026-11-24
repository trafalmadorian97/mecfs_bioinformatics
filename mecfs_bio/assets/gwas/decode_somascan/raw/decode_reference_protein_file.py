"""
A local copy of one deCODE SomaScan SMP sumstats file, AP4M1 (seq id 10076-1), kept as a
reference for exploratory checks: brute-force comparisons, row-set and allele-coding surveys.
Every file has (nearly) the same rows, so one file answers most questions about the layout.
"""

from mecfs_bio.asset_generator.decode_somascan_manifest import (
    load_decode_somascan_proteins,
)
from mecfs_bio.build_system.task.decode_database.download_decode_somascan_protein_task import (
    DownloadDecodeSomascanProteinTask,
)

_SEQ_ID = "10076-1"

(_PROTEIN,) = [p for p in load_decode_somascan_proteins("smp") if p.seq_id == _SEQ_ID]

DECODE_SMP_AP4M1_RAW_FILE = DownloadDecodeSomascanProteinTask.create(
    asset_id="decode_somascan_smp_ap4m1_10076_1_raw_file", protein=_PROTEIN
)
