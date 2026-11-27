"""
Entry point that builds the HapMap3 deCODE SomaScan database for the SMP group: the variant info
table, the shared index, every per-protein directory (4,907) and the hit table.

Each per-protein task downloads its ~970 MB file and discards it after one pass, so the full
build takes about 3 days at ~54 s per protein. The AWS profile decode-download must hold deCODE's
keys (aws configure --profile decode-download).
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_hit_table import (
    HAPMAP_3_DECODE_SMP_HIT_TABLE_5E8,
)


def main() -> None:
    DEFAULT_RUNNER.run([HAPMAP_3_DECODE_SMP_HIT_TABLE_5E8])


if __name__ == "__main__":
    main()
