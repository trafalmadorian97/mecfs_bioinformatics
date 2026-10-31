"""
Builds the gnomAD v4.1 genomes GRCh38 allele-frequency panel (about 524 GiB streamed, about
14 h). The per-chromosome parts are routed to the external drive by default_runner_config.yaml.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.gnomad.gnomad_allele_frequency_panels import (
    GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES,
)

if __name__ == "__main__":
    DEFAULT_RUNNER.run(
        [GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES],
        incremental_save=True,
    )
