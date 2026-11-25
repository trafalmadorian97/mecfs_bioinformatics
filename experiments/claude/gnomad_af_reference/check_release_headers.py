"""
Check every per-chromosome URL of both gnomAD releases against its release description.

Runs the gnomAD Task's own header check (INFO fields for every group, contig assembly,
contig declared) on each chromosome's live header, so a wrong URL template or contig name
fails here in minutes rather than hours into a build.

Usage:
    pixi r python -m experiments.claude.gnomad_af_reference.check_release_headers \
        2>&1 | tee experiments/claude/gnomad_af_reference/check_release_headers.log
"""

from mecfs_bio.assets.reference_data.gnomad.gnomad_releases import (
    GNOMAD_V2_1_1_GENOMES,
    GNOMAD_V4_1_GENOMES,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_chromosome_allele_frequency_task import (
    assert_gnomad_header,
    read_vcf_header,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_release import (
    gnomad_vcf_url,
)


def main() -> None:
    for release in [GNOMAD_V2_1_1_GENOMES, GNOMAD_V4_1_GENOMES]:
        for chrom in release.chromosomes:
            url = gnomad_vcf_url(release, chrom)
            assert_gnomad_header(read_vcf_header(url, max_attempts=3), release, chrom)
            print(f"ok {release.name} chromosome {chrom} {url}")


if __name__ == "__main__":
    main()
