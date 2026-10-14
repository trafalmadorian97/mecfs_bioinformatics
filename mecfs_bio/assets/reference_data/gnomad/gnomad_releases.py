"""
gnomAD genome sites-VCF releases used for allele-frequency panels.

Files are in the public bucket gnomad-public-us-east-1 (anonymous access, no egress charge
to us). Groups and header facts were measured on chr21; see
experiments/claude/design_specs/2026-10-02-gnomad-allele-frequency-panel-design.md.
"""

from typing import get_args

from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_release import (
    GnomadRelease,
)
from mecfs_bio.constants.allele_frequency_panel_constants import (
    GnomadV2Group,
    GnomadV4Group,
)

_AUTOSOMES = tuple(range(1, 23))
_X = 23
_Y = 24

GNOMAD_V2_1_1_GENOMES = GnomadRelease(
    name="v2_1_1_genomes",
    build="19",
    vcf_url_template=(
        "https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/vcf/genomes/"
        "gnomad.genomes.r2.1.1.sites.{chrom}.vcf.bgz"
    ),
    contig_prefix="",
    chromosomes=(*_AUTOSOMES, _X),
    main_groups=("afr", "amr", "asj", "eas", "fin", "nfe"),
    extra_groups=("oth", "nfe_nwe", "nfe_seu", "nfe_onf", "nfe_est"),
    header_assembly="gnomAD_GRCh37",
)

GNOMAD_V4_1_GENOMES = GnomadRelease(
    name="v4_1_genomes",
    build="38",
    vcf_url_template=(
        "https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/genomes/"
        "gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz"
    ),
    contig_prefix="chr",
    chromosomes=(*_AUTOSOMES, _X, _Y),
    main_groups=("afr", "amr", "asj", "eas", "fin", "nfe", "sas"),
    extra_groups=("ami", "mid", "remaining"),
    header_assembly="gnomAD_GRCh38",
)

# The release group lists and the per-release Literal types must describe the same groups.
for _release, _literal in [
    (GNOMAD_V2_1_1_GENOMES, GnomadV2Group),
    (GNOMAD_V4_1_GENOMES, GnomadV4Group),
]:
    assert set(_release.groups) == set(get_args(_literal)), (
        f"{_release.name} groups {sorted(_release.groups)} differ from "
        f"{sorted(get_args(_literal))}"
    )
