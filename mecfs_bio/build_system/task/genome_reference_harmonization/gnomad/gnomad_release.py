"""
Descriptor of one gnomAD sites-VCF release: where its per-chromosome files live and which
ancestry groups it carries. Instances live under mecfs_bio/assets/reference_data/gnomad/.
"""

from typing import Literal

from attrs import frozen

from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    gwaslab_code_to_contig_name,
)
from mecfs_bio.constants.allele_frequency_panel_constants import GnomadGroup
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_CODE_FOR_NAME

ContigPrefix = Literal["", "chr"]
CHROM_PLACEHOLDER = "{chrom}"
# Asset-store group shared by every gnomAD asset.
GNOMAD_GROUP = "gnomad"


@frozen(slots=True)
class GnomadRelease:
    """A gnomAD release's per-chromosome sites VCFs.

    vcf_url_template contains CHROM_PLACEHOLDER, replaced by the bare chromosome name
    (1-22, X, Y). contig_prefix is what the VCF prepends to contig names ("chr" in GRCh38
    releases). main_groups admit records (a record is kept when polymorphic in one of
    them); extra_groups are stored but admit nothing. header_assembly is the assembly
    string every ##contig header line must carry.
    """

    name: str
    build: GenomeBuild
    vcf_url_template: str
    contig_prefix: ContigPrefix
    chromosomes: tuple[int, ...]
    main_groups: tuple[GnomadGroup, ...]
    extra_groups: tuple[GnomadGroup, ...]
    header_assembly: str

    def __attrs_post_init__(self):
        assert CHROM_PLACEHOLDER in self.vcf_url_template, (
            f"vcf_url_template must contain {CHROM_PLACEHOLDER}"
        )
        assert self.main_groups, "a release needs at least one main group"
        assert len(set(self.groups)) == len(self.groups), (
            f"groups must be distinct, and main and extra disjoint: {self.groups}"
        )
        assert self.chromosomes, "a release needs at least one chromosome"
        assert list(self.chromosomes) == sorted(set(self.chromosomes)), (
            f"chromosomes must be ascending gwaslab codes without repeats: {self.chromosomes}"
        )
        assert GWASLAB_CHROM_CODE_FOR_NAME["MT"] not in self.chromosomes, (
            "gnomAD panels exclude MT"
        )

    @property
    def groups(self) -> tuple[GnomadGroup, ...]:
        return self.main_groups + self.extra_groups


def gnomad_contig_name(release: GnomadRelease, chrom: int) -> str:
    """The contig name the release's VCF uses for a gwaslab chromosome code."""
    return release.contig_prefix + gwaslab_code_to_contig_name(chrom)


def gnomad_vcf_url(release: GnomadRelease, chrom: int) -> str:
    return release.vcf_url_template.replace(
        CHROM_PLACEHOLDER, gwaslab_code_to_contig_name(chrom)
    )
