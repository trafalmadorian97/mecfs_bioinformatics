"""
Genome-reference harmonization of the Million Veteran Program EUR ankylosing spondylitis
sumstats (lifted over to build 37), against the UCSC hg19 genome and the gnomAD v2.1.1
genomes panel, so that NEA is the plus-strand reference allele. Ambiguous variants are
resolved with the gnomAD non-Finnish European (nfe) allele frequencies: US veterans of
European ancestry are not specific to any one nfe subgroup.

The source reports odds ratios and p-values, so BETA and SE are derived before
harmonization. Duplicated variant keys, likely a consequence of liftover, are reduced
to the row with the smallest SE.
"""

from mecfs_bio.assets.gwas.ankylosing_spondylitis.million_veterans.processed.mv_eur_ank_spond_sumstats_dump_to_parquet import (
    MILLION_VETERANS_ANK_SPOND_SUMSTATS_37_DUMP_TO_PARQUET,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.gnomad.gnomad_allele_frequency_panels import (
    GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES_REHOSTED,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)

MV_EUR_ANK_SPOND_GENOME_REFERENCE_HARMONIZED = GenomeReferenceHarmonizationTask.create(
    asset_id="million_veterans_ank_spond_genome_reference_harmonized",
    sumstats_task=MILLION_VETERANS_ANK_SPOND_SUMSTATS_37_DUMP_TO_PARQUET,
    fasta_task=UCSC_HG19_INDEXED_FASTA,
    panel_task=GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES_REHOSTED,
    panel_ancestry="nfe",
    pipe=CompositePipe(
        [
            ComputeBetaPipe(),
            ComputeSEPipe(out_backend="polars"),
            UniquePipe(
                by=[
                    GWASLAB_CHROM_COL,
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL,
                ],
                keep="first",
                order_by=[GWASLAB_SE_COL],
            ),
        ]
    ),
)
