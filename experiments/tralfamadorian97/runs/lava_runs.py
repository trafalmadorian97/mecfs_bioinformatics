from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.lava.initial_lava import BASIC_G100_LAVA_ANALYSIS
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.lava.initial_lava_ukbb import BASIC_UKBB_LAVA_ANALYSIS
from mecfs_bio.assets.reference_data.lava_ld_reference.g1000_eur.processed.lava_thousand_geomes_eur_ld_ref_extracted import \
    LAVA_G100_EUR_LD_REF_EXTRACTED
from mecfs_bio.assets.reference_data.lava_ld_reference.g1000_eur.raw.thousand_genomes_eur_ld_ref import \
    LAVA_G1000_EUR_REF_LD
from mecfs_bio.assets.reference_data.lava_ld_reference.ukbb.raw.ukbb_lava_ld_ref import LAVA_UKBB_LD_REF
from mecfs_bio.assets.reference_data.lava_locus_file.default.processed.default_lava_locus_file_fused import \
    DEFAULT_LAVA_LOCUS_FUSED
from mecfs_bio.assets.reference_data.lava_locus_file.default.raw.default_lava_locus_file import DEFAULT_LAVA_LOCUS_FILE


def run_lava():
    DEFAULT_RUNNER.run(
        [
            # LAVA_UKBB_LD_REF,
            BASIC_UKBB_LAVA_ANALYSIS
            # LAVA_G100_EUR_LD_REF_EXTRACTED,
            # LAVA_G1000_EUR_REF_LD,
            # DEFAULT_LAVA_LOCUS_FILE
            # BASIC_G100_LAVA_ANALYSIS
            # DEFAULT_LAVA_LOCUS_FUSED
        ],
        incremental_save=True,
        must_rebuild_transitive=[
            BASIC_UKBB_LAVA_ANALYSIS,
            DEFAULT_LAVA_LOCUS_FUSED
            # BASIC_G100_LAVA_ANALYSIS
        ]
    )


if __name__ == "__main__":
    run_lava()
