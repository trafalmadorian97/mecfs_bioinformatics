from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.gnomad.gnomad_allele_frequency_panels import \
    GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES, GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES


def go():
    DEFAULT_RUNNER.run(
        [
         GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES

         ]
    )

if __name__ == '__main__':
    go()