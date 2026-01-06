from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.raw.raw_bellenguez_data import BELLENGUEZ_ET_AL_ALZHIEMERS_RAW


def run_initial_alzhiemers_analysis():
    """
    Function to analyze GWAS of Alzheimer's phenotype.
    Includes:
    - S-LDSC analysis
    - MAGMA analysis (Using both GTEx and HBA reference data)
    """
    DEFAULT_RUNNER.run(
        [
            BELLENGUEZ_ET_AL_ALZHIEMERS_RAW
        ],
        )


if __name__ == "__main__":
    run_initial_alzhiemers_analysis()
