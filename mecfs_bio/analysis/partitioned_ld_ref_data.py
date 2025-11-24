from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.raw.partitioned_model_baseline_ld_scores import (
    THOUSAND_GENOME_BASE_MODEL_PARTITIONED_LD_SCORES,
)


def run_initial_analysis():
    DEFAULT_RUNNER.run([THOUSAND_GENOME_BASE_MODEL_PARTITIONED_LD_SCORES])


if __name__ == "__main__":
    run_initial_analysis()
