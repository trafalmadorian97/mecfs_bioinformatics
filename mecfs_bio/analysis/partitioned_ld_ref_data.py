from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_multi_tissue_gene_expr_ld_score_extracted import (
    PARTITIONED_MODEL_MULTI_TISSUE_GENE_EXPR_LD_SCORES_EXTRACTED,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_regression_weights_extracted import (
    PARTITIONED_MODEL_REGRESSION_WEIGHTS_EXTRACTED,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.thousand_genome_baseline_ld_extracted import (
    BASE_MODEL_PARTITIONED_LD_SCORES_EXTRACTED,
)


def run_initial_analysis():
    DEFAULT_RUNNER.run(
        [
            # THOUSAND_GENOME_BASE_MODEL_PARTITIONED_LD_SCORES,
            PARTITIONED_MODEL_MULTI_TISSUE_GENE_EXPR_LD_SCORES_EXTRACTED,
            BASE_MODEL_PARTITIONED_LD_SCORES_EXTRACTED,
            PARTITIONED_MODEL_REGRESSION_WEIGHTS_EXTRACTED,
        ]
    )


if __name__ == "__main__":
    run_initial_analysis()
