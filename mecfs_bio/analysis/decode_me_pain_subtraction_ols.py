from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me_minus_pain.analysis.residual_genetic_corr_ols import (
    DECODE_ME_MINUS_PAIN_OLS_GENETIC_CORR_GENERATOR,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_minus_pain.analysis.standard_analysis_decodeme_minus_pain_ols import (
    DECODE_ME_MINUS_PAIN_OLS_STANDARD_ANALYSIS,
)


def run_decode_me_minus_pain_ols():
    """
    Linear-scale (OLS) variant of GWAS-by-subtraction: DecodeME is treated as a
    quantitative/OLS trait so the subtraction coefficient and the per-SNP betas
    share a scale, making the remainder genuinely orthogonal to pain.

    Builds the genetic-correlation aggregation comparing the OLS remainder factor
    against DecodeME and multisite pain.
    """
    DEFAULT_RUNNER.run(
        DECODE_ME_MINUS_PAIN_OLS_GENETIC_CORR_GENERATOR.terminal_tasks()
        + DECODE_ME_MINUS_PAIN_OLS_STANDARD_ANALYSIS.get_terminal_tasks()
    )


if __name__ == "__main__":
    run_decode_me_minus_pain_ols()
