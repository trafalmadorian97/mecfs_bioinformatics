from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me_minus_pain.analysis.residual_genetic_corr import (
    DECODE_ME_MINUS_PAIN_GENETIC_CORR_GENERATOR,
)


def run_initial_decode_me_minus_pain():
    """
    Function to use GenomicSEM to subtract Johnston et al.'s pain GWAS from DecodeME.

    """
    DEFAULT_RUNNER.run(
        # [DECODE_ME_MINUS_PAIN_SUBTRACTION]
        # + DECODE_ME_MINUS_PAIN_STANDARD_ANALYSIS.get_terminal_tasks(),
        DECODE_ME_MINUS_PAIN_GENETIC_CORR_GENERATOR.terminal_tasks(),
        must_rebuild_transitive=[
            # DECODE_ME_MINUS_PAIN_STANDARD_ANALYSIS.magma_tasks.sumstats_task
        ],
    )


if __name__ == "__main__":
    run_initial_decode_me_minus_pain()
