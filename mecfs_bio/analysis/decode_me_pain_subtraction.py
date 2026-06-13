from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.multi_trait.genomic_sem.decode_me_minus_pain_subtraction import (
    DECODE_ME_MINUS_PAIN_SUBTRACTION,
)


def run_initial_decode_me_minus_pain():
    """ """
    DEFAULT_RUNNER.run(
        [DECODE_ME_MINUS_PAIN_SUBTRACTION],
        must_rebuild_transitive=[],
    )


if __name__ == "__main__":
    run_initial_decode_me_minus_pain()
