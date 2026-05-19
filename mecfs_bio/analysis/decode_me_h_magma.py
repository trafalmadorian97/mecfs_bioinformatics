from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.h_magma.decode_me_h_magma_asset_generator import (
    DECODE_ME_H_MAGMA_ASSET_GENERATOR,
)


def run_decode_me_h_magma():
    """
    Function to run H-MAGMA analysis of DecodeME data

    """
    DEFAULT_RUNNER.run(
        DECODE_ME_H_MAGMA_ASSET_GENERATOR.terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[],
    )


if __name__ == "__main__":
    run_decode_me_h_magma()
