from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.diagrams.ld_same_branch_correlation import LD_SAME_BRANCH_CORRELATION_DIAGRAM


def go():
    DEFAULT_RUNNER.run(
        [
            LD_SAME_BRANCH_CORRELATION_DIAGRAM
        ],
        must_rebuild_transitive=[LD_SAME_BRANCH_CORRELATION_DIAGRAM]
    )

if __name__ == "__main__":
    go()
