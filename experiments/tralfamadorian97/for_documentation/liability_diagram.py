from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.diagrams.liability_threshold_ascertained_sample import \
    LIABILITY_THRESHOLD_ASCERTAINED_SAMPLE_DIAGRAM


def go():
    DEFAULT_RUNNER.run([LIABILITY_THRESHOLD_ASCERTAINED_SAMPLE_DIAGRAM])

if __name__ == '__main__':
    go()