from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER


def go():
    DEFAULT_RUNNER.run(
        [
            KERREBIJN_ET_AL_FIBRO_EUR_RAW
        ]
    )