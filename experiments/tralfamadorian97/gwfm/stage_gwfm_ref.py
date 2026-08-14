from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.build_system.task.sbayesrc.stage_gwfm_reference_task import (
    StageGwfmReferenceTask,
)

BUCKET = "mecfs-bio-reference-data"


def main() -> None:
    task = StageGwfmReferenceTask.create(bucket=BUCKET)
    # A first run has no build trace, so the task executes and stages. To force a
    # re-run after the trace exists, pass must_rebuild_transitive=[task]; the
    # per-file S3 dedup still skips anything already uploaded, so this is cheap.
    DEFAULT_RUNNER.run([task])


if __name__ == "__main__":
    main()
