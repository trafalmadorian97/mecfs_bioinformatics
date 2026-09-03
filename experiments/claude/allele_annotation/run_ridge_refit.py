"""Refit the ridge annotation weights on the allele-bearing matrix with the
allele-aware snpvar join, then print the diagnostics. Verifies Task C at scale.
"""

import json
from pathlib import Path

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.polyfun.annotations.annotation_ridge_weights import (
    BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS,
)
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    DIAGNOSTICS_JSON_FILENAME,
)


def main() -> None:
    store = DEFAULT_RUNNER.run(
        [BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS], incremental_save=True
    )
    path = store[BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS.asset_id].path
    diag = json.loads((Path(path) / DIAGNOSTICS_JSON_FILENAME).read_text())
    print("OUTPUT:", path)
    print("DIAGNOSTICS:", json.dumps(diag, indent=2))


if __name__ == "__main__":
    main()
