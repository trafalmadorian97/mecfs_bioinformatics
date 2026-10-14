"""Build the single sorted allele-bearing annotation matrix from the extracted
members, then report its location. Verifies Task B at real (19.5M-row) scale.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.polyfun.annotations.baseline_lf_annotations import (
    BASELINE_LF_ANNOTATION_MATRIX,
)


def main() -> None:
    store = DEFAULT_RUNNER.run([BASELINE_LF_ANNOTATION_MATRIX], incremental_save=True)
    print("OUTPUT:", store[BASELINE_LF_ANNOTATION_MATRIX.asset_id].path)


if __name__ == "__main__":
    main()
