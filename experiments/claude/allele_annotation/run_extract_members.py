"""Run the streaming annotation-parquet extractor with DEFAULT_RUNNER and report
where the output landed. Verification step for the allele-aware annotation rework.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.polyfun.annotations.baseline_lf_annotations import (
    BASELINE_LF_ANNOTATION_PARQUET_MEMBERS,
)


def main() -> None:
    store = DEFAULT_RUNNER.run(
        [BASELINE_LF_ANNOTATION_PARQUET_MEMBERS], incremental_save=True
    )
    asset = store[BASELINE_LF_ANNOTATION_PARQUET_MEMBERS.asset_id]
    print("OUTPUT_DIR:", asset.path)


if __name__ == "__main__":
    main()
