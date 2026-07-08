"""
Utilities for invoking the POPs command-line scripts.

The POPs source is fetched as a pinned GitHub tarball and extracted to a
DirectoryAsset whose root contains munge_feature_directory.py, pops.py, and the
example directory. The canonical gene-annotation and control-feature files POPs
ships with live under example/data/utils and are used for real runs as well.

POPs depends only on numpy, scipy, pandas, and scikit-learn, all of which are in
the default pixi environment, so the scripts are run with pixi r python.
"""

from pathlib import Path, PurePath
from typing import Literal, Sequence

import structlog

from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()

# Valid values for munge_feature_directory.py's --nan_policy argument: what to do
# when a feature file is missing genes present in the gene annotation. "raise"
# errors, "ignore" writes nans, "mean"/"zero" impute.
NanPolicy = Literal["raise", "ignore", "mean", "zero"]

# Script names at the root of the extracted POPs source directory.
MUNGE_SCRIPT_NAME = "munge_feature_directory.py"
POPS_SCRIPT_NAME = "pops.py"

# Canonical utility files shipped inside the POPs source tarball. These are
# genome-wide (gene_annot_jun10.txt covers ~18k genes) and are the files the POPs
# authors use for their own runs, so they are appropriate for real analyses.
GENE_ANNOT_RELATIVE_PATH = (
    PurePath("example") / "data" / "utils" / "gene_annot_jun10.txt"
)
CONTROL_FEATURES_RELATIVE_PATH = (
    PurePath("example") / "data" / "utils" / "features_jul17_control.txt"
)

# Stem shared by the chunked feature matrix files produced by the munge step
# (produces <prefix>.mat.<i>.npy, <prefix>.cols.<i>.txt, and <prefix>.rows.txt).
MUNGED_FEATURES_PREFIX_NAME = "pops_features"

# Stem shared by the POPs output files (<stem>.preds, <stem>.coefs, <stem>.marginals).
POPS_OUTPUT_STEM_NAME = "pops_results"
POPS_PREDS_SUFFIX = ".preds"
POPS_COEFS_SUFFIX = ".coefs"
POPS_MARGINALS_SUFFIX = ".marginals"


def count_feature_chunks(munged_features_dir: Path, prefix_name: str) -> int:
    """Return the number of munged feature chunks in a directory, derived from the
    per-chunk column files the munge step writes (<prefix>.cols.<i>.txt). This is
    the value POPs' --num_feature_chunks argument expects."""
    chunk_files = list(munged_features_dir.glob(f"{prefix_name}.cols.*.txt"))
    assert len(chunk_files) > 0, (
        f"No munged feature chunk files ({prefix_name}.cols.*.txt) found in "
        f"{munged_features_dir}"
    )
    return len(chunk_files)


def invoke_pops_script(
    pops_source_dir: Path, script_name: str, args: Sequence[str]
) -> None:
    """Run a POPs script (by name, resolved against the extracted source directory)
    via pixi r python, streaming output and raising on a nonzero exit."""
    script_path = pops_source_dir / script_name
    assert script_path.is_file(), f"POPs script not found: {script_path}"
    cmd = ["pixi", "r", "python", str(script_path), *args]
    logger.debug(f"Running command: {' '.join(cmd)}")
    execute_command(cmd)
