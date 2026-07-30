from pathlib import Path

DOCS_DIRECTORY = Path("docs")

FIGURE_DIRECTORY = DOCS_DIRECTORY / "_figs"

FIGURE_GITHUB_RELEASE_TAG = "figures_v3"

FIGURES_ARCHIVE_TITLE = "figures_v3"

# Manifest mapping figure paths (relative to FIGURE_DIRECTORY) to SHA-256
# content hashes. The manifest is committed to git and is the source of truth
# for which figures are live; the GitHub release stores the blobs themselves.
FIGURE_MANIFEST_PATH = Path("mecfs_bio/figures/figures_manifest.json")
