"""
Constants for interacting with github.com.
"""

GH_OWNER = "trafalmadorian97"
GH_REPO_NAME = f"{GH_OWNER}/mecfs_bioinformatics"

# GitHub Container Registry (ghcr.io) namespace for images we publish, e.g. the
# public gctb image pulled by the remote fine-mapping runner.
GH_CONTAINER_REGISTRY = f"ghcr.io/{GH_OWNER}"
