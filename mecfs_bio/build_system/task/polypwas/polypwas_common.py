"""
Shared helpers for the polypwas tasks.

polypwas is installed as a CLI in the pixi env (git dependency).  Its train step
shells out to SBayesRC via Rscript, so we point it at a shim that forwards Rscript
calls into the SBayesRC Docker image (see sbayesrc_utils.write_docker_rscript_shim);
its assoc step needs no R.
"""

from pathlib import Path

from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.util.subproc.run_command import execute_command

_TRAIT_BEARING_METAS = (
    FilteredGWASDataMeta,
    GWASLabSumStatsMeta,
    GWASSummaryDataFileMeta,
    ResultTableMeta,
)


def derive_trait_project(meta: Meta) -> tuple[str, str]:
    """Pull (trait, project) from a dependency's meta, asserting it carries them."""
    assert isinstance(meta, _TRAIT_BEARING_METAS), (
        f"Cannot derive trait/project from meta {meta} of type {type(meta)}; "
        "the polypwas data source must wrap a task whose meta carries trait/project."
    )
    return meta.trait, meta.project


def run_polypwas(args: list[str], *, threads: int) -> str:
    """Run a polypwas CLI subcommand with OMP_NUM_THREADS set."""
    return execute_command(["export", f"OMP_NUM_THREADS={threads};", "polypwas", *args])


def polypwas_setup(rscript_path: Path, *, threads: int) -> str:
    """
    Point polypwas at the given Rscript (a shim forwarding into the SBayesRC Docker
    image) so its SBayesRC training runs in the container.  Idempotent.
    """
    return run_polypwas(["setup", "--rscript", str(rscript_path)], threads=threads)
