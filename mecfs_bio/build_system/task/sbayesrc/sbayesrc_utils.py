"""
Utilities for running SBayesRC via its Docker image.

SBayesRC is distributed as the Docker image zhiliz/sbayesrc, built on Ubuntu 24.04.
We use the image rather than a local R install for two reasons: the package's C++
does not compile against this environment's gcc 15 toolchain, and the image avoids
the broken Ubuntu 22.04 system OpenBLAS 0.3.20 that corrupts SBayesRC's
eigen-decomposition.  The image entrypoint is SBayesRC's own command-line wrapper
(see invoke_sbayesrc), and it also bundles Rscript + SBayesRC, which we expose to
polypwas through write_docker_rscript_shim.

This mirrors the MiXeR Docker helper in mixer_utils.py.
"""

import stat
from pathlib import Path
from typing import Mapping, Sequence

from mecfs_bio.util.subproc.run_command import execute_command

SBAYESRC_IMAGE = "zhiliz/sbayesrc:0.2.6"


def _get_docker_command(extra_mounts: Mapping[Path, Path]) -> list[str]:
    # Run as the host user so container outputs are not root-owned (this Docker is
    # not rootless); set HOME to the writable working dir for R.
    inner_docker_command = (
        "docker run --rm --shm-size=2g "
        "--user $(id -u):$(id -g) -e HOME=/home -v $PWD:/home"
    )
    for key, value in extra_mounts.items():
        inner_docker_command += f" -v {str(key)}:{str(value)}"
    inner_docker_command += " -w /home"
    return ["export", f'DOCKER_RUN="{inner_docker_command}";']


def invoke_sbayesrc(
    args: Sequence[str] | str,
    extra_mounts: Mapping[Path, Path],
) -> None:
    """
    Run the SBayesRC command-line wrapper inside the Docker image.

    The container working directory is the host cwd (mounted at /home), so input
    and output paths passed in args should be relative to the host cwd; reference
    directories are made available through extra_mounts.
    """
    if isinstance(args, str):
        args = [args]
    execute_command(
        _get_docker_command(extra_mounts=extra_mounts)
        + ["export", f'SBAYESRC="$DOCKER_RUN {SBAYESRC_IMAGE}";']
        + ["${SBAYESRC}"]
        + list(args)
    )


def write_docker_rscript_shim(
    extra_mount_dirs: Sequence[Path], shim_path: Path
) -> Path:
    """
    Write an executable shim that forwards Rscript calls into the SBayesRC image.

    polypwas train shells out to SBayesRC via Rscript; pointing polypwas at this
    shim lets it use the container's R + SBayesRC.  The host /tmp and the current
    working directory are bind-mounted at identical paths so the absolute file
    paths polypwas hands to Rscript resolve the same way inside the container; any
    additional directories (e.g. the LD reference) are mounted the same way.
    """
    mounts = ["-v /tmp:/tmp", '-v "$PWD":"$PWD"']
    for mount_dir in extra_mount_dirs:
        resolved = mount_dir.resolve()
        mounts.append(f'-v "{resolved}":"{resolved}"')
    mount_str = " ".join(mounts)
    # The image entrypoint is the SBayesRC wrapper, so override it to run Rscript.
    # Run as the host user (this Docker is not rootless) so SBayesRC's output files
    # are not root-owned, which would otherwise break polypwas's tempdir cleanup.
    shim_path.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        f'exec docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp '
        f'{mount_str} -w "$PWD" '
        f'--entrypoint Rscript {SBAYESRC_IMAGE} "$@"\n'
    )
    shim_path.chmod(shim_path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP)
    return shim_path
