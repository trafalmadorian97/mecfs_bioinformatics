"""
Utilities for interacting with the MiXeR docker images
Based on documentation here: https://github.com/precimed/mixer
"""

from pathlib import Path
from typing import Mapping, Sequence

from mecfs_bio.util.subproc.run_command import execute_command

MIXER_VERSION = "2.2.1"
SETUP_MIXER_DOCKER = [
    "export",
    f'MIXER_PY="$DOCKER_RUN ghcr.io/precimed/gsa-mixer:{MIXER_VERSION} python /tools/mixer/precimed/mixer.py";',
]

SETUP_MIXER_FIGURES_DOCKER = [
    "export",
    f'MIXER_FIGURES_PY="$DOCKER_RUN ghcr.io/precimed/gsa-mixer:{MIXER_VERSION} python /tools/mixer/precimed/mixer_figures.py";',
]


def _get_docker_command(extra_mounts: Mapping[Path, Path]) -> list[str]:
    inner_docker_command = "docker run --rm --shm-size=2g -v $PWD:/home"
    for key, value in extra_mounts.items():
        inner_docker_command += f" -v {str(key)}:{str(value)}"
    inner_docker_command += " -w /home"
    return ["export", f'DOCKER_RUN="{inner_docker_command}";']


def invoke_mixer(
    args: Sequence[str] | str,
    extra_mounts: Mapping[Path, Path],
):
    if isinstance(args, str):
        args = [args]
    execute_command(
        _get_docker_command(extra_mounts=extra_mounts)
        + SETUP_MIXER_DOCKER
        + ["${MIXER_PY}"]
        + list(args)
    )


def invoke_mixer_figures(
    args: Sequence[str] | str,
    extra_mounts: Mapping[Path, Path],
):
    if isinstance(args, str):
        args = [args]
    execute_command(
        _get_docker_command(extra_mounts=extra_mounts)
        + SETUP_MIXER_FIGURES_DOCKER
        + ["${MIXER_FIGURES_PY}"]
        + list(args)
    )
