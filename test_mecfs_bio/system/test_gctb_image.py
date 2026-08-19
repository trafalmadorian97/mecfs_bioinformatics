"""
System test for the self-built GCTB Docker image.

Builds the image from docker/gctb by compiling the pinned GCTB fork from source,
then runs the gctb binary inside it and confirms the GCTB banner appears in the
output. The binary is dynamically linked, so the runtime image carries libstdc++6
and libgomp1 for the C++/OpenMP runtime.
"""

from subprocess import CalledProcessError

from mecfs_bio.build_system.task.sbayesrc import gctb_gwfm_constants as c
from mecfs_bio.util.subproc.run_command import execute_command

_TEST_IMAGE_TAG = "gctb:test"


def test_gctb_image_builds_and_runs() -> None:
    execute_command(
        [
            "docker",
            "build",
            *c.gctb_image_build_args(),
            "-t",
            _TEST_IMAGE_TAG,
            "docker/gctb",
        ]
    )
    try:
        output = execute_command(["docker", "run", "--rm", _TEST_IMAGE_TAG, "gctb"])
    except CalledProcessError as error:
        # gctb with no analysis option prints its banner then exits non-zero.
        output = error.output
    assert "GCTB" in output
