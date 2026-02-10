import subprocess
import sys
from subprocess import SubprocessError

import structlog

logger = structlog.get_logger()


def execute_command(cmd: list[str]) -> str:
    """
    Execute and watch shell command, and return its output
    Source: https://stackoverflow.com/questions/29558074/output-of-subprocess-both-to-pipe-and-directly-to-stdout/77111541#77111541
    """
    output = ""
    str_cmd = " ".join(cmd)
    with subprocess.Popen(
        str_cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        encoding="utf-8",
        errors="replace",
    ) as p:
        while True:
            stdout = p.stdout
            assert stdout is not None
            line = stdout.readline()
            if line != "":
                output += line + "\n"
                print(line)
            elif p.poll() != None:
                break
        sys.stdout.flush()
        if p.returncode == 0:
            return output
        else:
            raise subprocess.CalledProcessError(cmd, p.returncode, output)
