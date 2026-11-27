"""
Entry point of the child process started by run_in_child_process: python -m this module with
the directory holding the pickled call (see child_call.run_call).
"""

import sys
from pathlib import Path

from mecfs_bio.util.subproc.child_call import run_call

if __name__ == "__main__":
    run_call(Path(sys.argv[1]))
