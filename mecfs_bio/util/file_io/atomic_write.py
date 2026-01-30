import os
import tempfile
from pathlib import Path
from typing import Any

import yaml


def atomic_write_yaml(
    path: str | Path,
    data: Any,
    *,
    sort_keys: bool = True,
) -> None:
    """
    The idea is that if the process is interrupted while we are writing a yaml file, either the complete old file or the complete new file should remain.  We should not get partially written new file.


    Prior to introducing this function, I encountered an issue where interrupting the build system mid-run could result in a corrupted build info yaml file.


    See https://en.wikipedia.org/wiki/Atomicity_(database_systems)

    This is partially from Chatgpt.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Create temp file in same directory for same-filesystem atomicity.
    # delete=False so we can close it before replace (important on Windows).
    tmp_file = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="\n",
        prefix=path.name + ".",
        suffix=".tmp",
        dir=path.parent,
        delete=False,
    )

    tmp_path = Path(tmp_file.name)

    try:
        with tmp_file as f:
            yaml.dump(
                data,
                f,
                sort_keys=sort_keys,
                default_flow_style=False,
                allow_unicode=True,
            )
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())

        # Atomic replace (overwrites existing).
        tmp_path.replace(path)

        # Best-effort directory fsync for POSIX durability.
        _fsync_dir(path.parent)

    finally:
        # If anything failed before the replace, remove the temp file.
        # If replace succeeded, tmp_path no longer exists at that location.
        try:
            tmp_path.unlink(missing_ok=True)
        except Exception:
            pass


def _fsync_dir(directory: Path) -> None:
    """
    Best-effort fsync of a directory for rename durability (POSIX).
    No-op on platforms/filesystems that don't support it.
    """
    try:
        fd = os.open(directory, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    except OSError:
        return
    try:
        os.fsync(fd)
    except OSError:
        pass
    finally:
        os.close(fd)
