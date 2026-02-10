import hashlib
from pathlib import Path

import structlog

logger = structlog.get_logger()


def verify_hash(downloaded_file: Path, expected_hash: str | None):
    if expected_hash is None:
        return
    logger.debug("Verifying MD5 hash of downloaded file...")
    hash_of_downloaded_file = calc_md5_checksum(downloaded_file)
    if hash_of_downloaded_file == expected_hash:
        logger.debug("Hash verified.")
        return
    head_file(downloaded_file)
    raise AssertionError(
        f"Expected has {hash_of_downloaded_file} of file {downloaded_file} to be equal to {expected_hash}"
    )


def head_file(filename: Path, n=10):
    """Prints the first n lines of a file, like the Unix head command."""
    logger.debug(f"Attempting first {n} lines of file {filename}:")
    try:
        with open(filename) as f:
            for i, line in enumerate(f):
                if i >= n:
                    break
                logger.debug(line.rstrip("\n"))  # rstrip to avoid double newlines
    except FileNotFoundError:
        logger.debug(f"Error: The file {filename} was not found.")
    except UnicodeDecodeError:
        logger.debug(f"file {filename} does not appear to be a text file")


def hash_matches(downloaded_file: Path, expected_hash: str | None) -> bool:
    if expected_hash is None:
        return True
    logger.debug("Verifying MD5 hash of downloaded file...")
    hash_of_downloaded_file = calc_md5_checksum(downloaded_file)
    if hash_of_downloaded_file == expected_hash:
        logger.debug("Hash verified.")
        return True
    return False


def calc_md5_checksum(filepath: Path, chunk_size: int = 8192) -> str:
    hasher = hashlib.md5()
    with open(filepath, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break  # End of file
            hasher.update(chunk)
    return hasher.hexdigest()
