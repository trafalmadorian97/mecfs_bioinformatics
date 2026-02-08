from contextlib import contextmanager
from pathlib import Path


@contextmanager
def log_on_error(filepath: Path):
    """

    Main use case:
    - A system test fails, and I want to spill the current state of the build system info store for debugging
    """
    try:
        print(f"Debug file path is {filepath}")
        yield
    except Exception:
        print(f"\n--- TRACE: Contents of {filepath} ---")
        try:
            with open(filepath) as f:
                print(f.read())
        except FileNotFoundError:
            print(f"Warning: Debug file '{filepath}' not found.")
        print("--- END OF TRACE ---\n")

        raise
