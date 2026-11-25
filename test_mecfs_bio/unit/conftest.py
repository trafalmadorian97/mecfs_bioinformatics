import gc

import pytest


def pytest_collection_finish(session: pytest.Session) -> None:
    """
    Freeze every object alive after test collection into gc's permanent generation.

    gwaslab, the runner's post-execute hook, and several tasks call gc.collect() explicitly. A full
    collection walks every tracked object, and after collection the test process holds every
    module the suite imports, so each call costs roughly 0.3s. Frozen objects are skipped by later
    collections, which made the full unit-test suite about 18% faster serially and 27% faster under
    xdist. Objects created during collection live for the whole session anyway, so nothing that
    would otherwise be freed is kept.
    """
    gc.freeze()
