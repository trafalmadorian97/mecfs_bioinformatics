import random
import time
from collections.abc import Callable

import structlog

logger = structlog.get_logger()


def call_with_retries[T](
    func: Callable[[], T],
    retry_on: tuple[type[Exception], ...],
    max_attempts: int = 6,
    base_backoff_seconds: float = 2.0,
    max_backoff_seconds: float = 60.0,
    sleep: Callable[[float], None] = time.sleep,
) -> T:
    """
    Call func, retrying with exponential backoff and jitter when it raises one
    of the exception types in retry_on.

    Intended for calls that hit flaky remote services, which intermittently
    fail in ways a later identical call would survive. Exceptions outside
    retry_on propagate immediately, so genuine bugs are not retried. The final
    failure is re-raised so a persistently broken call still fails loudly.

    func takes no arguments: bind its arguments at the call site with
    functools.partial. The sleep function is injectable so tests can run
    without real delays.

    NOTE: for discussion of the benefits of exponential backoff with jitter, see:
    https://sre.google/sre-book/addressing-cascading-failures/
    """
    assert max_attempts >= 1, f"max_attempts must be >= 1, got {max_attempts}"
    assert len(retry_on) >= 1, "retry_on must name at least one exception type"
    last_error: Exception | None = None
    for attempt in range(max_attempts):
        try:
            return func()
        except retry_on as e:
            last_error = e
            if attempt >= max_attempts - 1:
                break
            capped = min(base_backoff_seconds * 2**attempt, max_backoff_seconds)
            backoff = random.uniform(0, capped)
            logger.debug(
                f"Call failed (attempt {attempt + 1}/{max_attempts}, "
                f"{type(e).__name__}: {e}). Backing off {backoff:.1f}s before retry."
            )
            sleep(backoff)
    assert last_error is not None
    raise last_error
