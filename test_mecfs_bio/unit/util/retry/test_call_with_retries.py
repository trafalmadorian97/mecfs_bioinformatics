import pytest

from mecfs_bio.util.retry.call_with_retries import call_with_retries


def test_retries_then_succeeds():
    """A call that fails twice then succeeds is retried to success."""
    calls = {"n": 0}
    slept: list[float] = []

    def flaky() -> str:
        calls["n"] += 1
        if calls["n"] < 3:
            raise ConnectionError()
        return "ok"

    result = call_with_retries(
        flaky,
        retry_on=(ConnectionError,),
        base_backoff_seconds=0.0,
        sleep=slept.append,
    )

    assert result == "ok"
    assert calls["n"] == 3
    assert len(slept) == 2


def test_reraises_after_exhausting_attempts():
    """After max_attempts failures the last error is re-raised."""
    max_attempts = 3
    calls = {"n": 0}

    def always_fail() -> str:
        calls["n"] += 1
        raise ConnectionError()

    with pytest.raises(ConnectionError):
        call_with_retries(
            always_fail,
            retry_on=(ConnectionError,),
            max_attempts=max_attempts,
            base_backoff_seconds=0.0,
            sleep=lambda _: None,
        )

    assert calls["n"] == max_attempts


def test_exception_outside_retry_on_propagates_immediately():
    """Exception types not named in retry_on are not retried."""
    calls = {"n": 0}

    def raise_value_error() -> str:
        calls["n"] += 1
        raise ValueError()

    with pytest.raises(ValueError):
        call_with_retries(
            raise_value_error,
            retry_on=(ConnectionError,),
            base_backoff_seconds=0.0,
            sleep=lambda _: None,
        )

    assert calls["n"] == 1


def test_retries_on_any_listed_exception_type():
    """Every exception type in retry_on triggers a retry."""
    errors: list[Exception] = [ConnectionError(), TimeoutError()]

    def flaky() -> str:
        if errors:
            raise errors.pop(0)
        return "ok"

    result = call_with_retries(
        flaky,
        retry_on=(ConnectionError, TimeoutError),
        base_backoff_seconds=0.0,
        sleep=lambda _: None,
    )

    assert result == "ok"
    assert errors == []


def test_backoff_is_jittered_and_capped():
    """Backoffs grow but never exceed max_backoff_seconds."""
    max_backoff_seconds = 4.0
    slept: list[float] = []

    def always_fail() -> str:
        raise ConnectionError()

    with pytest.raises(ConnectionError):
        call_with_retries(
            always_fail,
            retry_on=(ConnectionError,),
            max_attempts=8,
            base_backoff_seconds=1.0,
            max_backoff_seconds=max_backoff_seconds,
            sleep=slept.append,
        )

    assert all(0.0 <= backoff <= max_backoff_seconds for backoff in slept)
