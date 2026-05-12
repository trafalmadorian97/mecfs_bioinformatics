"""
Streak-aware state tracking for the daily ``check_links`` workflow.

The daily lychee check fails frequently for spurious reasons (transient host
outages, anti-bot rate limiting, ...). To distinguish a permanently dead link
from a temporarily unreachable one, we record a per-URL "consecutive failure
streak" and only fail the workflow when a URL has been broken for ``threshold``
runs in a row.

The state file is a small JSON map that the GitHub Actions workflow round-trips
through an artifact between runs.

Schema (``link_check_state.json``)::

    {
        "https://example.org/foo": {
            "streak": 3,
            "last_status": "fail",
            "last_run_iso": "2026-05-12"
        },
        ...
    }
"""

from __future__ import annotations

import datetime
import enum
import json
from pathlib import Path
from typing import Any

import attrs


class UrlStatus(enum.StrEnum):
    """Outcome of the most recent lychee check for a URL."""

    OK = "ok"
    FAIL = "fail"


class LycheeReportSchemaError(ValueError):
    """
    Raised when a lychee ``--format json`` report does not match the schema we
    expect.

    We deliberately fail loudly on schema drift rather than silently skipping
    unrecognised entries: silent degradation would let the daily workflow
    report success while in fact checking nothing, defeating the purpose of
    the check. An unexpected schema is an unambiguous signal that lychee was
    upgraded and our parser needs updating.
    """


@attrs.frozen
class UrlState:
    """Tracked state for a single URL."""

    streak: int  # consecutive failing runs; 0 if last run succeeded
    last_status: UrlStatus
    last_run: datetime.date

    def to_dict(self) -> dict[str, int | str]:
        return {
            "streak": self.streak,
            "last_status": str(self.last_status),
            "last_run_iso": self.last_run.isoformat(),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> UrlState:
        return cls(
            streak=int(data["streak"]),
            last_status=UrlStatus(str(data["last_status"])),
            last_run=datetime.date.fromisoformat(str(data["last_run_iso"])),
        )


@attrs.frozen
class UrlStreak:
    """A URL paired with its current consecutive-failure streak."""

    url: str
    streak: int


@attrs.frozen
class StreakResult:
    """Outcome of merging a lychee report into the prior streak state."""

    new_state: dict[str, UrlState]
    persistent_failures: list[UrlStreak]
    new_failures: list[str]
    cleared: list[str]
    all_failing: list[UrlStreak]


@attrs.frozen
class LycheeReport:
    """
    A parsed lychee ``--format json --verbose`` report.

    Construct via ``read_lychee_report``. Each URL field is a deduplicated
    frozenset because lychee checks each unique URL exactly once even when it
    appears in multiple source files (we don't care which files referenced
    it for streak tracking).
    """

    success_urls: frozenset[str]
    redirect_urls: frozenset[str]
    error_urls: frozenset[str]
    timeout_urls: frozenset[str]
    excluded_urls: frozenset[str]
    total: int
    errors: int
    timeouts: int


# (json_key, field_holding_url) for each per-file map lychee emits.
# ``redirect_map`` entries use ``origin`` (the URL as it appears in source)
# plus a separate ``redirects`` list of hops; the other maps use ``url``.
_URL_MAPS: tuple[tuple[str, str], ...] = (
    ("success_map", "url"),
    ("error_map", "url"),
    ("timeout_map", "url"),
    ("excluded_map", "url"),
    ("redirect_map", "origin"),
)

_REQUIRED_TOP_LEVEL_KEYS: tuple[str, ...] = (
    "total",
    "errors",
    "timeouts",
) + tuple(name for name, _ in _URL_MAPS)


def read_lychee_report(path: Path) -> LycheeReport:
    """
    Load and validate a lychee JSON report.

    Raises ``LycheeReportSchemaError`` if the file's structure doesn't match
    what we expect from lychee 0.24 (top-level keys present, per-file maps
    are objects whose values are lists of entry objects, each entry has the
    expected URL field).
    """
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise LycheeReportSchemaError(
            f"Expected JSON object at top level of {path}, got {type(raw).__name__}"
        )

    missing = [key for key in _REQUIRED_TOP_LEVEL_KEYS if key not in raw]
    if missing:
        raise LycheeReportSchemaError(
            f"Lychee report at {path} is missing required keys: {sorted(missing)}. "
            f"This likely indicates a lychee schema change."
        )

    extracted: dict[str, frozenset[str]] = {
        map_name: _extract_urls(raw, map_name, url_field, path)
        for map_name, url_field in _URL_MAPS
    }

    return LycheeReport(
        success_urls=extracted["success_map"],
        error_urls=extracted["error_map"],
        timeout_urls=extracted["timeout_map"],
        excluded_urls=extracted["excluded_map"],
        redirect_urls=extracted["redirect_map"],
        total=_require_int(raw, "total", path),
        errors=_require_int(raw, "errors", path),
        timeouts=_require_int(raw, "timeouts", path),
    )


def _require_int(raw: dict, key: str, path: Path) -> int:
    value = raw[key]
    # bool is an int subclass; reject it to catch shape drift.
    if not isinstance(value, int) or isinstance(value, bool):
        raise LycheeReportSchemaError(
            f"Expected integer for {key!r} in {path}, got {type(value).__name__}"
        )
    return value


def _extract_urls(
    raw: dict, map_name: str, url_field: str, path: Path
) -> frozenset[str]:
    file_map = raw[map_name]
    if not isinstance(file_map, dict):
        raise LycheeReportSchemaError(
            f"Expected {map_name!r} to be a JSON object in {path}, "
            f"got {type(file_map).__name__}"
        )
    urls: set[str] = set()
    for file_path, entries in file_map.items():
        if not isinstance(entries, list):
            raise LycheeReportSchemaError(
                f"Expected {map_name}[{file_path!r}] to be a list in {path}, "
                f"got {type(entries).__name__}"
            )
        for entry in entries:
            if not isinstance(entry, dict):
                raise LycheeReportSchemaError(
                    f"Expected each entry in {map_name}[{file_path!r}] of {path} "
                    f"to be an object, got {type(entry).__name__}"
                )
            if url_field not in entry:
                raise LycheeReportSchemaError(
                    f"Entry in {map_name}[{file_path!r}] of {path} is missing "
                    f"required field {url_field!r}. Entry keys: "
                    f"{sorted(entry.keys())}. This likely indicates a lychee "
                    f"schema change."
                )
            urls.add(str(entry[url_field]))
    return frozenset(urls)


def update_streaks(
    report: LycheeReport,
    previous_state: dict[str, UrlState],
    *,
    run_date: datetime.date,
    threshold: int,
) -> StreakResult:
    """
    Merge a lychee report into the prior streak state.

    Rules:
      * Success (success_map or redirect_map): reset streak to 0.
      * Failure (error_map or timeout_map): increment streak by 1.
      * If a URL appears in both the success and failure sets in the same
        run, the **success wins**. Anti-bot rate limiting often causes
        subsequent checks of an already-checked URL to fail; a single
        successful response is sufficient proof the link is alive.
      * Excluded URLs: carry prior state forward unchanged. They are still
        referenced in the codebase, just intentionally skipped by lychee.
      * URLs not seen at all this run: dropped from state (no longer in
        the codebase).
    """
    if threshold < 1:
        raise ValueError(f"threshold must be >= 1, got {threshold}")

    successful = report.success_urls | report.redirect_urls
    # A single successful check proves the URL is alive even if other
    # attempts in the same run were rate-limited or otherwise failed.
    failed = (report.error_urls | report.timeout_urls) - successful
    excluded = report.excluded_urls

    new_state: dict[str, UrlState] = {}
    new_failures: list[str] = []
    cleared: list[str] = []

    for url in successful:
        prior = previous_state.get(url)
        if prior is not None and prior.streak > 0:
            cleared.append(url)
        new_state[url] = UrlState(streak=0, last_status=UrlStatus.OK, last_run=run_date)

    for url in failed:
        prior = previous_state.get(url)
        prior_streak = prior.streak if prior is not None else 0
        if prior_streak == 0:
            new_failures.append(url)
        new_state[url] = UrlState(
            streak=prior_streak + 1,
            last_status=UrlStatus.FAIL,
            last_run=run_date,
        )

    # Carry prior state forward for URLs that were excluded this run but were
    # previously tracked. They are not "gone from the codebase" so we should
    # not drop them.
    for url in excluded:
        if url in successful or url in failed:
            continue
        prior = previous_state.get(url)
        if prior is not None:
            new_state[url] = prior

    persistent_failures = sorted(
        (
            UrlStreak(url=url, streak=state.streak)
            for url, state in new_state.items()
            if state.streak >= threshold
        ),
        key=lambda s: (-s.streak, s.url),
    )
    all_failing = sorted(
        (
            UrlStreak(url=url, streak=state.streak)
            for url, state in new_state.items()
            if state.last_status == UrlStatus.FAIL
        ),
        key=lambda s: (-s.streak, s.url),
    )

    return StreakResult(
        new_state=new_state,
        persistent_failures=persistent_failures,
        new_failures=sorted(new_failures),
        cleared=sorted(cleared),
        all_failing=all_failing,
    )


def load_state(path: Path) -> dict[str, UrlState]:
    """Load streak state. A missing file is treated as an empty state."""
    if not path.exists():
        return {}
    raw = json.loads(path.read_text(encoding="utf-8"))
    return {url: UrlState.from_dict(entry) for url, entry in raw.items()}


def dump_state(path: Path, state: dict[str, UrlState]) -> None:
    """Write streak state as sorted, indented JSON for diff-friendliness."""
    serialisable = {url: state[url].to_dict() for url in sorted(state)}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(serialisable, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
