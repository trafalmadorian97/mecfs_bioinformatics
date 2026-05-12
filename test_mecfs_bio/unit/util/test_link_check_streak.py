"""Unit tests for ``mecfs_bio.util.link_check_streak``."""

import datetime
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from mecfs_bio.util.link_check_streak import (
    LycheeReport,
    LycheeReportSchemaError,
    UrlState,
    UrlStatus,
    UrlStreak,
    dump_state,
    load_state,
    read_lychee_report,
    update_streaks,
)

RUN_DATE = datetime.date(2026, 5, 12)
PREV_DATE = datetime.date(2026, 5, 11)
THRESHOLD = 7


def _make_report(
    *,
    success: list[str] | None = None,
    error: list[str] | None = None,
    timeout: list[str] | None = None,
    redirect: list[str] | None = None,
    excluded: list[str] | None = None,
) -> LycheeReport:
    """Build a ``LycheeReport`` for testing the streak logic in isolation."""
    success = success or []
    error = error or []
    timeout = timeout or []
    redirect = redirect or []
    excluded = excluded or []
    all_urls = set(success) | set(error) | set(timeout) | set(redirect) | set(excluded)
    return LycheeReport(
        success_urls=frozenset(success),
        error_urls=frozenset(error),
        timeout_urls=frozenset(timeout),
        redirect_urls=frozenset(redirect),
        excluded_urls=frozenset(excluded),
        total=len(all_urls),
        errors=len(error),
        timeouts=len(timeout),
    )


def _prior(streak: int, status: UrlStatus = UrlStatus.FAIL) -> UrlState:
    return UrlState(streak=streak, last_status=status, last_run=PREV_DATE)


# ---- update_streaks ---------------------------------------------------------


def test_first_run_with_no_prior_state() -> None:
    """An empty prior state should produce streak=0 on success and streak=1 on failure."""
    report = _make_report(
        success=["https://ok.example/"], error=["https://bad.example/"]
    )

    result = update_streaks(
        report, previous_state={}, run_date=RUN_DATE, threshold=THRESHOLD
    )

    assert result.new_state["https://ok.example/"] == UrlState(
        streak=0, last_status=UrlStatus.OK, last_run=RUN_DATE
    )
    assert result.new_state["https://bad.example/"] == UrlState(
        streak=1, last_status=UrlStatus.FAIL, last_run=RUN_DATE
    )
    assert result.new_failures == ["https://bad.example/"]
    assert result.cleared == []
    assert result.persistent_failures == []


def test_success_resets_streak_to_zero() -> None:
    report = _make_report(success=["https://flaky.example/"])
    prior = {"https://flaky.example/": _prior(streak=4)}

    result = update_streaks(report, prior, run_date=RUN_DATE, threshold=THRESHOLD)

    assert result.new_state["https://flaky.example/"].streak == 0
    assert result.cleared == ["https://flaky.example/"]


def test_failure_increments_streak() -> None:
    report = _make_report(error=["https://broken.example/"])
    prior = {"https://broken.example/": _prior(streak=3)}

    result = update_streaks(report, prior, run_date=RUN_DATE, threshold=THRESHOLD)

    assert result.new_state["https://broken.example/"].streak == 4
    # Not a brand-new failure (was already failing).
    assert result.new_failures == []


def test_url_dropped_when_no_longer_in_codebase() -> None:
    """A URL absent from the lychee report this run is removed from state."""
    report = _make_report(success=["https://still-here.example/"])
    prior = {
        "https://still-here.example/": _prior(streak=2),
        "https://gone.example/": _prior(streak=5),
    }

    result = update_streaks(report, prior, run_date=RUN_DATE, threshold=THRESHOLD)

    assert "https://gone.example/" not in result.new_state
    assert result.new_state["https://still-here.example/"].streak == 0


def test_excluded_url_carries_prior_state_forward() -> None:
    """If a URL is excluded by lychee config this run, its streak should not change."""
    report = _make_report(excluded=["https://excluded.example/"])
    prior_state = _prior(streak=3)
    prior = {"https://excluded.example/": prior_state}

    result = update_streaks(report, prior, run_date=RUN_DATE, threshold=THRESHOLD)

    # Carried forward unchanged - same date, same streak.
    assert result.new_state["https://excluded.example/"] == prior_state


def test_persistent_failure_at_threshold() -> None:
    """A URL whose streak reaches the threshold appears in persistent_failures."""
    report = _make_report(error=["https://dead.example/"])
    prior = {"https://dead.example/": _prior(streak=THRESHOLD - 1)}

    result = update_streaks(report, prior, run_date=RUN_DATE, threshold=THRESHOLD)

    assert result.persistent_failures == [
        UrlStreak(url="https://dead.example/", streak=THRESHOLD)
    ]


def test_persistent_failure_just_below_threshold() -> None:
    """Streak one short of threshold should NOT trigger a persistent failure."""
    report = _make_report(error=["https://flaky.example/"])
    prior = {"https://flaky.example/": _prior(streak=THRESHOLD - 2)}

    result = update_streaks(report, prior, run_date=RUN_DATE, threshold=THRESHOLD)

    assert result.persistent_failures == []
    assert result.new_state["https://flaky.example/"].streak == THRESHOLD - 1


def test_timeout_counts_as_failure() -> None:
    report = _make_report(timeout=["https://slow.example/"])
    prior = {"https://slow.example/": _prior(streak=2)}

    result = update_streaks(report, prior, run_date=RUN_DATE, threshold=THRESHOLD)

    assert result.new_state["https://slow.example/"].streak == 3
    assert result.new_state["https://slow.example/"].last_status == UrlStatus.FAIL


def test_redirect_counts_as_success() -> None:
    report = _make_report(redirect=["https://redirects.example/"])
    prior = {"https://redirects.example/": _prior(streak=2)}

    result = update_streaks(report, prior, run_date=RUN_DATE, threshold=THRESHOLD)

    assert result.new_state["https://redirects.example/"].streak == 0


def test_url_in_both_success_and_failure_is_treated_as_success() -> None:
    """
    Anti-bot defense: if a URL has even one successful check, treat it as a
    success regardless of any failure entries in the same run. A single
    successful response is proof the link is alive; the failures are likely
    rate-limiting noise.
    """
    report = _make_report(
        success=["https://ambiguous.example/"],
        error=["https://ambiguous.example/"],
    )
    prior = {"https://ambiguous.example/": _prior(streak=4)}

    result = update_streaks(report, prior, run_date=RUN_DATE, threshold=THRESHOLD)

    assert result.new_state["https://ambiguous.example/"].last_status == UrlStatus.OK
    assert result.new_state["https://ambiguous.example/"].streak == 0
    assert result.cleared == ["https://ambiguous.example/"]


def test_invalid_threshold_raises() -> None:
    with pytest.raises(ValueError):
        update_streaks(_make_report(), {}, run_date=RUN_DATE, threshold=0)


def test_all_failing_sorted_by_streak_desc() -> None:
    report = _make_report(
        error=["https://a.example/", "https://b.example/", "https://c.example/"]
    )
    prior = {
        "https://a.example/": _prior(streak=0),  # will become 1
        "https://b.example/": _prior(streak=4),  # will become 5
        "https://c.example/": _prior(streak=1),  # will become 2
    }

    result = update_streaks(report, prior, run_date=RUN_DATE, threshold=THRESHOLD)

    assert [s.url for s in result.all_failing] == [
        "https://b.example/",
        "https://c.example/",
        "https://a.example/",
    ]


# ---- state file round-trip --------------------------------------------------


def test_roundtrip_state_file(tmp_path: Path) -> None:
    """dump_state followed by load_state should reproduce the original state."""
    state = {
        "https://a.example/": UrlState(
            streak=0, last_status=UrlStatus.OK, last_run=RUN_DATE
        ),
        "https://b.example/": UrlState(
            streak=5, last_status=UrlStatus.FAIL, last_run=RUN_DATE
        ),
    }
    path = tmp_path / "state.json"

    dump_state(path, state)
    assert path.exists()
    reloaded = load_state(path)

    assert reloaded == state


def test_load_state_missing_file_returns_empty(tmp_path: Path) -> None:
    assert load_state(tmp_path / "does-not-exist.json") == {}


def test_state_file_human_readable_format(tmp_path: Path) -> None:
    """The state file uses the documented JSON shape (status as string, date as ISO)."""
    path = tmp_path / "state.json"
    dump_state(
        path,
        {
            "https://x.example/": UrlState(
                streak=2, last_status=UrlStatus.FAIL, last_run=RUN_DATE
            )
        },
    )

    on_disk = json.loads(path.read_text(encoding="utf-8"))
    assert on_disk == {
        "https://x.example/": {
            "streak": 2,
            "last_status": "fail",
            "last_run_iso": "2026-05-12",
        }
    }


# ---- read_lychee_report -----------------------------------------------------


def _write_report(path: Path, body: dict) -> Path:
    path.write_text(json.dumps(body), encoding="utf-8")
    return path


_MINIMAL_VALID_REPORT = {
    "total": 0,
    "errors": 0,
    "timeouts": 0,
    "success_map": {},
    "error_map": {},
    "timeout_map": {},
    "redirect_map": {},
    "excluded_map": {},
}


def test_read_lychee_report_parses_valid_file(tmp_path: Path) -> None:
    body = {
        **_MINIMAL_VALID_REPORT,
        "total": 2,
        "errors": 1,
        "success_map": {"f.md": [{"url": "https://ok.example/"}]},
        "error_map": {"f.md": [{"url": "https://bad.example/"}]},
    }
    report = read_lychee_report(_write_report(tmp_path / "lychee.json", body))

    assert report.success_urls == frozenset({"https://ok.example/"})
    assert report.error_urls == frozenset({"https://bad.example/"})
    assert report.total == 2
    assert report.errors == 1
    assert report.timeouts == 0


def test_read_lychee_report_handles_redirect_origin_field(tmp_path: Path) -> None:
    """
    Regression: lychee's redirect_map entries key the URL under "origin", not
    "url" (with a separate "redirects" list of hops). The parser must extract
    "origin" for that map specifically.
    """
    body = {
        **_MINIMAL_VALID_REPORT,
        "total": 1,
        "redirect_map": {
            "f.md": [
                {
                    "origin": "https://redirected.example/",
                    "redirects": [{"url": "https://final.example/", "code": 302}],
                }
            ]
        },
    }
    report = read_lychee_report(_write_report(tmp_path / "lychee.json", body))

    assert report.redirect_urls == frozenset({"https://redirected.example/"})


def test_read_lychee_report_dedupes_urls_across_files(tmp_path: Path) -> None:
    body = {
        **_MINIMAL_VALID_REPORT,
        "total": 1,
        "success_map": {
            "a.md": [{"url": "https://shared.example/"}],
            "b.md": [{"url": "https://shared.example/"}],
        },
    }
    report = read_lychee_report(_write_report(tmp_path / "lychee.json", body))

    assert report.success_urls == frozenset({"https://shared.example/"})


def test_read_lychee_report_fails_on_missing_top_level_key(tmp_path: Path) -> None:
    body = {"total": 0, "errors": 0}  # missing timeouts and every *_map
    with pytest.raises(LycheeReportSchemaError, match="missing required keys"):
        read_lychee_report(_write_report(tmp_path / "lychee.json", body))


def test_read_lychee_report_fails_on_entry_missing_url(tmp_path: Path) -> None:
    body = {
        **_MINIMAL_VALID_REPORT,
        "success_map": {"f.md": [{"weird": "shape"}]},
    }
    with pytest.raises(LycheeReportSchemaError, match="missing required field 'url'"):
        read_lychee_report(_write_report(tmp_path / "lychee.json", body))


def test_read_lychee_report_fails_on_redirect_entry_missing_origin(
    tmp_path: Path,
) -> None:
    body = {
        **_MINIMAL_VALID_REPORT,
        "redirect_map": {"f.md": [{"url": "https://x.example/", "redirects": []}]},
    }
    with pytest.raises(
        LycheeReportSchemaError, match="missing required field 'origin'"
    ):
        read_lychee_report(_write_report(tmp_path / "lychee.json", body))


def test_read_lychee_report_fails_on_non_object_map(tmp_path: Path) -> None:
    body = {**_MINIMAL_VALID_REPORT, "success_map": ["not", "an", "object"]}
    with pytest.raises(LycheeReportSchemaError, match="to be a JSON object"):
        read_lychee_report(_write_report(tmp_path / "lychee.json", body))


def test_read_lychee_report_fails_on_non_list_entries(tmp_path: Path) -> None:
    body = {**_MINIMAL_VALID_REPORT, "success_map": {"f.md": "not a list"}}
    with pytest.raises(LycheeReportSchemaError, match="to be a list"):
        read_lychee_report(_write_report(tmp_path / "lychee.json", body))


def test_read_lychee_report_fails_on_non_int_total(tmp_path: Path) -> None:
    body = {**_MINIMAL_VALID_REPORT, "total": "lots"}
    with pytest.raises(LycheeReportSchemaError, match="Expected integer for 'total'"):
        read_lychee_report(_write_report(tmp_path / "lychee.json", body))


def test_read_lychee_report_fails_on_non_object_top_level(tmp_path: Path) -> None:
    path = tmp_path / "lychee.json"
    path.write_text(json.dumps([1, 2, 3]), encoding="utf-8")
    with pytest.raises(LycheeReportSchemaError, match="Expected JSON object"):
        read_lychee_report(path)


# ---- end-to-end schema-drift guard ------------------------------------------


def test_parser_handles_real_lychee_output(tmp_path: Path) -> None:
    """
    Run the installed ``lychee`` binary on a tiny offline fixture and feed the
    resulting JSON through ``read_lychee_report``.

    The other parser tests use hand-crafted JSON, which means a future lychee
    upgrade that changes the report shape could land without any unit test
    failing - we'd only notice in the daily CI run. This test pins the parser
    to whatever lychee version is currently installed in the pixi environment,
    so schema drift surfaces locally on ``pixi r invoke test``.

    Uses ``--offline`` (lychee only resolves ``file://`` and relative paths),
    keeping the test hermetic and fast.
    """
    lychee_binary = shutil.which("lychee")
    # Fail loudly if lychee is missing rather than skip - the streak workflow
    # is meaningless without lychee, so a missing binary is a real problem,
    # not a reason to silently let CI go green.
    assert lychee_binary is not None, (
        "lychee binary not on PATH; run tests via `pixi r ...`"
    )

    # Two-link fixture: one resolvable relative path (success) and one
    # missing file (error). This exercises success_map and error_map together.
    target = tmp_path / "target.md"
    target.write_text("# target\n", encoding="utf-8")
    doc = tmp_path / "doc.md"
    doc.write_text(
        "- working: [t](target.md)\n- broken:  [m](no_such_file.md)\n",
        encoding="utf-8",
    )
    report_path = tmp_path / "lychee.json"

    # lychee exits non-zero when any link fails; that's expected here, so
    # don't set check=True. We only care that the JSON report is produced.
    subprocess.run(
        [
            lychee_binary,
            "--offline",
            "--format",
            "json",
            "--verbose",  # required so success_map gets populated
            "--no-progress",
            "--output",
            str(report_path),
            str(doc),
        ],
        cwd=tmp_path,
        check=False,
        capture_output=True,
    )
    assert report_path.exists(), "lychee did not produce a JSON report"

    # The parser must accept the real report without raising. If lychee
    # changes its schema, read_lychee_report will throw LycheeReportSchemaError
    # here and this test will fail loudly.
    report = read_lychee_report(report_path)

    # And the classification must match what we set up. lychee resolves the
    # relative paths to absolute file:// URLs, so we substring-match.
    assert any("target.md" in url for url in report.success_urls), (
        f"expected target.md in success_urls, got {sorted(report.success_urls)}"
    )
    assert any("no_such_file.md" in url for url in report.error_urls), (
        f"expected no_such_file.md in error_urls, got {sorted(report.error_urls)}"
    )
