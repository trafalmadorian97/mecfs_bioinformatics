"""
Check that every TeX equation in the docs site parses without error in MathJax.

MathJax does not fail a docs build on bad TeX; it renders the equation as a red error message
in the browser.  This module finds such equations without building the site or opening a
browser:

1. Extract each page's equations through the site's Markdown pipeline (extract_equations).
2. Parse them with the MathJax version the site loads, running under Node.js in the
   mathcheck-env pixi environment (check_tex.cjs), using the site's own MathJax tex
   configuration.
3. Report each equation that MathJax renders as an error, with its file and line.

The MathJax version installed in mathcheck-env must equal the version that mkdocs.yml loads
from the CDN; the check asserts this so that a passing check means the site renders.
"""

import json
import re
import shlex
import tempfile
from pathlib import Path

import yaml
from attrs import frozen
from mkdocs.utils.yaml import get_yaml_loader

from mecfs_bio.util.docs_math_check.extract_equations import (
    PageEquations,
    extract_page_equations,
    load_markdown_extension_config,
)
from mecfs_bio.util.subproc.run_command import execute_command

CHECK_TEX_SCRIPT = Path(__file__).parent / "check_tex.cjs"
MATHCHECK_PIXI_ENV = "mathcheck-env"
CDN_MATHJAX_VERSION_PATTERN = re.compile(r"mathjax@([0-9][^/]*)/")


@frozen(slots=True)
class MathError:
    md_path: Path
    line: int | None
    tex: str
    message: str


@frozen(slots=True)
class MathCheckReport:
    n_pages: int
    n_equations: int
    errors: tuple[MathError, ...]


def check_docs_math(
    docs_dir: Path, mkdocs_yml: Path, site_mathjax_config: Path
) -> MathCheckReport:
    """
    Check every equation in the markdown files under docs_dir.

    site_mathjax_config is the JavaScript file that sets window.MathJax for the site (listed
    in extra_javascript in mkdocs.yml).
    """
    cdn_version = _cdn_mathjax_version(mkdocs_yml)
    extension_config = load_markdown_extension_config(mkdocs_yml)
    pages = [
        extract_page_equations(md_path, extension_config)
        for md_path in sorted(docs_dir.rglob("*.md"))
    ]
    pages = [page for page in pages if page.equations]
    checked = _run_mathjax(pages, site_mathjax_config=site_mathjax_config)
    assert checked.mathjax_version == cdn_version, (
        f"{mkdocs_yml} loads MathJax {cdn_version} but {MATHCHECK_PIXI_ENV} has MathJax "
        f"{checked.mathjax_version}.  Pin both to the same version."
    )
    return MathCheckReport(
        n_pages=len(pages),
        n_equations=sum(len(page.equations) for page in pages),
        errors=_collect_errors(pages, checked.page_errors),
    )


def format_report(report: MathCheckReport) -> str:
    summary = (
        f"Checked {report.n_equations} equations on {report.n_pages} pages: "
        f"{len(report.errors)} error(s)."
    )
    details = [
        f"{error.md_path}:{error.line if error.line is not None else '?'}: "
        f"{error.message}\n    {_one_line(error.tex)}"
        for error in report.errors
    ]
    return "\n".join([*details, summary])


@frozen(slots=True)
class _MathJaxOutput:
    mathjax_version: str
    page_errors: list[list[str | None]]


def _cdn_mathjax_version(mkdocs_yml: Path) -> str:
    with open(mkdocs_yml) as infile:
        config = yaml.load(infile, Loader=get_yaml_loader())
    versions = [
        match.group(1)
        for url in config.get("extra_javascript", [])
        if isinstance(url, str)
        and (match := CDN_MATHJAX_VERSION_PATTERN.search(url)) is not None
    ]
    assert len(versions) == 1, (
        f"Expected exactly one versioned MathJax URL in extra_javascript of {mkdocs_yml}, "
        f"found {versions}"
    )
    return versions[0]


def _run_mathjax(
    pages: list[PageEquations], site_mathjax_config: Path
) -> _MathJaxOutput:
    payload = {
        "pages": [
            {
                "equations": [
                    {"tex": equation.tex, "display": equation.display}
                    for equation in page.equations
                ]
            }
            for page in pages
        ]
    }
    with tempfile.TemporaryDirectory() as tmp:
        input_path = Path(tmp) / "equations.json"
        output_path = Path(tmp) / "results.json"
        input_path.write_text(json.dumps(payload), encoding="utf-8")
        execute_command(
            [
                "pixi",
                "run",
                "-e",
                MATHCHECK_PIXI_ENV,
                "node",
                *(
                    shlex.quote(str(path))
                    for path in (
                        CHECK_TEX_SCRIPT,
                        site_mathjax_config,
                        input_path,
                        output_path,
                    )
                ),
            ]
        )
        result = json.loads(output_path.read_text(encoding="utf-8"))
    return _MathJaxOutput(
        mathjax_version=result["mathjax_version"], page_errors=result["pages"]
    )


def _collect_errors(
    pages: list[PageEquations], page_errors: list[list[str | None]]
) -> tuple[MathError, ...]:
    return tuple(
        MathError(
            md_path=page.md_path,
            line=equation.line,
            tex=equation.tex,
            message=message,
        )
        for page, messages in zip(pages, page_errors, strict=True)
        for equation, message in zip(page.equations, messages, strict=True)
        if message is not None
    )


def _one_line(tex: str) -> str:
    return " ".join(tex.split())
