from pathlib import Path

from mecfs_bio.util.docs_math_check.check_docs_math import check_docs_math

REPO_ROOT = Path(__file__).parents[4]
SITE_MATHJAX_CONFIG = REPO_ROOT / "docs" / "javascripts" / "mathjax.js"

GOOD_PAGE = r"""# Good page

Inline \(\hat{\beta}_j < \sum_{i=1}^{N} x_i\) and dollars $h^2_{\mathrm{SNP}}$.

$$
\begin{aligned}
\chi^2_j &= 1 + \frac{N h^2}{M} \ell_j \label{eq:ldsc}
\end{aligned}
$$

\begin{align}
a &= \cancel{b} \\
c &= \boldsymbol{d}
\end{align}

Price is \$5 and not math.
"""

BAD_PAGE = r"""# Bad page

Fine: \(x^2\)

Missing brace: \(\frac{a}{b\)

Undefined macro: $\alpah + 1$

$$
y = \label{dup} 1
$$

$$
z = \label{dup} 2
$$
"""


def test_check_docs_math_flags_each_kind_of_error(tmp_path: Path) -> None:
    docs = tmp_path / "docs"
    (docs / "sub").mkdir(parents=True)
    (docs / "good.md").write_text(GOOD_PAGE)
    (docs / "sub" / "bad.md").write_text(BAD_PAGE)
    (docs / "no_math.md").write_text("# Nothing here\n")

    report = check_docs_math(
        docs_dir=docs,
        mkdocs_yml=REPO_ROOT / "mkdocs.yml",
        site_mathjax_config=SITE_MATHJAX_CONFIG,
    )

    bad_path = docs / "sub" / "bad.md"
    assert report.n_pages == 2
    assert {(error.md_path, error.line) for error in report.errors} == {
        (bad_path, 5),
        (bad_path, 7),
        (bad_path, 14),
    }
