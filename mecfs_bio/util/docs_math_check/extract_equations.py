"""
Extract the TeX equations that MathJax will see from the markdown pages of the docs site.

Each page is converted with Python-Markdown using the markdown_extensions block of mkdocs.yml,
so equations are extracted exactly as pymdownx.arithmatex emits them, after any changes made
by Markdown itself.  Arithmatex (in generic mode) wraps inline math in a span with class
arithmatex and delimiters \\( \\), and display math in a div with delimiters \\[ \\].

The mkdocs-macros Jinja step is not applied, so equations inside files pulled in by Jinja
macros are not extracted.
"""

from pathlib import Path

import markdown
import yaml
from attrs import frozen
from mkdocs.utils.yaml import get_yaml_loader
from selectolax.lexbor import LexborHTMLParser

INLINE_DELIMITERS = (r"\(", r"\)")
DISPLAY_DELIMITERS = (r"\[", r"\]")


@frozen(slots=True)
class MarkdownExtensionConfig:
    """
    The extensions and per-extension options that mkdocs passes to Python-Markdown.
    """

    extensions: tuple[str, ...]
    extension_configs: dict[str, dict[str, object]]


@frozen(slots=True)
class TexEquation:
    """
    One equation as MathJax receives it, without delimiters.

    line is the 1-based line in the markdown source where the equation's TeX starts, or None
    when the TeX does not appear verbatim in the source (for example because Markdown
    changed it, which is itself a likely cause of a rendering error).
    """

    tex: str
    display: bool
    line: int | None


@frozen(slots=True)
class PageEquations:
    md_path: Path
    equations: tuple[TexEquation, ...]


def load_markdown_extension_config(mkdocs_yml: Path) -> MarkdownExtensionConfig:
    """
    Read the markdown_extensions block of mkdocs.yml.  Each entry is either an extension name
    or a single-key mapping from the extension name to its options.
    """
    with open(mkdocs_yml) as infile:
        config = yaml.load(infile, Loader=get_yaml_loader())
    entries = config.get("markdown_extensions", [])
    extensions: list[str] = []
    extension_configs: dict[str, dict[str, object]] = {}
    for entry in entries:
        if isinstance(entry, str):
            extensions.append(entry)
            continue
        assert isinstance(entry, dict) and len(entry) == 1, (
            f"Unexpected markdown_extensions entry in {mkdocs_yml}: {entry!r}"
        )
        ((name, options),) = entry.items()
        extensions.append(name)
        extension_configs[name] = options or {}
    return MarkdownExtensionConfig(
        extensions=tuple(extensions), extension_configs=extension_configs
    )


def extract_page_equations(
    md_path: Path, extension_config: MarkdownExtensionConfig
) -> PageEquations:
    """
    Convert one markdown page to HTML and return its equations in document order.
    """
    source = md_path.read_text(encoding="utf-8")
    html = _markdown_to_html(source, extension_config)
    tex_with_display = [
        _strip_delimiters(node.text(deep=True), display=node.tag == "div")
        for node in LexborHTMLParser(html).css(".arithmatex")
    ]
    lines = _locate_in_source(source, [tex for tex, _ in tex_with_display])
    return PageEquations(
        md_path=md_path,
        equations=tuple(
            TexEquation(tex=tex, display=display, line=line)
            for (tex, display), line in zip(tex_with_display, lines, strict=True)
        ),
    )


def _markdown_to_html(source: str, extension_config: MarkdownExtensionConfig) -> str:
    converter = markdown.Markdown(
        extensions=list(extension_config.extensions),
        extension_configs=extension_config.extension_configs,
    )
    return converter.convert(source)


def _strip_delimiters(text: str, display: bool) -> tuple[str, bool]:
    opening, closing = DISPLAY_DELIMITERS if display else INLINE_DELIMITERS
    assert text.startswith(opening) and text.endswith(closing), (
        f"arithmatex output lacks {opening} {closing} delimiters: {text!r}"
    )
    return text[len(opening) : -len(closing)], display


def _locate_in_source(source: str, tex_strings: list[str]) -> list[int | None]:
    """
    Find the source line of each equation.  Equations are searched for in order, starting
    after the previous match, so repeated identical equations map to successive occurrences.
    """
    lines: list[int | None] = []
    cursor = 0
    for tex in tex_strings:
        needle = tex.strip()
        position = source.find(needle, cursor)
        if position == -1:
            position = source.find(needle)
        if position == -1 or not needle:
            lines.append(None)
            continue
        lines.append(source.count("\n", 0, position) + 1)
        cursor = position + len(needle)
    return lines
