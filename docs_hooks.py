"""MkDocs build hooks for this site.

Registered via the hooks: key in mkdocs.yml. The only hook here disambiguates
the titles of analysis pages so that search results and the tag index are
readable without changing the (deliberately short) left-hand navigation labels.

Background: analysis pages are organized as Analysis/<Trait>/<Study>/<page>.md
and typically carry a short H1 such as "# LDSC", since the folder tree already
tells the reader which trait and study the page is about. That reads well in the
nav but poorly in search or a tag index, where many pages all show as "LDSC"
with no way to tell them apart.

To fix this without editing every page, on_page_markdown derives the trait and
study from the file path and sets a verbose page title of the form
"LDSC (DecodeME, ME/CFS)". The verbose title is what search, the browser tab,
and the tag index read. It also stashes the original short title in the
nav_title metadata key, which the overrides/partials/nav-item.html template
override prefers when rendering the sidebar, so the nav stays short.
"""

import re
from pathlib import PurePath

# Analysis pages live at Analysis/<Trait>/<Study>/<page>.md. Anything shallower
# (section landing pages) or under a different top-level folder is left alone.
ANALYSIS_ROOT = "Analysis"
MIN_ANALYSIS_PARTS = 4

# Metadata key read by overrides/partials/nav-item.html to keep nav labels short.
NAV_TITLE_KEY = "nav_title"

# Trait folders whose generic underscore-to-space rendering is wrong. Add an
# entry here (folder name -> display string) whenever a trait needs punctuation
# or capitalization that the folder name cannot carry.
TRAIT_DISPLAY = {
    "ME_CFS": "ME/CFS",
}

H1_PATTERN = re.compile(r"^\s*#\s+(.+?)\s*$", re.MULTILINE)


def _humanize(part: str) -> str:
    """Turn a path segment such as Low_Density_Lipoprotein into spaced words."""
    return part.replace("_", " ")


def _trait_display(part: str) -> str:
    """Render a trait folder name, honoring TRAIT_DISPLAY overrides."""
    return TRAIT_DISPLAY.get(part, _humanize(part))


def _trait_and_study(src_uri: str) -> tuple[str, str] | None:
    """Return the humanized (study, trait) for an analysis page, or None.

    Returns None when the page is not an Analysis/<Trait>/<Study>/... page, in
    which case its title should be left untouched.
    """
    parts = PurePath(src_uri).parts
    if len(parts) < MIN_ANALYSIS_PARTS or parts[0] != ANALYSIS_ROOT:
        return None
    trait = _trait_display(parts[1])
    study = _humanize(parts[2])
    return study, trait


def on_page_markdown(markdown: str, *, page, config, files) -> str:
    """Give analysis pages a verbose title for search, tags, and the browser tab.

    The short original title is preserved in page.meta under nav_title so the
    sidebar can stay concise; see the module docstring for the full rationale.
    """
    study_trait = _trait_and_study(page.file.src_uri)
    if study_trait is None:
        return markdown
    study, trait = study_trait

    match = H1_PATTERN.search(markdown)
    short = match.group(1) if match else str(page.title)
    verbose = f"{short} ({study}, {trait})"

    # Verbose title feeds the tag index, the browser tab, and the search
    # fallback; short title feeds the nav via the nav-item template override.
    page.meta["title"] = verbose
    page.meta[NAV_TITLE_KEY] = short

    # Search indexes the rendered H1, so qualify the H1 too when present.
    if match:
        markdown = markdown[: match.start()] + f"# {verbose}" + markdown[match.end() :]
    return markdown
