"""
Pragmatic validation of Docker image references.

The goal is to catch realistic mistakes (empty strings, embedded whitespace,
uppercase or otherwise malformed repository paths, a dangling tag colon) at
graph-construction time, not to reimplement the full distribution/reference
grammar. Standard references such as ghcr.io/trafalmadorian97/gctb:2.5.5,
debian:stable-slim, zhiliz/sbayesrc, and localhost:5000/foo are accepted.
"""

import re

# A single repository path component: lowercase alphanumerics, optionally split
# by a single separator (., _, __, or a run of -) between alphanumeric runs.
_NAME_COMPONENT_RE = re.compile(r"[a-z0-9]+(?:(?:[._]|__|-+)[a-z0-9]+)*")
# A tag: leading word char, then up to 127 more of word chars / . / -.
_TAG_RE = re.compile(r"[a-zA-Z0-9_][a-zA-Z0-9_.-]{0,127}")
# A digest: algorithm (e.g. sha256) then a colon and a hex string.
_DIGEST_RE = re.compile(
    r"[A-Za-z][A-Za-z0-9]*(?:[-_+.][A-Za-z][A-Za-z0-9]*)*:[0-9a-fA-F]{32,}"
)
# A registry host: dotted host labels, optional :port. localhost (no dot) is fine.
_REGISTRY_RE = re.compile(
    r"[a-zA-Z0-9](?:[a-zA-Z0-9-]*[a-zA-Z0-9])?"
    r"(?:\.[a-zA-Z0-9](?:[a-zA-Z0-9-]*[a-zA-Z0-9])?)*"
    r"(?::[0-9]+)?"
)


def _is_valid_name(name: str) -> bool:
    """Validate the registry-and-path portion of a reference (no tag/digest).

    The first slash-separated segment is treated as a registry host when it
    looks like one (contains a dot or colon, or is exactly localhost); otherwise
    it is an implicit-Docker-Hub path component. Every path component must match
    the repository grammar, and there must be at least one.
    """
    parts = name.split("/")
    first = parts[0]
    if len(parts) > 1 and ("." in first or ":" in first or first == "localhost"):
        if not _REGISTRY_RE.fullmatch(first):
            return False
        path_parts = parts[1:]
    else:
        path_parts = parts
    return len(path_parts) > 0 and all(
        _NAME_COMPONENT_RE.fullmatch(part) for part in path_parts
    )


def is_valid_docker_image_reference(ref: str) -> bool:
    """Return whether ref is a plausibly valid Docker image reference.

    Deliberately pragmatic (see module docstring): rejects the empty string,
    embedded whitespace, malformed digests/tags, and uppercase or empty
    repository components, while accepting the standard registry/path:tag@digest
    forms. Not a substitute for the daemon's own parsing.
    """
    if not ref or any(char.isspace() for char in ref):
        return False

    remainder = ref
    if "@" in remainder:
        remainder, _, digest = remainder.partition("@")
        if not _DIGEST_RE.fullmatch(digest):
            return False

    name = remainder
    last_slash = remainder.rfind("/")
    last_colon = remainder.rfind(":")
    # A colon after the last slash introduces a tag (a colon before it is a
    # registry :port, which belongs to the name).
    if last_colon > last_slash:
        name = remainder[:last_colon]
        tag = remainder[last_colon + 1 :]
        if not _TAG_RE.fullmatch(tag):
            return False

    return bool(name) and _is_valid_name(name)
