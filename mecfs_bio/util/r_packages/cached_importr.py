"""
Per-process cache for R package bindings created by rpy2's importr.

importr builds a Python wrapper for every function a package exports, on every call.
For a large package such as ggplot2 that takes roughly half a second, even though the
embedded R session has already loaded the package. A process has exactly one embedded
R session, so the bindings for a given package name can be built once and reused.
"""

import functools

from rpy2.robjects.packages import InstalledPackage, InstalledSTPackage, importr


@functools.cache
def cached_importr(name: str) -> InstalledSTPackage | InstalledPackage:
    """
    Return the rpy2 binding for the R package called name, creating it on first use.
    """
    package = importr(name)
    assert isinstance(package, InstalledSTPackage | InstalledPackage), (
        f"Unexpected binding type {type(package)} for R package {name}"
    )
    return package
