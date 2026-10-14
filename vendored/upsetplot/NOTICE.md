# Vendored: upsetplot

Vendored from **upsetplot 0.9.0** (https://github.com/jnothman/UpSetPlot), which is the
latest release. Licensed BSD-3-Clause; see `LICENSE` (copyright retained).

## Why vendored

upsetplot 0.9.0 crashes under numpy 2 when rendering count/percentage labels
(`show_counts=True`, our usage). In `plotting.py._label_sizes`:

    margin = 0.01 * abs(np.diff(ax.get_xlim()))

`np.diff` of the 2-element axis-limit tuple returns a length-1 array, so `margin` is a
1-d array and the `ax.text(width + margin, ...)` coordinate is a 1-d array. numpy 2 raises
`TypeError: only 0-dimensional arrays can be converted to Python scalars` when matplotlib
does `float(...)` on it (numpy 1.x only warned). 0.9.0 is the last release, so there is no
upstream version to upgrade to.

## Local modifications

Only `plotting.py._label_sizes` is changed: `.item()` is appended to the three `margin`
assignments (the `right`, `left`, `top` branches) so `margin` is a Python float. Search for
`VENDOR FIX` in `plotting.py`. Nothing else is modified.

Report upstream if the project is ever revived.
