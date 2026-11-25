# Documentation site

The docs are a Material for MkDocs site built from `docs/` and configured in `mkdocs.yml`. CI
builds with `--strict`, so a broken reference, a missing nav target, or a disallowed tag fails CI.

## Building locally

- **Fetch the figures before building the docs: `pixi r invoke pfig`.** Figures are stored on
  GitHub, not in the repo, and the site won't build without them. Run `pfig` again after pulling
  changes that add figures.
- Then run `pixi r invoke serve-docs` (strict, like CI) or `pixi r mkdocs build --strict`.
  `pixi r invoke sdocs` runs `pfig` and then serves the docs in one step.

## Who edits `docs/`

- **AI agents do not edit prose in `docs/`.** It is human-written. When a code change should
  change what a page says (new expected values, new figures), propose the edit in general terms to a human
  instead of making it. Specs and plans never go in `docs/` (see [workflow.md](workflow.md)).

## Analysis pages carry a technique tag

- Each page under `docs/Analysis` has at most one front-matter `tags:` entry naming its main
  technique. `docs/Analysis/Index_of_Techniques.md` lists pages grouped by tag. Overview pages and
  other pages not about a technique stay untagged.
- Tags must come from `tags_allowed` in `mkdocs.yml`. A tag outside that list (including a typo)
  fails the strict build. A *missing* tag fails nothing, so an untagged page silently drops out of
  the index. For a new technique, extend `tags_allowed` first.

## `.nav.yml` rest-patterns are bare strings

In mkdocs-awesome-nav `.nav.yml` files, write `- "*"`, never `- glob: "*"`. The one-key mapping is
parsed as a link titled "glob" pointing at a page named `*`, which fails `--strict`. The mapping
form is valid only when it also carries options such as `sort:` or `ignore:`.

## Embedding Plotly figures

Use the `plotly_embed` macro (`docs_macros/macros.py`):
`{{ plotly_embed("path/plot.html", id="unique-id") }}`, with an `id` unique within the page.
Don't wrap Plotly iframes in GLightbox; they render black.

## Parquet tables for the docs

Write docs-facing parquet tables with zstd compression and byte-stream-split floats. Use
`write_parquet_table(..., compression="zstd", byte_stream_split_columns=float_column_names(table))`
from `mecfs_bio/build_system/task/dataframe_output.py`, or
`ParquetOutFormat(ParquetWriteOptions(byte_stream_split_floats=True))` where a Task takes an
`OutFormat`. Don't call `pl.DataFrame.write_parquet` directly for these.
