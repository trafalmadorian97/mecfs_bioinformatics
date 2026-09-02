# Documentation site conventions

Docs are a Material for MkDocs site built from `docs/` (config `mkdocs.yml`); serve/build with
`pixi r mkdocs serve` / `build`. CI runs `mkdocs build --strict`, so any broken reference, missing
nav target, or disallowed tag fails the build.

## Analysis pages carry a technique tag

Every page under `docs/Analysis` carries at most one `tags:` entry in its front matter naming its
main analysis technique, and `docs/Analysis/Index_of_Techniques.md` renders them grouped by tag.
The allow list lives in `mkdocs.yml` under the `tags` plugin (`tags_allowed`) — e.g. LDSC, S-LDSC,
CT-LDSC, MAGMA, H-MAGMA, Mendelian randomization, MiXeR, SuSiE, LCV. One `MAGMA` tag covers
gene-level, GTEx, HBA, and gene-set analyses; `H-MAGMA` is separate. Pages that are not about a
technique (overviews, gene lists, Manhattan plots) stay untagged on purpose.

- **A tag outside the allow list — including a typo — fails `mkdocs build --strict`, and therefore
  CI.** A *missing* tag fails nothing, so a new analysis page just silently never appears in the
  index.
- When adding an analysis page, add the tag. When adding a genuinely new technique, extend
  `tags_allowed` first.
- `docs_hooks.py` supplies the disambiguated `<technique> (<Study>, <Trait>)` titles the index
  displays (the study segment comes from path part 2), so page H1s can stay short.

## `.nav.yml` rest-patterns must be bare strings

In mkdocs-awesome-nav `.nav.yml` files, write a rest-pattern as a bare string, never as a
single-key mapping:

```yaml
nav:
  - Index_of_Techniques.md
  - "*"        # correct
  # - glob: "*"  # WRONG — silently misparsed
```

Writing `- glob: "*"` builds but does the wrong thing: pydantic resolves the one-key mapping as a
title→target link, producing a nav entry *titled* "glob" pointing at a page named `*`. The symptom
is `WARNING - A reference to '*' is included in the 'nav' configuration, which is not found ...`,
which fails `--strict`. The dict form is only valid when it also carries pattern options (`sort:`,
`ignore:`) that disambiguate it.
