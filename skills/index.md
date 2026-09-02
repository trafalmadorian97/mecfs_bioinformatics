# Contributor conventions index

Durable coding preferences, style guidelines, architecture rules, and hard-won gotchas for this
bioinformatics repo. These were distilled from working sessions on the codebase into a portable
form so every contributor works from the same conventions. Read the page relevant to what you're
about to do; each is short.

Baseline that everything below assumes: commands run through pixi (`pixi r <command>`,
`pixi r python <script>`); bioinformatics workflows are built by the `Task`-based build system
(`docs/Codebase_Concepts/Build_System.md`); after any significant change run `pixi r invoke green`.
See also the repo-root `CLAUDE.md` for the shortest form of these rules.

| Page | What it covers | Reach for it when |
|------|----------------|-------------------|
| [python-style.md](python-style.md) | Paths (`Path`/`PurePath`), `Literal` types, attrs invariants & fail-fast, no bare-tuple returns, named kwargs, free-function helpers, polars/narwhals, column constants, `execute_command`/`call_with_retries` | Writing or refactoring any Python |
| [testing.md](testing.md) | Test at the Task level; cost/benefit for every assertion; no message-text asserts; no monkeypatch/mock (inject instead); no library-presence skipif; Docker system-test pattern | Adding or changing tests |
| [build-system-and-tasks.md](build-system-and-tasks.md) | `.create()` derives meta from deps; domain `Meta` (not `SimpleFileMeta`); `scan_dataframe_asset`; explicit dependency acquisition; the cache doesn't see code changes; asset-store split; the 16GB-runner constraint | Writing a `Task`, or reasoning about rebuilds/caching |
| [workflow.md](workflow.md) | Running & verifying `invoke green` (log-and-grep, testmon caveat); `experiments/claude/` for probes; `experiments/claude/design_specs/` for specs; reading papers with `pdf-env` | Running the toolchain or doing exploratory/design work |
| [documentation-site.md](documentation-site.md) | MkDocs technique tags on Analysis pages (`tags_allowed`, `--strict`); `.nav.yml` bare-string rest-patterns | Editing the `docs/` site |
| [environment-gotchas.md](environment-gotchas.md) | Why the rpy2/pysam/setuptools/numpy pins exist; gwaslab rg-error masking & harmonize OOM; DrvFs `copy2` breakage; SBayesRC-via-Docker | A dependency bump breaks the build, or a tool misbehaves |
