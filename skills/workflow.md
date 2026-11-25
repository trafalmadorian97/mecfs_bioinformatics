# Workflow

How to run and verify changes, manage dependencies, and where working files belong.

## Running and verifying

- Run everything through pixi: `pixi r <command>`, `pixi r python <script>`. Repo commands live in
  `tasks.py` and run as `pixi r invoke <name>`.
- **After any significant change, run `pixi r invoke green`** (lint, format, spellcheck,
  link-check, import-check, typecheck, test). Redirect the output to a log file and grep for the
  verdict; don't pipe it through `tail`, because tests that call R through rpy2 flood the end of the
  output with console noise:
  ```bash
  pixi r invoke green > green.log 2>&1; echo "EXIT=$?"
  grep -aE 'passed|failed|error|All checks passed' green.log | tail
  ```
- **`EXIT=0` does not prove tests ran.** `green` runs pytest with `--testmon`, which selects tests
  by *source* changes. After an environment-only change (`pyproject.toml`, `pixi.lock`), it runs
  nothing and still exits 0. To verify such a change, run the full unit suite:
  `pixi r python -m pytest --typeguard-packages=mecfs_bio test_mecfs_bio/unit -q`.

## Dependencies

- **Add a new dependency to the `analysis` feature in `pyproject.toml`**:
  `[tool.pixi.feature.analysis.pypi-dependencies]` if it is on PyPI, otherwise
  `[tool.pixi.feature.analysis.dependencies]` (conda). Put it elsewhere only for a stated reason.
  `pixi add --pypi` defaults to `[project].dependencies`, which is wrong here, so check where the
  entry landed. Match the neighboring pin style (`">=3.8.0,<4"`), run `pixi lock`, and review
  `git diff pixi.lock`.

## Where working files go

- **Agent-generated working files go under `experiments/<agent>/`**, named for the coding agent
  in lowercase (`experiments/claude/`, `experiments/codex/`, `experiments/gemini/`). Human
  contributors use their own folder (e.g. `experiments/tralfamadorian97/`).
- **Exploratory analysis that informs a decision goes in a committed script under
  `experiments/<agent>/`**, not in throwaway one-liners. Use a subdirectory for a multi-step probe.
  Assert the expected conclusion so the script fails loudly if upstream data changes. Tee output
  to `experiments/<agent>/logs/`.
- **Design specs and plans go in `experiments/<agent>/design_specs/`**, named
  `YYYY-MM-DD-<topic>-design.md`. Never put them in `docs/`.
- **Write equation-heavy notes in Typst** (`.typ`, compiled with `pixi r typst compile <file>.typ`),
  not Markdown math. Render to PNG (`--format png`) and inspect the pages; the compiler misses
  layout bugs. Common pitfall: a one-character subscript followed by `(` swallows the parenthesis
  (`K_t(x)` renders as K_{t(x)}); write `K_t (x)`.

## Committing

- **Stage explicit paths**, never `git add -A` or a whole directory.
- **Never commit log files or PDFs.** Commit the source (`.typ`, the script), not its output.

## Reading papers

- **When a decision depends on what a paper says, read the paper**, not an abstract, a summary, or
  your memory of it. Supplementary tables (.xlsx, readable with openpyxl) often hold the definitive
  cohort and sample-size definitions. Say which conclusions come from the text and which are
  inference.
- The default environment has no PDF tools; use `pdf-env`:
  `pixi r -e pdf-env pdftotext "<pdf>" -`, or PyMuPDF (`import fitz`) for more faithful text.
