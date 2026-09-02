# Workflow conventions

How to run commands, verify changes, and where different kinds of working files belong.

## Running commands

- All commands run through the pixi package manager: `pixi r <command>`; Python scripts via
  `pixi r python <script>`. Repo-level commands (running tasks, typechecking, docs, lint) live in
  `tasks.py` and are invoked via `pixi r invoke <name>`.
- After any significant change, run `pixi r invoke green` (lint, format, spellcheck, link-check,
  import-check, typecheck, test). Note it runs pytest with `--testmon`, which skips tests unaffected
  by the latest *source* change.

## Verifying `invoke green` — check the summary line, not just the exit code

- **Never inspect `pixi r invoke green` with `| tail -N`.** R-backed tests (TwoSampleMR,
  GenomicSEM, susieR via rpy2) emit hundreds of `R callback write-console:` warning lines, so the
  tail is almost always R noise rather than the pytest verdict. Redirect the whole run to a logfile
  and grep it:
  ```bash
  LOG=<scratchpad>/green.log
  pixi r invoke green > "$LOG" 2>&1; echo "EXIT=$?"
  grep -aE 'passed|failed|error|All checks passed' "$LOG" | tail
  ```
  Use `grep -a` — the log has ANSI escapes / progress-bar control chars that can make grep treat it
  as binary.
- **`EXIT=0` alone does not prove tests ran.** Because green uses `--testmon` (selects tests by
  *source* change), an environment-only change (editing `pyproject.toml`, re-solving with
  `pixi install`) makes testmon report "no tests ran" and green still exits 0. To actually verify an
  environment change, bypass testmon with a full run:
  ```bash
  pixi r python -m pytest --typeguard-packages=mecfs_bio test_mecfs_bio/unit -q
  ```

## Where working files go

- **Exploratory analysis that establishes facts for a design decision → a committed script under
  `experiments/claude/`** (its own subdirectory for a multi-step probe), not throwaway terminal
  one-liners. Tee output to `experiments/claude/logs/`. A committed script is reproducible and
  auditable; add asserts on the expected conclusion so it fails loudly if upstream data changes.
  Follow the existing convention there: module docstring stating what/why (no backticks/RST per
  CLAUDE.md), `from __future__ import annotations`, run via
  `pixi r python experiments/claude/<dir>/<script>.py`.
- **Design / brainstorming spec documents → `experiments/claude/design_specs/`**, named
  `YYYY-MM-DD-<topic>-design.md`. Do NOT put them under `docs/` — that is the published Material for
  MkDocs site (see [documentation-site.md](documentation-site.md)), and working specs would pollute
  the built site.

## Reading papers before deciding

Whenever an important decision hinges on what a paper actually says, **download the paper and read
it** with PDF tools — do not infer its contents from abstracts, memory, or web summaries. (An
indirect argument once produced a confidently wrong conclusion about a cohort definition that the
paper's own text plainly contradicted.)

- The default pixi env has no PDF libraries — use the dedicated **`pdf-env`** environment (feature
  `pdf` in pyproject.toml):
  - `pymupdf` (import `fitz`) — preferred for finding facts; preserves phrases better than
    pdftotext. `pixi r -e pdf-env python -c "import fitz,sys; d=fitz.open(sys.argv[1]); ..." "<pdf>"`.
  - `poppler`'s `pdftotext` CLI for quick dumps: `pixi r -e pdf-env pdftotext "<pdf>" -`.
- Pass PDF paths as an argument and quote paths with spaces; `pixi r` may not resolve a bare
  relative filename.
- Supplementary tables (.xlsx) often hold the definitive cohort/sample-size definitions — read them
  with openpyxl (available in the default env).
- State clearly which conclusions come from the primary text vs. inference, and correct the record
  when a source overturns an earlier claim.
