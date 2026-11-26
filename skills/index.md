# Contributor conventions

Durable conventions for this repo, written for human and AI contributors alike. Read the page that
matches what you are about to do; each is short.

Baseline assumed throughout: run commands through pixi (`pixi r <command>`); workflows are built
by the `Task`-based build system (`docs/Codebase_Concepts/Build_System.md`); run
`pixi r invoke green` after any significant change. The repo-root `AGENTS.md` has the shortest
form of these rules.

| Page | Read it when |
|------|--------------|
| [python-style.md](python-style.md) | Writing or refactoring any Python: paths, types and interfaces, fail-fast invariants, function signatures, comments and docstrings, dataframes, subprocesses |
| [build-system-and-tasks.md](build-system-and-tasks.md) | Writing a `Task`, deciding which layer code belongs in, adding external data, or reasoning about caching and the asset store |
| [genomic-data.md](genomic-data.md) | Defining columns, joining gene or protein tables, harmonizing alleles, choosing significance thresholds |
| [testing.md](testing.md) | Adding, changing, or running tests |
| [workflow.md](workflow.md) | Verifying changes, adding dependencies, placing experiments, specs, and notes, committing, reading papers |
| [documentation-site.md](documentation-site.md) | Touching anything under `docs/` or producing figures and tables for it |
