# Python style

Rules for new Python code. The common thread: let the type checker and runtime asserts catch
mistakes, and make call sites self-documenting.

## Paths

- **Use `pathlib.Path`, not `str`**, for every filesystem path: parameters, return types, dict
  keys, fields. Convert to `str` only when serializing (JSON, text, a CLI argument), and convert
  back with `Path(s)` at the deserialization site.
- **Use `PurePath` for a path relative to a base that is not yet known** (e.g. a file inside an
  installed package). It withholds filesystem methods until you concretize it with
  `Path(base) / relpath`.
- Test containment with `child == parent or parent in child.parents`, not string prefixes.

## Types

- **Type a string parameter with a fixed set of values as a `Literal` alias**, defined once near
  the values and reused everywhere: `NanPolicy = Literal["raise", "ignore", "mean", "zero"]`.
- **Avoid `Any` and `Callable[...]`.** Give an injected callable its full signature
  (`Callable[[LdscSumstats, str], LdscFit]`). Instead of `Any` for an object you use through a few
  members, define an interface.
- **For an interface, prefer an `ABC` over a `Protocol`.** Implementations subclass it explicitly,
  so an IDE's go-to-base finds the interface. Use a `Protocol` only for what an ABC cannot do, such
  as describing a third-party class you cannot subclass.
- **Narrow with `assert isinstance(x, T)`, not `typing.cast`.** The assert fails loudly when the
  assumption is wrong; `cast` silently lies. Reserve `cast` for a third-party value whose
  annotation is wrong, with a comment saying so.

## Invariants: fail fast

- **Assert preconditions at the boundary**, with a message naming what was expected and what was
  found. Don't let invalid input flow downstream to fail confusingly or produce wrong results.
- **Make invalid states unrepresentable.** If an attrs field is meaningful only in some
  configurations, type it as optional (`str | None`) and assert the cross-field rule in
  `__attrs_post_init__`.
- **Assert numpy invariants in `__attrs_post_init__`** for attrs classes holding arrays (`ndim`,
  shared length, `dtype.kind`), including private helper classes. A shape comment is not enforced;
  an assert is.

## Functions

- **Don't return bare tuples.** Return a small `@frozen` attrs class with named attributes,
  defined next to the function. Exception: a local helper whose tuple is unpacked a line or two
  later.
- **Pass same-typed arguments by keyword** (`f(source_url=..., dest=...)`). The type checker
  cannot catch two transposed `str` or `Path` arguments.
- **Prefer module-level free functions to helper methods.** A free function takes only what it
  uses, so its signature shows its real dependencies; a method silently receives all of `self`.
- **Inject dependencies instead of hard-wiring them**: make the collaborator (a command runner, a
  downloader, `time.sleep`) a parameter with the production default. This is also how tests
  replace it (see [testing.md](testing.md)).

## Imports, comments, docstrings

- **Put all imports at the top of the module.** Never import inside a function, method, or test.
  The only exception is `tasks.py`.
- **Comments describe the current code, not its history.** Don't write "mirrors the old ...",
  "extracted from ...", or "replaces the previous ...": a reader of the current tree can't see the
  old one. Migration rationale goes in the commit message or PR description.
- **Docstrings follow AGENTS.md**: no backticks around inline code, no RST.
- **In numerical or statistical code, define every symbol.** Each dimension letter used in a
  shape such as `(S, q)` is defined in a nearby docstring (a module-level dimension glossary works
  well). Define acronyms and terse library names in every file that uses them. For code that
  assembles linear-algebra objects, write out the equations and map each parameter to its term.

## Dataframes

- **Use polars for new dataframe code.** Existing pandas code need not be rewritten.
- **Convert a narwhals frame with `.to_polars()`, never `.to_native()`.** `to_native()` returns
  whatever backend narwhals happens to wrap.
- **Use named constants for column names**, not repeated string literals. Each constant lives in
  `mecfs_bio/constants/` in the file for the tool or source that defines the column (e.g.
  `gwaslab_constants.py` for gwaslab's own columns, `ldsc_constants.py` for LDSC's). Add a new file
  for a new source. Only a genuinely local, one-off column may stay a literal.

## Subprocesses and retries

- **Run external commands with `execute_command`** (`mecfs_bio/util/subproc/run_command.py`), not
  raw `subprocess.run` or `Popen`. It logs output and raises on non-zero exit.
- **Retry flaky remote calls with `call_with_retries`** (`mecfs_bio/util/retry/call_with_retries.py`),
  the repo's single backoff loop. Pass an explicit `retry_on` tuple of exception types so that
  genuine bugs are never retried. Don't write another backoff loop.
